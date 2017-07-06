#!/usr/bin/env python
import sys
import os
import cPickle as pickle
import logging
import re
from collections import Counter
from WF_FinalizeNode import get_alignment
import argparse


def usage():
	print "USAGE: ClusterPostProcessing.py [path to genomes/ in synergy2 working dir] [path to locus_mapping.pkl of node to be processed] [number of leaf genomes descendent of node to be processed]"


def main():
	usage = "usage: ClusterPostProcessing [options] genomes_directory specific_node_directory number_of_genomes"
	parser = argparse.ArgumentParser(usage)
	parser.add_argument('-a', dest="alignment", type=str, choices=["none", "scc", "all"], default="none", help="Whether to output whole cluster alignments (slow) (Default=False)")
	parser.add_argument('folders', nargs=3, help="genome_directory specific_node_directory/locus_mappings.pkl number_of_genomes (Required)")
	args = parser.parse_args()

	FORMAT = "%(asctime)-15s %(levelname)s %(module)s.%(name)s.%(funcName)s at %(lineno)d :\n\t%(message)s\n"
	logger = logging.getLogger()
	logging.basicConfig(filename='ClusterPostProcessing.log', format=FORMAT, filemode='w', level=logging.DEBUG)
	# add a new Handler to print all INFO and above messages to stdout
	ch = logging.StreamHandler(sys.stdout)
	ch.setLevel(logging.INFO)
	logger.addHandler(ch)
	logger.info('Started')

	# node_path = argv[0]
	locus_mapping = args.folders[1]
	genome_path = args.folders[0]
	num_genomes = int(args.folders[2])

	with open(locus_mapping, 'rb') as pklFile:
		locusMap = pickle.load(pklFile)

	l_t = {}  # locus to transcript
	l_s = {}  # locus to sequence
	t_n = {}  # transcript to IDs: mRNA, alias, name
	tagToGenome = {}  # encoded genome name (tag) to genome name
	genomeToAnnot = {}  # encoded genome name (tag) to annotation file
	# genomes = os.listdir(genome_path)
	nwksMap = {}
	tmp = os.path.realpath(locus_mapping)
	tmp2 = tmp.rfind('/')
	current_root = tmp[tmp.rfind("/", 0, tmp2) + 1:tmp2]
	nodes_path = tmp[:tmp.rfind("/", 0, tmp2) + 1]
	# nodes = os.listdir(nodes_path)

	with open(genome_path + "locus_tag_file.txt", 'r') as f:
		for t in f:
			t.rstrip()
			line = t.split()
			tagToGenome[line[1]] = line[0]

	genomes = []
	nodes = [current_root]
	inter_nodes = [current_root]
	while inter_nodes:
		current_node = inter_nodes.pop()
		children_nodes = tagToGenome[current_node].split(";")
		for child_node in children_nodes:
			if child_node[0] == "N":
				inter_nodes.append(child_node)
			else:  # child_node[0] == "L"
				genomes.append(tagToGenome[child_node])
			nodes.append(child_node)

	for g in genomes:
		if os.path.isdir(genome_path + g):  # shouldn't be needed anymore after changing how the genomes list is created
			with open(genome_path + g + "/annotation.txt", "r") as f:
				d = f.readline().rstrip().split("\t")
				genomeToAnnot[d[0]] = d[1]
				for d in f:
					d = d.rstrip()
					line = d.split("\t")
					l_t[line[1]] = line[0]
					l_s[line[1]] = line[10]
					t_n[line[0]] = [line[7], line[8], line[9]]
		# else: other files in the genome directory (nwk tree with tags and any other user generated file)

	distrib_out = open(nodes_path + current_root + "/cluster_dist_per_genome.txt", "w")
	clusters_out = open(nodes_path + current_root + "/clusters.txt", "w")
	if args.alignment == "all":
		alignment_all_out = open(nodes_path + current_root + "/alignments_all.txt", "w")
		alignment_scc_out = open(nodes_path + current_root + "/alignments_scc.txt", "w")
	elif args.alignment == "scc":
		alignment_scc_out = open(nodes_path + current_root + "/alignments_scc.txt", "w")
	distrib_out.write("#cluster_id\tname")
	leaves = []
	for n in nodes:
		if n[0] == "N":
			with open(nodes_path + n + "/clusters_newick.pkl", "r") as f:
				nwksMap[n] = pickle.load(f)
		else:
			distrib_out.write("\t" + tagToGenome[n])
			leaves.append(n)
	distrib_out.write("\n")

	query = re.compile("[a-zA-Z0-9-_]+")
	modified = True
	while(modified):
		modified = False
		for k in nwksMap[current_root].keys():
			if k != "children":
				for s in query.findall(nwksMap[current_root][k][0]):
					if len(s) == 39:
						if s[0] == "N":
							if nwksMap[s[:32]][s][0].count(",") > 0:  # first line is tree with homology distances, second is with synteny distances
								nwksMap[current_root][k][0] = nwksMap[current_root][k][0].replace(s, "(" + nwksMap[s[:32]][s][0] + ")")
								nwksMap[current_root][k][1] = nwksMap[current_root][k][1].replace(s, "(" + nwksMap[s[:32]][s][1] + ")")
							else:
								nwksMap[current_root][k][0] = nwksMap[current_root][k][0].replace(s, nwksMap[s[:32]][s][0])
								nwksMap[current_root][k][1] = nwksMap[current_root][k][1].replace(s, nwksMap[s[:32]][s][1])
						elif s[0] == "L":
							nwksMap[current_root][k][0] = nwksMap[current_root][k][0].replace(s, tagToGenome[s[:32]] + "_" + l_t[s])
							nwksMap[current_root][k][1] = nwksMap[current_root][k][1].replace(s, tagToGenome[s[:32]] + "_" + l_t[s])
						else:
							continue
						modified = True

	cluster_out = nodes_path + current_root + "/final_clusters.txt"
	cTt_out = nodes_path + current_root + "/clust_to_trans.txt"
	cout = open(cluster_out, 'w')
	ct_out = open(cTt_out, 'w')
	nwk_out = open(nodes_path + current_root + "/newicks_full.txt", "w")
	totalGenes = 0
	scc_count = 0
	mcc_count = 0
	cluster_noOrphan = 0
	for l in locusMap:
		counter = l[33:]
		cid = "Cluster" + counter
		cout_buffer = ""
		ct_out_buffer = ""
		stdin_data = ""
		names = []
		genomes = []
		leafKids = locusMap[l]
		for k in leafKids:
			prefix = "_".join(k.split("_")[:-1])
			cout_buffer += l_t[k] + " "
			clusters_out.write("\t".join([counter, tagToGenome[prefix], genomeToAnnot[tagToGenome[prefix]], t_n[l_t[k]][0], l_t[k], t_n[l_t[k]][1], t_n[l_t[k]][2] + "\n"]))  # STORE DATA CATALOG INFO: genome name and translation to encoded (locus_tag_file?), annotation file name
			if args.alignment == "all" or (args.alignment == "scc" and len(leafKids) == num_genomes):
				stdin_data += ">" + tagToGenome[prefix] + "_" + l_t[k] + "\n" + l_s[k] + "\n"
			if t_n[l_t[k]][2] is not "None":
				names.append(t_n[l_t[k]][2])
			ct_out_buffer += cid + "\t" + l_t[k] + "\n"
			genomes.append(prefix)
		clusters_out.write("\n")
		if args.alignment == "scc" and len(leafKids) == num_genomes:
			ali = get_alignment(stdin_data)
			alignment_scc_out.write(cid + "\n" + ali + "\n")
			if args.alignment == "all":
				alignment_all_out.write(cid + "\n" + ali + "\n")
		elif args.alignment == "all":
			alignment_all_out.write(cid + "\n" + get_alignment(stdin_data) + "\n")

		prefix_count = Counter(genomes)
		distrib_buffer = ""
		if names == []:
			distrib_buffer = counter + "\tNone"  # max count for choosing the name to print when there are no names
		else:
			distrib_buffer = counter + "\t" + Counter(names).most_common(1)[0][0]  # max count for choosing the name to print
		for leaf in leaves:
			distrib_buffer += "\t" + str(prefix_count[leaf])
		distrib_buffer += "\n"
		if len(prefix_count) == num_genomes:
			if prefix_count.most_common(1)[0][1] > 1:
				mcc_count += 1
			else:
				scc_count += 1
		cout.write(cid + " (taxa: " + str(len(prefix_count)) + ", genes: " + str(len(leafKids)) + ")\t" + cout_buffer[:-1] + "\n")  # removing trailing space from cout_buffer
		ct_out.write(ct_out_buffer)
		distrib_out.write(distrib_buffer)
		nwk_out.write(cid + ": (" + nwksMap[current_root][current_root + "_" + counter][0] + ");\n")
		nwk_out.write(cid + ": (" + nwksMap[current_root][current_root + "_" + counter][1] + ");\n")
		totalGenes += len(leafKids)
		if len(leafKids) > 1:
			cluster_noOrphan += 1

	orphan_count = len(locusMap) - cluster_noOrphan
	aux_count = cluster_noOrphan - mcc_count - scc_count
	print "total genes: ", totalGenes
	print "scc:", scc_count
	print "mcc:", mcc_count
	print "aux:", aux_count
	print "orphans:", orphan_count
	print "non-orphan clusters:", cluster_noOrphan  # count

	cout.close()
	ct_out.close()
	nwk_out.close()
	distrib_out.close()
	clusters_out.close()
	if args.alignment == "all":
		alignment_all_out.close()
		alignment_scc_out.close()
	elif args.alignment == "scc":
		alignment_scc_out.close()

	# unreferencing
	l_t = None
	l_s = None
	t_n = None
	tagToGenome = None
	genomeToAnnot = None
	nwksMap = None
	locusMap = None

	print "Now adding inparalogs"
	with open(nodes_path + current_root + "/current_inparalogs.pkl") as f:
		inparalogs = pickle.load(f)

	inpar = {}
	for k in inparalogs:
		inpar[k.split("_")[-1]] = inparalogs[k][0].split("_")[-1]
	inparalogs = None

	changed = True
	while changed:
		changed = False
		for t in inpar:
			if inpar[t] in inpar and inpar[t] != inpar[inpar[t]]:
				inpar[t] = inpar[inpar[t]]
				changed = True

	cluster_to_genes = {}
	cluster_to_genomes = {}
	cluster_to_inparalog_output = {}
	with open(nodes_path + current_root + "/clusters.txt", "r") as f:
		for line in f:
			if line != "\n":
				fields = line.rstrip().split()
				if fields[0] in inpar:
					fields[0] = inpar[fields[0]]
				if fields[0] not in cluster_to_inparalog_output:
					cluster_to_inparalog_output[fields[0]] = []
				cluster_to_inparalog_output[fields[0]].append("\t".join(fields) + "\n")
				if fields[0] not in cluster_to_genes:
					cluster_to_genes[fields[0]] = [fields[4]]
					cluster_to_genomes[fields[0]] = set([fields[1]])
				else:
					cluster_to_genes[fields[0]].append(fields[4])
					cluster_to_genomes[fields[0]].add(fields[1])

	with open(nodes_path + current_root + "/clusters_with_inparalogs.txt", "w") as o:
		for cluster in cluster_to_inparalog_output:
				for line in cluster_to_inparalog_output[cluster]:
					o.write(line)
				o.write("\n")

	scc = 0
	mcc = 0
	orphans = 0
	auxilary = 0
	for cluster in cluster_to_genomes:
		if len(cluster_to_genomes[cluster]) == num_genomes:
			if len(cluster_to_genes[cluster]) == num_genomes:
				scc += 1
			else:
				mcc += 1
		elif len(cluster_to_genes[cluster]) == 1:
			orphans += 1
		else:
			auxilary += 1
	print "with inparalogs scc = " + str(scc)
	print "with inparalogs mcc = " + str(mcc)
	print "with inparalogs orphans = " + str(orphans)
	print "with inparalogs auxilary = " + str(auxilary)
	print "with inparalogs non-orphans = " + str(scc + mcc + auxilary)


if __name__ == "__main__":
	main()
