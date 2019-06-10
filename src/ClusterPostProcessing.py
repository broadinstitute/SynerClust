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
	# nwksMap = {}
	tmp = os.path.realpath(locus_mapping)
	tmp2 = tmp.rfind('/')
	current_root = tmp[tmp.rfind("/", 0, tmp2) + 1:tmp2]
	nodes_path = tmp[:tmp.rfind("/", 0, tmp2) + 1]
	main_path = tmp[:tmp.rfind("/", 0, len(nodes_path) - 1) + 1]
	# nodes = os.listdir(nodes_path)

	if 'without_inparalogs' not in os.listdir(nodes_path + current_root):
		os.mkdir(nodes_path + current_root + "/without_inparalogs")

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

	distrib_out = open(nodes_path + current_root + "/without_inparalogs/cluster_dist_per_genome_without_inparalogs.txt", "w")
	distrib_inparalogs_out = open(nodes_path + current_root + "/cluster_dist_per_genome.txt", "w")

	if args.alignment == "all":
		alignment_all_out = open(nodes_path + current_root + "/alignments_all.txt", "w")
		alignment_scc_out = open(nodes_path + current_root + "/alignments_scc.txt", "w")
	elif args.alignment == "scc":
		alignment_scc_out = open(nodes_path + current_root + "/alignments_scc.txt", "w")
	distrib_out.write("#cluster_id\tname")
	distrib_inparalogs_out.write("#cluster_id\tname")
	leaves = []
	for n in nodes:
		# if n[0] == "N":  #### EDIT
		# 	with open(nodes_path + n + "/clusters_newick.pkl", "r") as f:
		# 		nwksMap[n] = pickle.load(f)
		# else:
		if n[0] == "L":  # leaf
			distrib_out.write("\t" + tagToGenome[n])
			distrib_inparalogs_out.write("\t" + tagToGenome[n])
			leaves.append(n)
	distrib_out.write("\n")
	distrib_inparalogs_out.write("\n")

	totalGenes = 0
	totalGenes_inparalogs = 0
	scc_count = 0
	scc_inparalogs_count = 0
	mcc_count = 0
	mcc_inparalogs_count = 0
	cluster_noOrphan = 0
	cluster_inparalogs_noOrphan = 0

	with open(nodes_path + current_root + "/current_inparalogs.pkl") as f:
		inparalogs = pickle.load(f)

	inpar = {}
	for k in inparalogs:
		# inpar[k.split("_")[-1]] = inparalogs[k][0].split("_")[-1]
		# if inparalogs[k][0] not in inpar:  # to prevent reciprocal swaps
		inpar[k] = inparalogs[k][0]
	inparalogs = None

	changed = True
	while changed:
		changed = False
		for t in inpar:
			if inpar[t] in inpar and inpar[t] != inpar[inpar[t]]:
				inpar[t] = inpar[inpar[t]]
				changed = True

	cluster_to_genes_inparalogs = {}
	for cluster in locusMap:
		cluster_inparalog = cluster
		# if cluster not in inpar:
		# 	cluster_to_genes_inparalogs[cluster] = []
		# else:
		if cluster in inpar:
			cluster_inparalog = inpar[cluster]
		if cluster_inparalog not in cluster_to_genes_inparalogs:
			cluster_to_genes_inparalogs[cluster_inparalog] = []
		for locus in locusMap[cluster]:
			cluster_to_genes_inparalogs[cluster_inparalog].append(locus)

	clusters_out = open(nodes_path + current_root + "/without_inparalogs/clusters_without_inparalogs.txt", "w")
	clusters_inparalogs_out = open(nodes_path + current_root + "/clusters.txt", "w")
	final_clusters_out = open(nodes_path + current_root + "/without_inparalogs/final_clusters_without_inparalogs.txt", 'w')
	final_clusters_inparalogs_out = open(nodes_path + current_root + "/final_clusters.txt", 'w')
	clust_to_trans_out = open(nodes_path + current_root + "/without_inparalogs/clust_to_trans.txt_without_inparalogs", 'w')
	clust_to_trans_inparalogs_out = open(nodes_path + current_root + "/clust_to_trans.txt", 'w')

	for cluster in cluster_to_genes_inparalogs:
		counter = cluster[33:]
		cid = "Cluster" + counter
		stdin_data = ""
		clust_to_trans_buffer = ""
		clust_to_trans_inparalogs_additional_buffer = ""
		final_clusters_buffer = ""
		final_clusters_inparalogs_additional_buffer = ""
		names = []
		names_inparalogs = []
		genomes = []
		genomes_inparalogs = []
		for locus in cluster_to_genes_inparalogs[cluster]:
			prefix = "_".join(locus.split("_")[:-1])
			out_buffer = "\t".join([counter, tagToGenome[prefix], genomeToAnnot[tagToGenome[prefix]], t_n[l_t[locus]][0], l_t[locus], t_n[l_t[locus]][1], t_n[l_t[locus]][2] + "\n"])
			clusters_inparalogs_out.write(out_buffer)
			genomes_inparalogs.append(prefix)
			if t_n[l_t[locus]][2] is not "None":
				names_inparalogs.append(t_n[l_t[locus]][2])
			if locus in locusMap[cluster]:  # not an added inparalog but an 'ortholog'
				clusters_out.write(out_buffer)
				genomes.append(prefix)
				if t_n[l_t[locus]][2] is not "None":
					names.append(t_n[l_t[locus]][2])
				clust_to_trans_buffer += cid + "\t" + l_t[locus] + "\n"
				final_clusters_buffer += l_t[locus] + " "
			else:
				clust_to_trans_inparalogs_additional_buffer += cid + "\t" + l_t[locus] + "\n"
				final_clusters_inparalogs_additional_buffer += l_t[locus] + " "
			if args.alignment == "all" or (args.alignment == "scc" and len(cluster_to_genes_inparalogs[cluster]) == num_genomes):
				stdin_data += ">" + tagToGenome[prefix] + "_" + l_t[locus] + "\n" + l_s[locus] + "\n"
		clusters_out.write("\n")
		clusters_inparalogs_out.write("\n")

		if args.alignment == "scc" and len(cluster_to_genes_inparalogs[cluster]) == num_genomes and len(set(genomes_inparalogs)) == num_genomes:
			ali = get_alignment(stdin_data)
			alignment_scc_out.write(cid + "\n" + ali + "\n")
			if args.alignment == "all":
				alignment_all_out.write(cid + "\n" + ali + "\n")
		elif args.alignment == "all":
			alignment_all_out.write(cid + "\n" + get_alignment(stdin_data) + "\n")

		clust_to_trans_out.write(clust_to_trans_buffer)
		clust_to_trans_inparalogs_out.write(clust_to_trans_buffer + clust_to_trans_inparalogs_additional_buffer)

		prefix_count = Counter(genomes)
		prefix_count_inparalogs = Counter(genomes_inparalogs)

		final_clusters_out.write(cid + " (taxa: " + str(len(prefix_count)) + ", genes: " + str(len(locusMap[cluster])) + ")\t" + final_clusters_buffer[:-1] + "\n")
		final_clusters_inparalogs_out.write(cid + " (taxa: " + str(len(prefix_count_inparalogs)) + ", genes: " + str(len(cluster_to_genes_inparalogs[cluster])) + ")\t" + final_clusters_buffer + final_clusters_inparalogs_additional_buffer[:-1] + "\n")

		distrib_buffer = ""
		if names == []:
			distrib_buffer = counter + "\tNone"  # max count for choosing the name to print when there are no names
		else:
			distrib_buffer = counter + "\t" + Counter(names).most_common(1)[0][0]  # max count for choosing the name to print
		for leaf in leaves:
			distrib_buffer += "\t" + str(prefix_count[leaf])
		distrib_out.write(distrib_buffer + "\n")

		distrib_buffer = ""
		if names_inparalogs == []:
			distrib_buffer = counter + "\tNone"  # name to print when there are no names
		else:
			distrib_buffer = counter + "\t" + Counter(names_inparalogs).most_common(1)[0][0]  # max count for choosing the name to print
		for leaf in leaves:
			distrib_buffer += "\t" + str(prefix_count_inparalogs[leaf])
		distrib_inparalogs_out.write(distrib_buffer + "\n")

		if len(prefix_count) == num_genomes:
			if prefix_count.most_common(1)[0][1] > 1:
				mcc_count += 1
			else:
				scc_count += 1

		if len(prefix_count_inparalogs) == num_genomes:
			if prefix_count_inparalogs.most_common(1)[0][1] > 1:
				mcc_inparalogs_count += 1
			else:
				scc_inparalogs_count += 1

		if len(locusMap[cluster]) > 1:
			cluster_noOrphan += 1
		totalGenes += len(locusMap[cluster])

		if len(cluster_to_genes_inparalogs[cluster]) > 1:
			cluster_inparalogs_noOrphan += 1
		totalGenes_inparalogs += len(cluster_to_genes_inparalogs[cluster])

	for cluster in locusMap:
		if cluster in cluster_to_genes_inparalogs:
			continue  # already processed in inparalogs loop
		else:
			counter = cluster[33:]
			cid = "Cluster" + counter
			final_clusters_buffer = ""
			clust_to_trans_buffer = ""
			stdin_data = ""
			names = []
			genomes = []
			for locus in locusMap[cluster]:
				prefix = "_".join(locus.split("_")[:-1])
				final_clusters_buffer += l_t[locus] + " "
				clusters_out.write("\t".join([counter, tagToGenome[prefix], genomeToAnnot[tagToGenome[prefix]], t_n[l_t[locus]][0], l_t[locus], t_n[l_t[locus]][1], t_n[l_t[locus]][2] + "\n"]))
				if t_n[l_t[locus]][2] is not "None":
					names.append(t_n[l_t[locus]][2])
				clust_to_trans_buffer += cid + "\t" + l_t[locus] + "\n"
				genomes.append(prefix)
			clusters_out.write("\n")

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

			final_clusters_out.write(cid + " (taxa: " + str(len(prefix_count)) + ", genes: " + str(len(locusMap[cluster])) + ")\t" + final_clusters_buffer[:-1] + "\n")  # removing trailing space from cout_buffer
			clust_to_trans_out.write(clust_to_trans_buffer)
			distrib_out.write(distrib_buffer)

			if len(locusMap[cluster]) > 1:
				cluster_noOrphan += 1
			totalGenes += len(locusMap[cluster])

	orphan_count = len(locusMap) - cluster_noOrphan
	aux_count = cluster_noOrphan - mcc_count - scc_count
	print "total genes: ", totalGenes
	print "scc:", scc_count
	print "mcc:", mcc_count
	print "aux:", aux_count
	print "orphans:", orphan_count
	print "non-orphan clusters:", cluster_noOrphan  # count

	orphan_inparalogs_count = len(cluster_to_genes_inparalogs) - cluster_inparalogs_noOrphan
	aux_inparalogs_count = cluster_inparalogs_noOrphan - mcc_inparalogs_count - scc_inparalogs_count
	print "total genes with inparalogs: ", totalGenes_inparalogs
	print "scc with inparalogs:", scc_inparalogs_count
	print "mcc with inparalogs:", mcc_inparalogs_count
	print "aux with inparalogs:", aux_inparalogs_count
	print "orphans with inparalogs:", orphan_inparalogs_count
	print "non-orphan clusters with inparalogs:", cluster_inparalogs_noOrphan  # count

	clusters_out.close()
	clusters_inparalogs_out.close()
	final_clusters_out.close()
	final_clusters_inparalogs_out.close()
	clust_to_trans_out.close()
	clust_to_trans_inparalogs_out.close()
	distrib_out.close()
	distrib_inparalogs_out.close()

	if args.alignment == "all":
		alignment_all_out.close()
		alignment_scc_out.close()
	elif args.alignment == "scc":
		alignment_scc_out.close()

	if "results" not in os.listdir(main_path):
		os.mkdir(main_path + "results")
	os.system("ln -sf " + nodes_path + current_root + "/clusters.txt " + main_path + "results/clusters.txt")
	os.system("ln -sf " + nodes_path + current_root + "/final_clusters.txt " + main_path + "results/final_clusters.txt")
	os.system("ln -sf " + nodes_path + current_root + "/clust_to_trans.txt " + main_path + "results/clust_to_trans.txt")
	os.system("ln -sf " + nodes_path + current_root + "/cluster_dist_per_genome.txt " + main_path + "results/cluster_dist_per_genome.txt")
	os.system("ln -sf " + nodes_path + current_root + "/cluster_dist_per_genome.txt " + main_path + "results/cluster_dist_per_genome.txt")
	os.system("ln -sf " + nodes_path + current_root + "/cluster_dist_per_genome.txt " + main_path + "results/cluster_dist_per_genome.txt")
	if args.alignment == "all":
		os.system("ln -sf " + nodes_path + current_root + "/alignments_all.txt " + main_path + "results/alignments_all.txt")
		os.system("ln -sf " + nodes_path + current_root + "/alignments_scc.txt " + main_path + "results/alignments_scc.txt")
	elif args.alignment == "scc":
		os.system("ln -sf " + nodes_path + current_root + "/alignments_scc.txt " + main_path + "results/alignments_scc.txt")

	nwksMap = {}
	for n in nodes:
		if n[0] == "N":
			with open(nodes_path + n + "/clusters_newick.pkl", "r") as f:
				nwksMap[n] = pickle.load(f)
	# unreferencing
	# l_t = None
	l_s = None
	t_n = None
	# tagToGenome = None
	genomeToAnnot = None
	# locusMap = None
	query = re.compile("[a-zA-Z0-9-_]+")
	modified = True

	while(modified):
		modified = False
		for k in nwksMap[current_root].keys():
			if k != "children":
				for s in query.findall(nwksMap[current_root][k][0]):
					if len(s) == 39 and s[:32] in nodes:
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

	nwk_out = open(nodes_path + current_root + "/without_inparalogs/newicks_full.txt", "w")
	for l in locusMap:
		counter = l[33:]
		cid = "Cluster" + counter
		nwk_out.write(cid + ": (" + nwksMap[current_root][current_root + "_" + counter][0] + ");\n")
		nwk_out.write(cid + ": (" + nwksMap[current_root][current_root + "_" + counter][1] + ");\n")

	nwk_out.close()


if __name__ == "__main__":
	main()
