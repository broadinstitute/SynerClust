#!/usr/bin/env python
import sys
import os
import pickle
import logging
import re
from collections import Counter
from WF_FinalizeNode_threaded import get_alignement


def usage():
	print "USAGE: ClusterPostProcessing.py [path to genomes/ in synergy2 working dir] [path to locus_mapping.pkl of node to be processed] [number of leaf genomes descendent of node to be processed]"


def main(argv):
	FORMAT = "%(asctime)-15s %(levelname)s %(module)s.%(name)s.%(funcName)s at %(lineno)d :\n\t%(message)s\n"
	logger = logging.getLogger()
	logging.basicConfig(filename='ClusterPostProcessing.log', format=FORMAT, filemode='w', level=logging.DEBUG)
	# add a new Handler to print all INFO and above messages to stdout
	ch = logging.StreamHandler(sys.stdout)
	ch.setLevel(logging.INFO)
	logger.addHandler(ch)
	logger.info('Started')

	# node_path = argv[0]
	locus_mapping = argv[1]
	genome_path = argv[0]
	num_genomes = int(argv[2])

	pklFile = open(locus_mapping, 'rb')
	locusMap = pickle.load(pklFile)
	pklFile.close()

	l_t = {}  # locus to transcript
	l_s = {}  # locus to sequence
	t_n = {}  # transcript to IDs: mRNA, alias, name
	tagToGenome = {}  # encoded genome name (tag) to genome name
	genomeToAnnot = {}  # encoded genome name (tag) to annotation file
	genomes = os.listdir(genome_path)
	for g in genomes:
		if g.find("locus_tag_file.txt") > -1:
			with open(genome_path + "locus_tag_file.txt", 'r') as f:
				for t in f:
					t.rstrip()
					line = t.split()
					tagToGenome[line[1]] = line[0]
		else:
			with open(genome_path + g + "/annotation.txt", "r") as f:
				d = f.readline().rstrip().split("\t")
				genomeToAnnot[d[0]] = d[1]
				for d in f:
					# dataFile = genome_path + g + "/annotation.txt"
					# data = open(dataFile, 'r').readlines()
					# for d in data:
					d = d.rstrip()
					line = d.split()
					l_t[line[1]] = line[0]
					l_s[line[1]] = line[10]
					t_n[line[0]] = [line[7], line[8], line[9]]

	nwksMap = {}
	tmp = os.path.realpath(locus_mapping)
	tmp2 = tmp.rfind('/')
	current_root = tmp[tmp.rfind("/", 0, tmp2) + 1:tmp2]
	nodes_path = tmp[:tmp.rfind("/", 0, tmp2) + 1]
	nodes = os.listdir(nodes_path)
	distrib_out = open(nodes_path + current_root + "/cluster_dist_per_genome.txt", "w")
	clusters_out = open(nodes_path + current_root + "/clusters.txt", "w")
	alignement_out = open(nodes_path + current_root + "/alignments.txt", "w")
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
							if nwksMap[s[:32]][s][0].count(",") > 0:
								nwksMap[current_root][k][0] = nwksMap[current_root][k][0].replace(s, "(" + nwksMap[s[:32]][s][0] + ")")
							else:
								nwksMap[current_root][k][0] = nwksMap[current_root][k][0].replace(s, nwksMap[s[:32]][s][0])
						elif s[0] == "L":
							nwksMap[current_root][k][0] = nwksMap[current_root][k][0].replace(s, tagToGenome[s[:32]] + "_" + l_t[s])
						else:
							continue
						modified = True

	cluster_out = nodes_path + current_root + "/final_clusters.txt"
	cTt_out = nodes_path + current_root + "/clust_to_trans.txt"
	cout = open(cluster_out, 'w')
	ct_out = open(cTt_out, 'w')
	nwk_out = open(nodes_path + current_root + "/newicks_full.txt", "w")
	totalGenes = 0
	pairs = set([])
	scc_count = 0
	mcc_count = 0
	cluster_noOrphan = 0
	for l in locusMap:
		# counter = "_".join(l.split("_")[:-1])
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
			stdin_data += ">" + tagToGenome[prefix] + "_" + l_t[k] + "\n" + l_s[k] + "\n"
			if t_n[l_t[k]][2] is not "None":
				names.append(t_n[l_t[k]][2])
			ct_out_buffer += cid + "\t" + l_t[k] + "\n"
			genomes.append(prefix)
		clusters_out.write("\n")
		alignement_out.write(cid + "\n" + get_alignement(stdin_data) + "\n")

		for i in xrange(len(leafKids)):
			for j in xrange(i + 1, len(leafKids)):
				if l_t[leafKids[i]] < l_t[leafKids[j]]:
					pairs.add((l_t[leafKids[i]], l_t[leafKids[j]]))
				else:
					pairs.add((l_t[leafKids[j]], l_t[leafKids[i]]))

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
		totalGenes += len(leafKids)
		if len(leafKids) > 1:
			cluster_noOrphan += 1

	orphan_count = len(locusMap) - cluster_noOrphan
	aux_count = cluster_noOrphan - mcc_count
	print "total genes: ", totalGenes
	print "pairs:", len(pairs)
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
	alignement_out.close()

	ds = locus_mapping.split("/")
	ds.pop()
	mydir = "/".join(ds) + "/"
	pair_pkl = mydir + "tuple_pairs.pkl"

	pdat = open(pair_pkl, 'wb')
	pickle.dump(pairs, pdat)
	pdat.close()


if __name__ == "__main__":
	if len(sys.argv) == 1:
		sys.exit(usage())
	else:
		main(sys.argv[1:])
