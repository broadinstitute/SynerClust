#!/usr/bin/env python
import sys
import os
import pickle
import logging
import re
from collections import Counter


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

	totalGenes = 0
	# counter = 1
	clusters = {}
	cluster_noOrphan = 0
	for l in locusMap:
		# counter = "_".join(l.split("_")[:-1])
		counter = l[33:]
		clusters[counter] = {'leaves': {}, 'transcripts': []}
		leafKids = locusMap[l]
		for k in leafKids:
			prefix = "_".join(k.split("_")[:-1])
			clusters[counter]['transcripts'].append([l_t[k], prefix])  # [gene_id, genome_tag]
			if prefix not in clusters[counter]['leaves']:
				clusters[counter]['leaves'][prefix] = 0
			clusters[counter]['leaves'][prefix] += 1
		totalGenes += len(clusters[counter]['transcripts'])
		if len(clusters[counter]['transcripts']) > 1:
			cluster_noOrphan += 1
		# counter += 1
	print "total genes: ", totalGenes

	pairs = set([])
	scc_count = 0
	mcc_count = 0
	count = 0
	dir_split = locus_mapping.split("/")
	dir_split.pop()
	dir_split.append("final_clusters.txt")
	cluster_out = "/".join(dir_split)
	dir_split.pop()
	dir_split.append("clust_to_trans.txt")
	cTt_out = "/".join(dir_split)
	cout = open(cluster_out, 'w')
	ct_out = open(cTt_out, 'w')
	nwk_out = open(nodes_path + current_root + "/newicks_full.txt", "w")
	for c in clusters:
		names = []
		distrib_buffer = ""
		for transcript in clusters[c]['transcripts']:
			clusters_out.write("\t".join([c, tagToGenome[transcript[1]], genomeToAnnot[tagToGenome[transcript[1]]], t_n[transcript[0]][0], transcript[0], t_n[transcript[0]][1], t_n[transcript[0]][2] + "\n"]))  # STORE DATA CATALOG INFO: genome name and translation to encoded (locus_tag_file?), annotation file name
			if t_n[transcript[0]][2] is not "None":
				names.append(t_n[transcript[0]][2])
		clusters_out.write("\n")
		distrib_buffer += c + "\t" + Counter(names).most_common(1)[0][0]  # max count for choosing the name to print
		for leaf in leaves:
			if leaf in clusters[c]['leaves']:
				distrib_buffer += "\t" + str(clusters[c]['leaves'][leaf])
			else:
				distrib_buffer += "\t0"
		distrib_buffer += "\n"
		distrib_out.write(distrib_buffer)
		transcripts = clusters[c]['transcripts']
		# cid = "Cluster" + str(c)
		cid = "Cluster" + c
		t_out = cid + " (taxa: " + str(len(clusters[c]['leaves'])) + ", genes: " + str(len(transcripts)) + ")\t" + " ".join(t[0] for t in transcripts) + "\n"
		cout.write(t_out)
		nwk_out.write(cid + ": (" + nwksMap[current_root][current_root + "_" + c][0] + ");\n")
		for t in transcripts:
			ct_out.write(cid + "\t" + t[0] + "\n")
		if len(transcripts) == 1:
			continue
		else:
			count += 1
			for t in transcripts:
				# tNum = int(t)
				for s in transcripts:
					if s[0] == t[0]:
						continue
					# sNum = int(s)
					# tlist = [sNum, tNum]
					tlist = [s[0], t[0]]
					tlist.sort()
					tup = (tlist[0], tlist[1])
					pairs.add(tup)
			if not (len(clusters[c]['leaves']) == num_genomes):
				continue
			else:
				if not (len(transcripts) == num_genomes):
					mcc_count += 1
					continue
				else:
					scc_count += 1
	orphan_count = len(clusters) - cluster_noOrphan
	aux_count = cluster_noOrphan - mcc_count
	print "pairs:", len(pairs)
	print "scc:", scc_count
	print "mcc:", mcc_count
	print "aux:", aux_count
	print "orphans:", orphan_count
	print "non-orphan clusters:", count
	cout.close()
	ct_out.close()
	nwk_out.close()
	distrib_out.close()
	clusters_out.close()

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
