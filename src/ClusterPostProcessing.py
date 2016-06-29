#!/usr/bin/env python
import sys
import os
import pickle
import logging
import re


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
	genomes = os.listdir(genome_path)
	for g in genomes:
		if g.find("locus_tag_file.txt") > -1:
			continue
		dataFile = genome_path + g + "/annotation.txt"
		data = open(dataFile, 'r').readlines()
		for d in data:
			d = d.rstrip()
			line = d.split()
			l_t[line[1]] = line[0]
			l_s[line[1]] = line[7]

	nwksMap = {}
	tmp = os.path.realpath(locus_mapping)
	tmp2 = tmp.rfind('/')
	current_root = tmp[tmp.rfind("/", 0, tmp2) + 1:tmp2]
	nodes_path = genome_path[:genome_path.rfind('/') + 1]
	nodes = os.listdir(nodes_path)
	for n in nodes:
		if n[0] == "N":
			with open(nodes_path + n + "/clusters_newick.pkl", "r") as f:
				nwksMap[n] = pickle.load(f)

	query = re.compile("[a-zA-Z0-9-_]+")
	modified = True
	while(modified):
		modified = False
		for k in nwksMap[current_root].keys():
			for s in query.findall(nwksMap[current_root][k]):
				if len(s) == 39:
					if s[0] == "N":
						if nwksMap[s[:32]][s].count(",") > 0:
							nwksMap[current_root][k] = nwksMap[current_root][k].replace(s, "(" + nwksMap[s[:32]][s] + ")")
						else:
							nwksMap[current_root][k] = nwksMap[current_root][k].replace(s, nwksMap[s[:32]][s])
					elif s[0] == "L":
						nwksMap[current_root][k] = nwksMap[current_root][k].replace(s, l_t[s])
					else:
						continue
					modified = True

	totalGenes = 0
	# counter = 1
	clusters = {}
	cluster_noOrphan = 0
	for l in locusMap:
# 		counter = "_".join(l.split("_")[:-1])
		counter = l[33:]
		clusters[counter] = {'leaves': {}, 'transcripts': []}
		leafKids = locusMap[l]
		for k in leafKids:
			clusters[counter]['transcripts'].append(l_t[k])
			prefix = "_".join(k.split("_")[:-1])
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
	nwk_out = open(nodes_path + "/newicks_full.txt", "w")
	for c in clusters:
		transcripts = clusters[c]['transcripts']
		# cid = "Cluster" + str(c)
		cid = "Cluster" + c
		t_out = cid + " (taxa: " + str(len(clusters[c]['leaves'])) + ", genes: " + str(len(transcripts)) + ")\t" + " ".join(transcripts) + "\n"
		cout.write(t_out)
		nwk_out.write(cid + ": " + nwksMap[current_node][curent_node + "_" + c])
		for t in transcripts:
			ct_out.write(cid + "\t" + t + "\n")
		if len(transcripts) == 1:
			continue
		else:
			count += 1
			for t in transcripts:
				# tNum = int(t)
				for s in transcripts:
					if s == t:
						continue
					# sNum = int(s)
					# tlist = [sNum, tNum]
					tlist = [s, t]
					tlist.sort()
					tup = (tlist[0], tlist[1])
					pairs.add(tup)
			if not (len(clusters[c]['leaves']) == num_genomes):
				continue
			else:
				mcc_count += 1
				if not (len(transcripts) == num_genomes):
					continue
				else:
					scc_count += 1
	orphan_count = len(clusters) - cluster_noOrphan
	aux_count = cluster_noOrphan - mcc_count
	print "pairs:", len(pairs)
	print "scc:", scc_count
	print "mcc:", mcc_count - scc_count
	print "aux:", aux_count
	print "orphans:", orphan_count
	print "non-orphan clusters:", count
	cout.close()
	ct_out.close()
	nwk_out.close()

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
