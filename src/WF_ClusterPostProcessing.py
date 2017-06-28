#!/usr/bin/env python

import sys
import os
import pickle
import logging


def usage():
	print """Checks to make sure that all of the genes associated with the leaf genomes of this node are represented in a transcript.

	WF_ClusterPostProcessing.py [genome dir] [node locus mapping]
	where [genome dir] is "genomes/" contained in your working directory
	where [node locus mapping] is locus_mappings.pkl found in this node's directory
	"""


def main(argv):
	FORMAT = "%(asctime)-15s %(levelname)s %(module)s.%(name)s.%(funcName)s at %(lineno)d :\n\t%(message)s\n"
	logger = logging.getLogger()
	logging.basicConfig(filename="/".join(argv[1].split("/")[:-1]) + '/ClusterPostProcessing.log', format=FORMAT, filemode='w', level=logging.DEBUG)
	# add a new Handler to print all INFO and above messages to stdout
	ch = logging.StreamHandler(sys.stdout)
	ch.setLevel(logging.INFO)
	logger.addHandler(ch)
	logger.info('Started')

	# node_path = argv[0]
	locus_mapping = argv[1]
	genome_path = argv[0]

	locus_dir = locus_mapping.split("locus")[0]
	locus = locus_dir.split("/")[-2]
	locus_files = os.listdir(locus_dir)
	print locus_dir
	print locus
	complete = 0
	locus_pkl = ""
	for lf in locus_files:
		if lf.find("NODE_COMPLETE") > -1 or lf.find("PICKLES_COMPLETE") > -1:
			complete = 1
		elif lf.find(locus + ".pkl") > -1:
			locus_pkl = lf
	if complete == 0:  # give the error for missing children data
		if len(locus_pkl) > 0:
			if "PICKLES_COMPLETE" not in os.listdir(locus_dir):
				exit("Error: Leaf data extraction in " + locus_dir + " not completed.")
		else:
			if "NODE_COMPLETE" not in os.listdir(locus_dir):
				exit("Error: Node " + locus_dir + " computation not completed.")

	pklFile = open(locus_mapping, 'rb')
	locusMap = pickle.load(pklFile)
	pklFile.close()

	l_t = {}  # locus to transcript
	lp_t = {}  # locus prefix to transcript
	genomes = os.listdir(genome_path)
	for g in genomes:
		if not os.path.isdir(genome_path + g):
			continue
		dataFile = genome_path + g + "/annotation.txt"
		data = open(dataFile, 'r').readlines()
		prefix = "_".join((data[1].split("\t")[1]).split("_")[:-1])
		logger.debug("%s splitted to %s" % (data[1].split("\t")[1], prefix))
		lp_t[prefix] = []
		for d in data[1:]:
			line = d.split("\t")
			lp_t[prefix].append(line[0])
			l_t[line[1]] = line[0]

	totalGenes = 0
	counter = 1
	clusters = {}
	tIDs = set([])
	# scc = 0
	clusterDist = {}
	prefixes = set([])
	logger.info(len(locusMap))
	# print len(locusMap)
	for l in locusMap:
		clusters[counter] = {'leaves': {}, 'transcripts': []}
		leafKids = locusMap[l]
		for k in leafKids:
			clusters[counter]['transcripts'].append(l_t[k])
			tIDs.add(l_t[k])
			prefix = "_".join(k.split("_")[:-1])
			logger.debug("%s splitted to %s" % (k, prefix))
			prefixes.add(prefix)
			if prefix not in clusters[counter]['leaves']:
				clusters[counter]['leaves'][prefix] = 0
			clusters[counter]['leaves'][prefix] += 1
		if not len(leafKids) in clusterDist:
			clusterDist[len(leafKids)] = 0
		clusterDist[len(leafKids)] += 1

		totalGenes += len(clusters[counter]['transcripts'])
		# totalGenes+= len(leafKids)
		counter += 1
	print clusterDist
	lp_t_count = 0
	for p in prefixes:
		lp_t_count += len(lp_t[p])
	print totalGenes, len(tIDs), lp_t_count
	if totalGenes - len(tIDs) == 0 and totalGenes - lp_t_count == 0:
		logger.info("ok!")
		# print "ok!"
	else:
		logger.error("error! %s %d %s" % (totalGenes, len(tIDs), lp_t_count))
		# print "error!", totalGenes, len(tIDs), lp_t_count
		sys.exit(1)


if __name__ == "__main__":
	if len(sys.argv) == 1:
		sys.exit(usage())
	else:
		main(sys.argv[1:])
