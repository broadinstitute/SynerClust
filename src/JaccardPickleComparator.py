#!/usr/bin/env python

import sys, os, pickle

def usage():
	print """
	USAGE: JaccardPickleComparator.py [tuple_pairs.pkl] [a different tuple_pairs.pkl file]
	Prints the intersection, union, and intersection/union or Jaccard Coefficient for two tuple_pairs.pkl
	files.  These files are generated with ClusterPostProcessing.py, and should compare clustering results
	generated for the same set of genomes.
	"""

def main(argv):
	pickle1 = argv[0]
	pickle2 = argv[1]
	
	pklFile = open(pickle1,'rb')
	setA = pickle.load(pklFile)
	pklFile.close()

	pklFile = open(pickle2,'rb')
	setB = pickle.load(pklFile)
	pklFile.close()
	
	#jaccard index is the intersection over the union of the cluster pairs that share an ortholog group in each set
	inter_set = setA & setB
	intersect = float(len(inter_set))
	union = float(len(setA | setB))
	
	index = intersect/union
	print "intersecting:",intersect
	print "all:",union
	print "Jaccard Coefficient:",index

if __name__ == '__main__':
	if len(sys.argv) == 1:
		usage()
		sys.exit(1)
	else:
		main(sys.argv[1:])