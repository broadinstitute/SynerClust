#!/usr/bin/env python

import sys
import os
import BlastHandling_bis
# import pickle
import logging
import shutil
import argparse
# import networkx as nx

# def usage():
# 	print """Creates rough clusters based on the BLAST output.
#
# 	WF_MakeRoughClusters.py [node dir] [node name] [min best hit] [child1] [child2]
#
# 	where [node dir] contains data pertinent to [node name], which has children [child1] and [child2].
# 	Tightness of intra-rough-cluster relationships can be regulated by [min best hit].
#
# 	[node dir], dir path
# 	full path to directory
# 	[node name], string
# 	name of node that [node dir] refers to
# 	[child1|2] string
# 	immediate children of [node name]
# 	[min best hit] float
# 	range is (0.0,1.0], 0.0+ is least stringent, 1.0 is most stringent
# 	"""
# 	sys.exit(1)


def main():
	usage = "usage: WF_MakeRoughClusters.py [options]"
	parser = argparse.ArgumentParser(usage)
	parser.add_argument('-dir', dest="node_dir", required=True, help="Path to the \"nodes\" folder. (Required)")
	parser.add_argument('-node', dest="node", required=True, help="Current node name. (Required)")
	parser.add_argument('-m', '--min_best_hit', type=float, dest="min_best_hit", required=True, help="Minimal % of match length for Blastp hits compared to best one.")
	parser.add_argument('-F', '--min_syntenic_fraction', type=float, dest="minSynFrac", required=True, help="Minimal syntenic fraction.")
	parser.add_argument('-diff', '--max_size_diff', type=float, default=3.0, dest="max_size_diff", required=False, help="Maximal ratio difference in size between query and target sequence for Blast.")
	parser.add_argument('children', nargs=2, help="Children nodes. (Required)")
	args = parser.parse_args()

	my_dir = args.node_dir + args.node + "/"
	blast_out = my_dir + "blast.m8"
	n_head = my_dir + "blast_headers.txt"

	if "TREES_FINISHED" in os.listdir(my_dir):
		sys.exit(0)

	FORMAT = "%(asctime)-15s %(levelname)s %(module)s.%(name)s.%(funcName)s at %(lineno)d :\n\t%(message)s\n"
	logger = logging.getLogger()
	logging.basicConfig(filename=my_dir + 'MakeRoughClusters.log', format=FORMAT, filemode='w', level=logging.DEBUG)
	# add a new Handler to print all INFO and above messages to stdout
	ch = logging.StreamHandler(sys.stdout)
	ch.setLevel(logging.INFO)
	logger.addHandler(ch)
	logger.info('Started')

	# # get synteny data
	# pickleSyn = {}
	# for c in args.children:
	# 	synFile = args.node_dir + c + "/synteny_data.pkl"
	# 	pklFile = open(synFile, 'rb')
	# 	pickleSyn[c] = pickle.load(pklFile)
	# 	pklFile.close()
	# logger.debug("Loaded synteny_data")

	# Create rough clusters with trees
	bp = BlastHandling.BlastParse(args.max_size_diff, args.node_dir + args.node + "/")
	logger.debug("Parsed Blast")
	hits = BlastHandling.BlastParse.readBlastM8FromFile(blast_out)
	logger.debug("Read Blast")
	# hits = bp.readBlat()

	bestReciprocalHits = bp.scoreHits(hits, n_head, args.min_best_hit, args.minSynFrac)
	logger.debug("Scored Hits")
	tree_dir = my_dir + "trees"
	if os.path.exists(tree_dir):
		if "old" not in os.listdir(my_dir):
			os.mkdir(os.path.join(my_dir, "old"))
		else:
			shutil.rmtree(os.path.join(my_dir, "old", "*"))
		shutil.move(tree_dir, os.path.join(my_dir, "old"))
	os.mkdir(tree_dir)
	tree_dir = tree_dir + os.sep
	retval = bp.makePutativeClusters(tree_dir, bestReciprocalHits)
	logger.debug("Made Putative Clusters")
	if retval > 0:
		sys.exit(retval)

	trees_done_file = my_dir + "TREES_FINISHED"
	tf = open(trees_done_file, 'w')
	tf.write("Way to go!\n")
	tf.close()

	sys.exit(0)


if __name__ == "__main__":
	main()
# 	if len(sys.argv) == 1:
# 		usage()
# 	else:
# 		main(sys.argv[1:])
