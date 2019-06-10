#!/usr/bin/env python

import argparse
import sys
import os
import cPickle as pickle
import logging
import subprocess
import math
import multiprocessing
import platform
import networkx as nx
from WF_RunBlast import Blast
from WF_RunBlast import Blaster
import BlastHandling
import pdb


DEVNULL = open(os.devnull, 'w')


# ~/SynerClust/bin/Post_clustering_high_similarity.py -dir ../ -node N_0000006_UxeaMVtp_D72DKHI3DrnFA -t 4
def usage():
	print """Launches BLAST after some brief pre-processing
	Post_clustering_high_similarity.py -dir [node directory] -node [node name] -B [e value] -t [number of cores]
	"""
	sys.exit(1)


def main(argv):
	parser = argparse.ArgumentParser()
	parser.add_argument('-dir', dest="node_dir", required=True, help="Path to the \"nodes\" folder. (Required)")
	parser.add_argument('-node', dest="root_node", required=True, help="Current node name. (Required)")
	parser.add_argument('-t', type=int, dest="numThreads", default=4, help="Number of threads to use. (default = 4)")
	# parser.add_argument('-B', '--blast_eval', type=float, dest="blast_eval", default=#BLAST_EVAL_DEFAULT, help="Minimal e-value for Blastp hits. (default = #BLAST_EVAL_DEFAULT)")
	parser.add_argument('-perc', '--min-best-hit', type=float, dest="min_best_hit", default=0.95, help="Minimal percent of identity for the Blast+ match compared to the best. (default = 0.95)")
	parser.add_argument('-id', '--min-perc-identity', type=float, dest="min_percent_identity", default=0.95, help="Minimal percent of identity for the Blast+ match. (default = 0.95)")
	parser.add_argument('-c', '--min-match-coverage', type=float, dest="min_match_coverage", default=0.95, help="Minimal fraction of the query length that the protein Blast+ match needs to cover. (default = 0.95)")
	parser.add_argument('-diff', '--max_size_diff', type=float, default=3.0, dest="max_size_diff", required=False, help="Maximal ratio difference in size between query and target sequence for Blast.")
	args = parser.parse_args()

	c = args.root_node
	my_dir = args.node_dir + args.root_node + "/"

	if "NODE_COMPLETE" not in os.listdir(my_dir) and "PICKLES_COMPLETE" not in os.listdir(my_dir):
		sys.exit("Error: Missing input for children " + c)

	FORMAT = "%(asctime)-15s %(levelname)s %(module)s.%(name)s.%(funcName)s at %(lineno)d :\n\t%(message)s\n"
	logger = logging.getLogger()
	logging.basicConfig(filename=my_dir + 'Post_clustering.log', format=FORMAT, filemode='w', level=logging.DEBUG)
	# add a new Handler to print all INFO and above messages to stdout
	ch = logging.StreamHandler(sys.stdout)
	ch.setLevel(logging.INFO)
	logger.addHandler(ch)
	logger.info('Started')

	c_head = my_dir + c + ".self_blast_post_processing_headers.txt"
	self_blast_out = my_dir + c + "_self_postclustering.blast.m8"

	if "POSTBLAST_COMPLETE" not in os.listdir(my_dir):
		blast_queue = multiprocessing.JoinableQueue()
		results_queue = multiprocessing.Queue()

		blasters = [Blaster(blast_queue, results_queue, logger) for i in xrange(args.numThreads)]
		for b in blasters:
			b.start()

		combined_orphans_header = my_dir + c + ".combined_orphans_headers.txt"
		pfile = my_dir + c + ".pep"
		c_fasta = my_dir + c + ".blast.fa"
		# fastas.append(c_fasta)
		# heads.append(c_head)
		# blast_out = my_dir + c + ".blast.m8"
		# m8s.append(blast_out)
		formatDB_log = my_dir + "formatdb.log"
		os.system("cp " + pfile + " " + c_fasta)
		os.system("grep '>' " + c_fasta + "| cut -f2 -d '>' | grep -v '^combined' > " + c_head)
		if c[0] == 'N':  # if children is a node, so has the combined orphans header
			os.system("cat " + combined_orphans_header + " >> " + c_head)
		# os.system("#BLAST_PATHmakeblastdb -in " + c_fasta + " -dbtype prot -logfile " + formatDB_log)
		os.system("makeblastdb -in " + c_fasta + " -dbtype prot -logfile " + formatDB_log)

		wc_cmd = ["wc", "-l", c_fasta]
		process = subprocess.Popen(wc_cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=DEVNULL)
		output = process.communicate()
		if output[1] is not None:
			sys.exit("error")
		line_count = float(output[0].split()[0])
		chunk_size = int(math.ceil(line_count / args.numThreads))
		if chunk_size % 2 != 0:
			chunk_size += 1
		if platform.system() == "Linux":
			split_cmd = ["split", "-d", "-a", "4", "-l", str(chunk_size), c_fasta, c_fasta + "."]
		elif platform.system() == "Darwin":
			split_cmd = ["gsplit", "-d", "-a", "4", "-l", str(chunk_size), c_fasta, c_fasta + "."]
		else:
			exit("Error: OS is neither Linux nor Mac")
		process = subprocess.Popen(split_cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=DEVNULL)
		output = process.communicate()
		if output[1] is not None:
			sys.exit("error")

		combine_queue = []
		for i in xrange(args.numThreads):
			# blast_queue.put(Blast(["#BLAST_PATHblastp", "-outfmt", "6", "-evalue", "0.00005", "-num_threads", "1", "-db", c_fasta, "-query", c_fasta + "." + "%04d" % (i), "-out", self_blast_out + "." + "%04d" % (i)]))
			blast_queue.put(Blast(["blastp", "-outfmt", "6", "-evalue", "0.00005", "-num_threads", "1", "-db", c_fasta, "-query", c_fasta + "." + "%04d" % (i), "-out", self_blast_out + "." + "%04d" % (i)]))
		combine_queue.append("cat " + self_blast_out + ".* >" + self_blast_out)

		for i in xrange(args.numThreads):
			blast_queue.put(None)

		blast_queue.join()

		for i in xrange(args.numThreads):
			if results_queue.get():
				exit("\n\nError while running BLAST+. Please make sure that the path is correct.\n\n")

		for combine in combine_queue:
			os.system(combine)

		with open(my_dir + "POSTBLAST_COMPLETE", "w") as cr:
			cr.write("Way to go!\n")

	if args.root_node[0] == "N":  # node and not leaf
		with open(args.node_dir + args.root_node + "/combined_orphans_translation_table.pkl") as f:
			translation_table = pickle.load(f)

	bp = BlastHandling.BlastParse(args.max_size_diff, args.node_dir + args.root_node + "/", translation_table)

	bestReciprocalHits = bp.prepareDiGraph(c_head)

	hits = BlastHandling.BlastParse.readBlastM8FromFile(self_blast_out)
	logger.debug("Read Blast")
	bestReciprocalHits = bp.scoreHits(hits, bestReciprocalHits, args.min_best_hit, 1, args.min_percent_identity, args.min_match_coverage)

	subs = list(nx.weakly_connected_component_subgraphs(bestReciprocalHits))

	# pdb.set_trace()

	# [len(subs[i].nodes()) > 1 for i in xrange(len(subs))]

	with open(args.node_dir + args.root_node + "/highly_similar_clusters.txt", "w") as f:
		f.write("//\n")
		for sub in subs:
			if len(sub.nodes()) > 1:
				for i in sub.nodes():
					f.write(i + "\n")
				f.write("//\n")


if __name__ == "__main__":
	if len(sys.argv) == 1:
		usage()
	else:
		main(sys.argv[1:])
