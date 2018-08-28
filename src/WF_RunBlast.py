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


DEVNULL = open(os.devnull, 'w')
QUEUE_ERROR = False
# LOGGER = None


def usage():
	print """Launches BLAST after some brief pre-processing

	WF_RunBlast.py [node directory] [node name] [e value] [number of cores] [child 1] [child 2]
	"""
	sys.exit(1)


class Blaster(multiprocessing.Process):
	def __init__(self, blast_queue, results_queue, LOGGER):
		multiprocessing.Process.__init__(self)
		self.blast_queue = blast_queue
		self.results_queue = results_queue
		self.LOGGER = LOGGER

	def run(self):
		self.LOGGER.info("Running blast")
		while True:
			next_task = self.blast_queue.get()
			if next_task is None:
				self.blast_queue.task_done()
				break
			err = next_task(self.LOGGER)
			if err:
				self.results_queue.put(True)
			else:
				self.results_queue.put(False)
			self.blast_queue.task_done()


class Blast(object):
	def __init__(self, cmd):
		self.cmd = cmd

	def __call__(self, LOGGER):
		try:
			process = subprocess.Popen(self.cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			(out, err) = process.communicate()
			LOGGER.debug("\n\tBlast+ stdout:\n" + out + "\n\tBlast+ stderr:\n" + err)
		except OSError:
			return True
		return False


def main(argv):
	parser = argparse.ArgumentParser()
	parser.add_argument('-dir', dest="node_dir", required=True, help="Path to the \"nodes\" folder. (Required)")
	parser.add_argument('-node', dest="node", required=True, help="Current node name. (Required)")
	parser.add_argument('-t', type=int, dest="numThreads", default=4, help="Number of threads to use. (default = 4)")
	parser.add_argument('-B', '--blast_eval', type=float, dest="blast_eval", default=#BLAST_EVAL_DEFAULT, help="Minimal e-value for Blastp hits. (default = #BLAST_EVAL_DEFAULT)")
	parser.add_argument('-c', '--min-match-coverage', type=float, dest="min_match_coverage", default=0.5, help="Minimal fraction of the query length that the protein Blast+ match needs to cover. (default = 0.5)")
	parser.add_argument('children', nargs=2, help="Children nodes. (Required)")
	args = parser.parse_args()

	# node_dir = argv[0]
	# node = argv[1]
	# children = argv[4:]
	# evalue = argv[2]
	# cores = int(argv[3])
	my_dir = args.node_dir + args.node + "/"

	FORMAT = "%(asctime)-15s %(levelname)s %(module)s.%(name)s.%(funcName)s at %(lineno)d :\n\t%(message)s\n"
	LOGGER = logging.getLogger()
	logging.basicConfig(filename=my_dir + 'RunBlast.log', format=FORMAT, filemode='w', level=logging.DEBUG)
	# add a new Handler to print all INFO and above messages to stdout
	ch = logging.StreamHandler(sys.stdout)
	ch.setLevel(logging.INFO)
	LOGGER.addHandler(ch)
	LOGGER.info('Started')

	# run BLAST if BLAST_FINISHED doesn't exist in the node dir
	if "BLAST_FINISHED" not in os.listdir(my_dir):
		blast_queue = multiprocessing.JoinableQueue()
		results_queue = multiprocessing.Queue()

		blasters = [Blaster(blast_queue, results_queue, LOGGER) for i in xrange(args.numThreads)]
		for b in blasters:
			b.start()

		combine_queue = []
		fastas = []
		heads = []
		my_head = my_dir + "blast_headers.txt"
		# m8 = my_dir + "blast.m8"
		m8s = []
		for c in args.children:
			cpath = args.node_dir + c + "/"
			if "NODE_COMPLETE" not in os.listdir(cpath) and "PICKLES_COMPLETE" not in os.listdir(cpath):
				sys.exit("Error: Missing input for children " + c)
			combined_orphans_header = cpath + c + ".combined_orphans_headers.txt"
			pfile = cpath + c + ".pep"
			c_fasta = my_dir + c + ".blast.fa"
			fastas.append(c_fasta)
			c_head = my_dir + c + ".blast_headers.txt"
			heads.append(c_head)
			blast_out = my_dir + c + ".blast.m8"
			m8s.append(blast_out)
			formatDB_log = my_dir + "formatdb.log"
			os.system("cp " + pfile + " " + c_fasta)
			os.system("grep '>' " + c_fasta + "| cut -f2 -d '>' | grep -v '^combined' > " + c_head)
			if c[0] == 'N':  # if children is a node, so has the combined orphans header
				os.system("cat " + combined_orphans_header + " >> " + c_head)
			os.system("#BLAST_PATHmakeblastdb -in " + c_fasta + " -dbtype prot -logfile " + formatDB_log)

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

			l_map_file = cpath + "/locus_mappings.pkl"
			pklFile = open(l_map_file, 'rb')
			locusMap = pickle.load(pklFile)
			pklFile.close()
			strains = set([])
			strains.add(c)
			for lm in locusMap:
				for lmm in locusMap[lm]:
					strains.add("_".join(lmm.split("_")[:-1]))
			LOGGER.info(strains)

			if len(strains) == 1:
				self_blast_out = my_dir + c + "_self.blast.m8"
				for i in xrange(args.numThreads):
					blast_queue.put(Blast(["#BLAST_PATHblastp", "-outfmt", "6", "-evalue", str(args.blast_eval), "-num_threads", "1", "-db", c_fasta, "-query", c_fasta + "." + "%04d" % (i), "-out", self_blast_out + "." + "%04d" % (i)]))
				combine_queue.append("cat " + self_blast_out + ".* >" + self_blast_out)

		cat_head_cmd = "cat " + heads[0] + " " + heads[1] + " > " + my_head
		LOGGER.info(cat_head_cmd)
		os.system(cat_head_cmd)

		for i in xrange(args.numThreads):
			blast_queue.put(Blast(["#BLAST_PATHblastp", "-outfmt", "6", "-evalue", str(args.blast_eval), "-num_threads", "1", "-db", fastas[0], "-query", fastas[1] + "." + "%04d" % (i), "-out", m8s[1] + "." + "%04d" % (i)]))
			blast_queue.put(Blast(["#BLAST_PATHblastp", "-outfmt", "6", "-evalue", str(args.blast_eval), "-num_threads", "1", "-db", fastas[1], "-query", fastas[0] + "." + "%04d" % (i), "-out", m8s[0] + "." + "%04d" % (i)]))
		combine_queue.append("cat " + m8s[1] + ".* > " + m8s[1])
		combine_queue.append("cat " + m8s[0] + ".* > " + m8s[0])

		for i in xrange(args.numThreads):
			blast_queue.put(None)

		blast_queue.join()

		for i in xrange(args.numThreads):
			if results_queue.get():
				exit("\n\nError while running BLAST+. Please make sure that the path is correct.\n\n")

		for combine in combine_queue:
			os.system(combine)

		if not os.stat(m8s[0]).st_size or not os.stat(m8s[1]).st_size:
			exit("Error, at least one of the Blast+ result files is empty. Please make sure your install is working.")

		# my_m8s = my_dir + "*m8"
		# os.system("cat " + my_m8s + " > " + m8)

		blast_finished_file = my_dir + "BLAST_FINISHED"
		bf = open(blast_finished_file, 'w')
		bf.write("Way to go!\n")
		bf.close()


if __name__ == "__main__":
	if len(sys.argv) == 1:
		usage()
	else:
		main(sys.argv[1:])
