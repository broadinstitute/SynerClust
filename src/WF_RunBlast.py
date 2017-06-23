#!/usr/bin/env python
import sys
import os
import pickle
import logging
import subprocess
import math
from multiprocessing import Process, Queue
from Queue import Empty
import platform


DEVNULL = open(os.devnull, 'w')
QUEUE_ERROR = False
LOGGER = None


def usage():
	print """Launches BLAST after some brief pre-processing

	WF_RunBlast.py [node directory] [node name] [e value] [number of cores] [child 1] [child 2]
	"""
	sys.exit(1)


def run_blast(blast_queue):
	print "running blast"
	while True:
		try:
			blast_cmd = blast_queue.get(block=False)
			process = subprocess.Popen(blast_cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			(out, err) = process.communicate()
			LOGGER.debug(out)
			LOGGER.warning(err)
		except Empty:
			break
		except:
			QUEUE_ERROR = True


def main(argv):
	node_dir = argv[0]
	node = argv[1]
	children = argv[4:]
	evalue = argv[2]
	cores = int(argv[3])
	my_dir = node_dir + node + "/"

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
		blast_queue = Queue(0)
		combine_queue = []
		fastas = []
		heads = []
		my_head = my_dir + "blast_headers.txt"
		m8 = my_dir + "blast.m8"
		m8s = []
		processes = []
		for c in children:
			cpath = node_dir + c + "/"
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
			# os.system("#BLAST_PATHformatdb -i " + c_fasta + " -l " + formatDB_log)
			os.system("#BLAST_PATHmakeblastdb -in " + c_fasta + " -dbtype prot -logfile " + formatDB_log)

			wc_cmd = ["wc", "-l", c_fasta]
			process = subprocess.Popen(wc_cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=DEVNULL)
			output = process.communicate()
			if output[1] is not None:
				sys.exit("error")
			line_count = float(output[0].split()[0])
			chunk_size = int(math.ceil(line_count / cores))
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
				for i in xrange(cores):
					blast_queue.put(["#BLAST_PATHblastp", "-outfmt", "6", "-evalue", evalue, "-qcov_hsp_perc", "50", "-num_threads", "1", "-db", c_fasta, "-query", c_fasta + "." + "%04d" % (i), "-out", self_blast_out + "." + "%04d" % (i)])
				combine_queue.append("cat " + self_blast_out + ".* >" + self_blast_out)

		cat_head_cmd = "cat " + heads[0] + " " + heads[1] + " > " + my_head
		LOGGER.info(cat_head_cmd)
		os.system(cat_head_cmd)

		for i in xrange(cores):
			blast_queue.put(["#BLAST_PATHblastp", "-outfmt", "6", "-evalue", evalue, "-qcov_hsp_perc", "50", "-num_threads", "1", "-db", fastas[0], "-query", fastas[1] + "." + "%04d" % (i), "-out", m8s[1] + "." + "%04d" % (i)])
			blast_queue.put(["#BLAST_PATHblastp", "-outfmt", "6", "-evalue", evalue, "-qcov_hsp_perc", "50", "-num_threads", "1", "-db", fastas[1], "-query", fastas[0] + "." + "%04d" % (i), "-out", m8s[0] + "." + "%04d" % (i)])
		combine_queue.append("cat " + m8s[1] + ".* > " + m8s[1])
		combine_queue.append("cat " + m8s[0] + ".* > " + m8s[0])

		processes = [Process(target=run_blast, args=(blast_queue,)) for i in xrange(cores)]

		for p in processes:
			p.start()  # wait for execution to finish

		for p in processes:
			p.join()

		if QUEUE_ERROR:
			sys.exit("Error in queue")

		for combine in combine_queue:
			os.system(combine)

		my_m8s = my_dir + "*m8"
		os.system("cat " + my_m8s + " > " + m8)

		blast_finished_file = my_dir + "BLAST_FINISHED"
		bf = open(blast_finished_file, 'w')
		bf.write("Way to go!\n")
		bf.close()
	# sys.exit(0)


if __name__ == "__main__":
	if len(sys.argv) == 1:
		usage()
	else:
		main(sys.argv[1:])
