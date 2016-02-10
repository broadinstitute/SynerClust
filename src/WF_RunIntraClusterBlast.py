#!/usr/bin/env python

import sys
import os
import pickle
import logging
import shutil


def usage():
	print """Launches BLAST after some brief pre-processing

	WF_RunBlast.py [node directory] [node name] [e value] [number of cores] [child 1] [child 2]
	"""
	sys.exit(1)


def main(argv):
	node_dir = argv[0]
	node = argv[1]
	# children = argv[4:]
	evalue = argv[2]
	cores = argv[3]
	my_dir = node_dir + node + "/"

	FORMAT = "%(asctime)-15s %(levelname)s %(module)s.%(name)s.%(funcName)s at %(lineno)d :\n\t%(message)s\n"
	logger = logging.getLogger()
	logging.basicConfig(filename=my_dir + 'RunInnerBlast.log', format=FORMAT, filemode='w', level=logging.DEBUG)
	# add a new Handler to print all INFO and above messages to stdout
	ch = logging.StreamHandler(sys.stdout)
	ch.setLevel(logging.INFO)
	logger.addHandler(ch)
	logger.info('Started')

	# run BLAST if BLAST_FINISHED doesn't exist in the node dir
	if "INNER_BLAST_FINISHED" not in os.listdir(my_dir):
		# fastas = []
		heads = []
		# my_head = my_dir + "blast_headers.txt"
		m8 = my_dir + "blast.m8"
		m8s = []
		# lengths = []
		new_cluster_dir = node_dir + node + "/clusters/"
		old_cluster_dir = node_dir + node + "/old_clusters/"
		blast_dir = node_dir + node + "/cluster_blast/"
		shutil.move(new_cluster_dir, old_cluster_dir)  # current "new" is the "old"
		try:
			os.mkdir(blast_dir)
		except OSError:  # folder already exists
			pass
			# shutil.rmtree(blast_dir)
			# os.mkdir(blast_dir)
		cluster_files = os.listdir(old_cluster_dir)
		cluster_files.remove("singletons.cons.pep")
		for cluster in cluster_files:
			# old_cluster_dir + cluster
			shutil.copy(old_cluster_dir + cluster, new_cluster_dir + cluster)
			# pfile = node_dir + node + "/" + cluster
			c_fasta = new_cluster_dir + cluster
			# fastas.append(c_fasta)
			c_head = new_cluster_dir + cluster + ".blast_headers.txt"
			heads.append(c_head)
			blast_out = new_cluster_dir + cluster + ".blast.m8"
			m8s.append(blast_out)
			formatDB_log = new_cluster_dir + cluster + "formatdb.log"
			# os.system("cp " + pfile + " " + c_fasta)
			os.system("grep '>' " + c_fasta + "| cut -f2 -d '>' > " + c_head)
			os.system("#BLAST_PATHformatdb -i " + c_fasta + " -l " + formatDB_log)

			# header_file = open(c_head, 'r').readlines()
			# headers = {}
			# for head in header_file:
			# 	head = head.rstrip()
			# 	line = head.split(";")
			# 	headers[line[0]] = int(line[1])
			# effective_length = 0
			# for head in headers:
			# 	effective_length += headers[head]
			# lengths.append(effective_length)

			# l_map_file = node_dir + c + "/locus_mappings.pkl"
			# pklFile = open(l_map_file, 'rb')
			# locusMap = pickle.load(pklFile)
			# pklFile.close()
			# strains = set([])
			# strains.add(c)
			# for lm in locusMap:
			# 	for lmm in locusMap[lm]:
			# 		strains.add("_".join(lmm.split("_")[:-1]))
			# 		logger.debug("%s splitted to %s" % (lmm, "_".join(lmm.split("_")[:-1])))
			# logger.info(strains)

			# if len(strains) == 1:
			self_blast_out = new_cluster_dir + cluster + "_self.blast.m8"
			os.system("#BLAST_PATHblastall -p blastp -m8 -e " + evalue + " -a " + cores + " -d " + c_fasta + " -i " + c_fasta + " -o " + self_blast_out)

		# cat_head_cmd = "cat " + heads[0] + " " + heads[1] + " > " + my_head
		# logger.info(cat_head_cmd)
		# # os.system(cat_head_cmd)
		# os.system("#BLAST_PATHblastall -p blastp -m8 -e " + evalue + " -a " + cores + " -d " + fastas[0] + " -i " + fastas[1] + " -o " + m8s[1])
		# os.system("#BLAST_PATHblastall -p blastp -m8 -e " + evalue + " -a " + cores + " -d " + fastas[1] + " -i " + fastas[0] + " -o " + m8s[0])

		my_m8s = my_dir + "*m8"
		os.system("cat " + my_m8s + " > " + m8)

		blast_finished_file = my_dir + "INNER_BLAST_FINISHED"
		bf = open(blast_finished_file, 'w')
		bf.write("Way to go!\n")
		bf.close()
	# sys.exit(0)


if __name__ == "__main__":
	if len(sys.argv) == 1:
		usage()
	else:
		main(sys.argv[1:])
