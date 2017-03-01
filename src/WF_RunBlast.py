#!/usr/bin/env python
import sys
import os
import pickle
import logging


def usage():
	print """Launches BLAST after some brief pre-processing

	WF_RunBlast.py [node directory] [node name] [e value] [number of cores] [child 1] [child 2]
	"""
	sys.exit(1)


def main(argv):
	node_dir = argv[0]
	node = argv[1]
	children = argv[4:]
	evalue = argv[2]
	cores = argv[3]
	my_dir = node_dir + node + "/"

	FORMAT = "%(asctime)-15s %(levelname)s %(module)s.%(name)s.%(funcName)s at %(lineno)d :\n\t%(message)s\n"
	logger = logging.getLogger()
	logging.basicConfig(filename=my_dir + 'RunBlast.log', format=FORMAT, filemode='w', level=logging.DEBUG)
	# add a new Handler to print all INFO and above messages to stdout
	ch = logging.StreamHandler(sys.stdout)
	ch.setLevel(logging.INFO)
	logger.addHandler(ch)
	logger.info('Started')

	# run BLAST if BLAST_FINISHED doesn't exist in the node dir
	if "BLAST_FINISHED" not in os.listdir(my_dir):
		fastas = []
		heads = []
		my_head = my_dir + "blast_headers.txt"
		m8 = my_dir + "blast.m8"
		m8s = []
		lengths = []
		for c in children:
			pfile = node_dir + c + "/" + c + ".pep"
			c_fasta = my_dir + c + ".blast.fa"
			fastas.append(c_fasta)
			c_head = my_dir + c + ".blast_headers.txt"
			heads.append(c_head)
			blast_out = my_dir + c + ".blast.m8"
			m8s.append(blast_out)
			formatDB_log = my_dir + "formatdb.log"
			os.system("cp " + pfile + " " + c_fasta)
			os.system("grep '>' " + c_fasta + "| cut -f2 -d '>' > " + c_head)
			# os.system("#BLAST_PATHformatdb -i " + c_fasta + " -l " + formatDB_log)
			os.system("#BLAST_PATHmakeblastdb -in " + c_fasta + " -dbtype prot -logfile " + formatDB_log)

			header_file = open(c_head, 'r').readlines()
			headers = {}
			for head in header_file:
				head = head.rstrip()
				line = head.split(";")
				headers[line[0]] = int(line[1])
			effective_length = 0
			for head in headers:
				effective_length += headers[head]
			lengths.append(effective_length)

			l_map_file = node_dir + c + "/locus_mappings.pkl"
			pklFile = open(l_map_file, 'rb')
			locusMap = pickle.load(pklFile)
			pklFile.close()
			strains = set([])
			strains.add(c)
			for lm in locusMap:
				for lmm in locusMap[lm]:
					strains.add("_".join(lmm.split("_")[:-1]))
			logger.info(strains)

			if len(strains) == 1:
				self_blast_out = my_dir + c + "_self.blast.m8"
				# os.system("#BLAST_PATHblastall -p blastp -b 6 -m8 -e " + evalue + " -a " + cores + " -d " + c_fasta + " -i " + c_fasta + " -o " + self_blast_out)
				os.system("#BLAST_PATHblastp -num_alignments 6 -outfmt 6 -evalue " + evalue + " -num_threads " + cores + " -db " + c_fasta + " -query " + c_fasta + " -out " + self_blast_out)

		cat_head_cmd = "cat " + heads[0] + " " + heads[1] + " > " + my_head
		logger.info(cat_head_cmd)
		os.system(cat_head_cmd)
		# os.system("#BLAST_PATHblastall -p blastp -m8 -e " + evalue + " -a " + cores + " -d " + fastas[0] + " -i " + fastas[1] + " -o " + m8s[1])
		os.system("#BLAST_PATHblastp -outfmt 6 -evalue " + evalue + " -num_threads " + cores + " -db " + fastas[0] + " -query " + fastas[1] + " -out " + m8s[1])
		# os.system("#BLAST_PATHblastall -p blastp -m8 -e " + evalue + " -a " + cores + " -d " + fastas[1] + " -i " + fastas[0] + " -o " + m8s[0])
		os.system("#BLAST_PATHblastp -outfmt 6 -evalue " + evalue + " -num_threads " + cores + " -db " + fastas[1] + " -query " + fastas[0] + " -out " + m8s[0])

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
