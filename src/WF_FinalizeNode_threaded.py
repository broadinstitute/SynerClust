#!/usr/bin/env python

import sys
import os
import pickle
import logging
import threading
from multiprocessing import Process, Queue
from Queue import Empty
import scipy.spatial.distance as distance
import subprocess

DEVNULL = open(os.devnull, 'w')

def usage():
	print "After consensus sequences are generated, concatenates them into one big file, NODE.pep, and creates a NODE_COMPLETE file."
	sys.exit(1)

def run_muscle(cmd, stdin_data):
	process = subprocess.Popen(cmd, stdin = subprocess.PIPE, stdout = subprocess.PIPE, stderr = DEVNULL)
	output = process.communicate(stdin_data)[0]
	return output

def makeConsensus(tq, hamm_dist, consensus_pep, output_lock):
	while True:
		try:
			nos = tq.get(block=True, timeout=3)
			pep_data = nos[0]
			clusterID = nos[1]
			muscle_cmd =["#MUSCLE_PATH", "-maxiters", "2", "-diags", "-sv", "-distance1", "kbit20_3", "-quiet"]
# 			muscle_cmd =["/home/kamigiri/tools/muscle3.8.31_i86linux64", "-maxiters", "2", "-diags", "-sv", "-distance1", "kbit20_3", "-quiet"]  # kept for debugging
			mus_out = run_muscle(muscle_cmd, "".join(pep_data)).split("\n")
			mus_seqs = {}
			mus_str_seqs = {}
			total_length = 0
			unrep = []
			for mo in mus_out:
				l = mo.rstrip()
				if l.find('>') > -1:
					line = l.split()
					seqID = line[0][1:]
					trailer = 1
					while seqID in mus_seqs:
						seqID = line[0][1:] + "." + str(trailer)
						trailer += 1
					mus_seqs[seqID] = []
					unrep.append(seqID)
					mus_str_seqs[seqID] = ""
				else:
					mus_str_seqs[seqID] += l
					for i in l:
						mus_seqs[seqID].append(i)
					total_length += len(l)
			output_lock.acquire()
			cons_out = open(consensus_pep, "a")
			while len(unrep) > 0:
				hams = {}
				minDist = float(len(unrep))
				minDistSeq = ""
				for ms in unrep:
					hams[ms] = {}
					m_m = mus_str_seqs[ms]
					m_len = len(m_m)
					m_start = m_len - len(m_m.lstrip("-"))
					m_stop = len(m_m.rstrip("-"))
					dist = 0
					for ns in unrep:
						hams[ms][ns] = distance.hamming(mus_seqs[ms][m_start:m_stop], mus_seqs[ns][m_start:m_stop])
						dist += distance.hamming(mus_seqs[ms], mus_seqs[ns])
					if dist < minDist or len(minDistSeq) == 0:
						minDist = dist
						minDistSeq = ms
				cons_seq = "".join(mus_seqs[minDistSeq])
				cons_seq = cons_seq.replace("-", "")
				cons_out.write(">" + clusterID + ";" + str(len(cons_seq)) + "\n")
				cons_out.write(cons_seq + "*\n")
				repped = []
				for ms in unrep:
					if hams[ms][minDistSeq] < hamm_dist:
						repped.append(ms)
				for r in repped:
					unrep.remove(r)
			cons_out.close()
			output_lock.release()
		except Empty:
			break

if __name__ == "__main__":
	argv = []
	if len(sys.argv) == 1:
		usage()
	else:
		argv = sys.argv[1:]
	# after consensus sequences are generated, concatenates them into one big file, NODE.pep, and creates a NODE_COMPLETE file
	node_dir = argv[0]
	node = argv[1]
	fileDir = node_dir + "clusters/"
	nThreads = int(argv[2])
	numThreads = min(nThreads, 4)
	chunk = 25
	h_dist = 0.6
	if len(argv) > 3:
		h_dist = float(argv[3])

	FORMAT = "%(asctime)-15s %(levelname)s %(module)s.%(name)s.%(funcName)s at %(lineno)d :\n\t%(message)s\n"
	logger = logging.getLogger()
	logging.basicConfig(filename=node_dir + 'FinalizeNode_threaded.log', format=FORMAT, filemode='w', level=logging.DEBUG)
	# add a new Handler to print all INFO and above messages to stdout
	ch = logging.StreamHandler(sys.stdout)
	ch.setLevel(logging.INFO)
	logger.addHandler(ch)
	logger.info('Started')

	pklPep = node_dir + "pep_data.pkl"
	sdat = open(pklPep, 'rb')
	pickleSeqs = pickle.load(sdat)
	sdat.close()

	notOKQ = Queue(0)
	output_lock = threading.Lock()
	consensus_pep = node_dir + node + ".pep"
	os.system("rm " + consensus_pep)  # since the file is opened in append mode every time, we need to delete anything from a previous run

	processes = [Process(target=makeConsensus, args=(notOKQ, h_dist, consensus_pep, output_lock)) for i in range(numThreads)]

	for clusterID in pickleSeqs:
		notOKQ.put((pickleSeqs[clusterID], clusterID))

	for p in processes:
		p.start()
		logger.info("Starting %s" % (p.pid))
		# print "Starting",p.pid
	for p in processes:
		logger.info("Starting %s" % (p.pid))
		# print "Stopping",p.pid
		p.join()

	logger.info("All consensus sequences accounted for.")
	os.system("cat " + node_dir + "clusters/singletons.cons.pep >> " + consensus_pep)  # ">>" to append 
	logger.info("Made pep file")

	node_done_file = node_dir + "NODE_COMPLETE"
	cr = open(node_done_file, 'w')
	cr.write("Way to go!\n")
	cr.close()
