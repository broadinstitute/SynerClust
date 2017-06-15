#!/usr/bin/env python

import sys
import os
import pickle
import logging
from multiprocessing import Process, Queue, RLock
from Queue import Empty
import subprocess
import networkx as nx
# import median_of_medians
import numpy
import argparse
import time

DEVNULL = open(os.devnull, 'w')
# DIST_THRESHOLD = 0.7
# REGEX = re.compile("(\([a-zA-Z0-9-_:;.]+,[a-zA-Z0-9-_:;.]+\))")
MUSCLE_CMD = ["#MUSCLE_PATH", "-maxiters", "2", "-diags", "-sv", "-distance1", "kbit20_3", "-quiet"]
# MUSCLE_CMD = ["/Users/cgeorges/Work/Tools/muscle3.8.31_i86darwin64", "-maxiters", "2", "-diags", "-sv", "-distance1", "kbit20_3", "-quiet"]  # kept for debugging
FASTTREE_CMD = ["#FASTTREE_PATH", "-quiet", "-nosupport"]
# FASTTREE_CMD = ["/Users/cgeorges/Work/Tools/FastTreeDouble", "-quiet", "-nosupport"]
OUTPUT_LOCK = RLock()
QUEUE_ERROR = False


def get_alignement(stdin_data):
	process = subprocess.Popen(MUSCLE_CMD, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=DEVNULL)
	output1 = process.communicate(stdin_data)[0]
	return output1


def get_fasttree(stdin_data):
	output1 = get_alignement(stdin_data)
	process = subprocess.Popen(FASTTREE_CMD, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=DEVNULL)
	output2 = process.communicate(output1)[0]
	return (output1, output2)


def makeConsensus(tq, resultsQueue, dist_threshold, consensus_pep):
	logger = logging.getLogger()
	cons_res = {}
	while True:
		try:
			nos = tq.get(block=False)
			pep_data = nos[0]
			clusterID = nos[1]
			logger.debug("Processing cluster " + clusterID)
			(mus_out, output) = get_fasttree("".join(pep_data))

			lengths = {}
			seqs = {}
			leaves = []
			for line in pep_data:
				if line[0] == ">":
					pep = line[1:].split("\n")
					seqs[pep[0]] = pep[1]  # could most likely use this dict instead of recreating one from mus_out
					lengths[pep[0]] = int(pep[0].split(";")[-1])
					leaves.append(pep[0])
			graph = nx.Graph()
			counter = 1
			representative_sequences = []
			logger.debug(output)
			while(True):
				r = output.find(")")
				l = output[:r].rfind("(")

				children = output[l + 1:r].split(",")
				if len(children) == 1:
					if len(graph.nodes()) == 0:
						representative_sequences.append(output.split(":")[0][1:])
					break
				group = "node" + str(counter)
				counter += 1

				if group not in graph.nodes():  # isn't it always a new node?
					graph.add_node(group)
				for child in children:
					child = child.split(":")
					if child[0] not in graph.nodes():
						graph.add_node(child[0])
					graph.add_edge(group, child[0], dist=float(child[1]))
				if "node" not in children[0] and "node" not in children[1] and children[0].split(":")[1] == 0.0 and seqs[children[0].split(":")[0]] != seqs[children[1].split(":")[0]]:
					logger.critical("0.0 distance in fasttree:\n" + output + "\n" + children[0].split(":")[0] + "\n" + seqs[children[0].split(":")[0]] + "\n" + children[1].split(":")[0] + "\n" + seqs[children[1].split(":")[0]])
				output = output[:l] + group + output[r + 1:]

			logger.debug("Built graph for cluster " + clusterID)
			matrix_size = len(leaves) * (len(leaves) - 1) / 2
			leaves_length = len(leaves)
			# nodes are headers/gene references
			dist_matrix = numpy.empty(matrix_size, float)  # diagonal matrix stored as an array
			j = 0
			i = 1
			for n in leaves[1:]:
				for m in leaves[:i]:
					dist_matrix[j] = nx.shortest_path_length(graph, n, m, "dist")
					j += 1
				i += 1

			while (matrix_size > 0):
				# select longest by reading dict once and keeping track of k,v for longest
				max_len = 0
				max_key = None
				for key, value in lengths.iteritems():
					if value > max_len:
						max_len = value
						max_key = key
					elif value == max_len:
						max_key = min(max_key, key)
				representative_sequences.append(max_key)
				index = leaves.index(max_key)
				# representative_sequences.append(leaves[index])
				logger.debug("Added representative to cluster " + clusterID)
				to_remove = numpy.full(leaves_length, -1, int)
				i = 0
				for j in xrange(leaves_length):
					if index == j:
						to_remove[i] = j
						del lengths[leaves[j - i]]
						leaves.remove(leaves[j - i])
						i += 1
					elif index < j:
						if dist_matrix[(j * (j - 1) / 2) + index] < dist_threshold:
							to_remove[i] = j
							del lengths[leaves[j - i]]
							leaves.remove(leaves[j - i])
							i += 1
					else:
						if dist_matrix[(index * (index - 1) / 2) + j] < dist_threshold:
							to_remove[i] = j
							del lengths[leaves[j - i]]
							leaves.remove(leaves[j - i])
							i += 1
				logger.debug("Distance check with " + str(i) + " represented sequences done for " + clusterID)

				leaves_length -= i
				k = 0
				l = 1
				offset = 0
				for pos in xrange(matrix_size):
					if k in to_remove or l in to_remove:
						offset += 1
					else:
						dist_matrix[pos - offset] = dist_matrix[pos]
					k += 1
					if not k < l:
						k = 0
						l += 1

				matrix_size = leaves_length * (leaves_length - 1) / 2
				logger.debug("Matrix updated for " + clusterID)
			# remove from matrix data that is not needed anymore
			# redo on remaining sequences using the matrix with only their sequences (prune the big matrix)

			# add single non represented sequence left, if any
			if leaves:
				representative_sequences.append(leaves.pop())

			mus_out = mus_out.split("\n")
			# mus_seqs = {}
			mus_str_seqs = {}
			# total_length = 0
			# unrep = []
			for mo in mus_out:
				l = mo.rstrip()  # probably not needed since split was used before to create it
				if l.find('>') > -1:
					seqID = l[1:]
					# line = l.split()   # why split??
					# seqID = line[0][1:line[0].rfind(";")]
					# trailer = 1
					# while seqID in mus_seqs:
					# 	seqID = line[0][1:] + "." + str(trailer)
					# 	trailer += 1
					# mus_seqs[seqID] = []
					# unrep.append(seqID)
					mus_str_seqs[seqID] = ""
				else:
					mus_str_seqs[seqID] += l
					# for i in l:
					# 	mus_seqs[seqID].append(i)
					# total_length += len(l)

			logger.debug("Trying to acquire lock for " + clusterID)
			OUTPUT_LOCK.acquire()
			logger.debug("Acquired lock for " + clusterID)
			cons_res[clusterID] = []
			for s in representative_sequences:
				# cons_seq = "".join(mus_seqs[s])
				out_buffer = ""
				cons_seq = mus_str_seqs[s]
				cons_seq = cons_seq.replace("-", "")
				out_buffer += ">" + clusterID + ";" + str(len(cons_seq)) + "\n"
				out_buffer += cons_seq + "*\n"
				cons_res[clusterID].append(out_buffer)
				cons_out.write(out_buffer)
				logger.debug("Wrote sequence for " + clusterID)
			cons_out.flush()
			logger.debug("Trying to release lock for " + clusterID)
			OUTPUT_LOCK.release()
			logger.debug("Released lock for " + clusterID)

		except Empty:
			resultsQueue.put(cons_res)
			break

		except Exception as e:
			QUEUE_ERROR = True
			print "Unexpected error in queue."
			raise e


if __name__ == "__main__":
	usage = "usage: WF_FinalizeNode_threaded.py [options]"
	parser = argparse.ArgumentParser(usage)
	parser.add_argument('-dir', dest="node_dir", required=True, help="Path to the \"nodes\" folder. (Required)")
	parser.add_argument('-node', dest="node", required=True, help="Current node name. (Required)")
	parser.add_argument('-t', type=int, dest="numThreads", default=4, help="Number of threads to use. (default = 4)")
	parser.add_argument('-dist', type=float, dest="dist", required=True, help="Branch distance threshold. (Required)")

	args = parser.parse_args()

	# after consensus sequences are generated, concatenates them into one big file, NODE.pep, and creates a NODE_COMPLETE file
	fileDir = args.node_dir + "clusters/"

	if "NODE_COMPLETE" in os.listdir(args.node_dir):
		sys.exit(0)

	FORMAT = "%(asctime)-15s %(levelname)s %(module)s.%(name)s.%(funcName)s at %(lineno)d :\n\t%(message)s\n"
	logger = logging.getLogger()
	logging.basicConfig(filename=args.node_dir + 'FinalizeNode_threaded.log', format=FORMAT, filemode='w', level=logging.DEBUG)
	# add a new Handler to print all INFO and above messages to stdout
	ch = logging.StreamHandler(sys.stdout)
	ch.setLevel(logging.INFO)
	logger.addHandler(ch)
	logger.info('Started')

	pklPep = args.node_dir + "pep_data.pkl"
	sdat = open(pklPep, 'rb')
	pickleSeqs = pickle.load(sdat)
	sdat.close()

	try:
		notOKQ = Queue(0)
		resultsQueue = Queue(0)
		consensus_pep = args.node_dir + args.node + ".pep"
		os.system("rm " + consensus_pep)  # since the file is opened in append mode every time, we need to delete anything from a previous run

		cons_out = open(consensus_pep, "a")
		with open(args.node_dir + "consensus_data.pkl", "r") as f:
			cons_pkl = pickle.load(f)

		processes = [Process(target=makeConsensus, args=(notOKQ, resultsQueue, args.dist, cons_out)) for i in xrange(args.numThreads)]

		for clusterID in pickleSeqs:
			notOKQ.put((pickleSeqs[clusterID], clusterID))
			time.sleep(0.01)

		for p in processes:
			p.start()
			logger.info("Starting %s" % (p.pid))

		for i in xrange(args.numThreads):
			cons_pkl.update(resultsQueue.get())

		cons_out.close()
	except Exception as e:
		logger.critical("Unexpected error:", sys.exc_info()[0])
		sys.exit(1)

	if QUEUE_ERROR:
		sys.exit("Error in python queue")
	# if len(pickleSeqs) != len(cons_pkl):  # shouldn't be equal because cons_pkl contains orphans that don't need representative selection
	# 	sys.exit("Error, missing results.")

	with open(args.node_dir + "consensus_data.pkl", "w") as f:
		pickle.dump(cons_pkl, f)

	logger.info("All consensus sequences accounted for.")
	os.system("cat " + args.node_dir + "clusters/singletons.cons.pep >> " + consensus_pep)  # ">>" to append
	logger.info("Made pep file")

	node_done_file = args.node_dir + "NODE_COMPLETE"
	cr = open(node_done_file, 'w')
	cr.write("Way to go!\n")
	cr.close()
