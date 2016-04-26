#!/usr/bin/env python

import sys
import os
import pickle
import logging
from multiprocessing import Process, Queue, RLock
from Queue import Empty
import scipy.spatial.distance as distance
import subprocess
import networkx as nx
import median_of_medians
import numpy

DEVNULL = open(os.devnull, 'w')
DIST_THRESHOLD = 0.7
# REGEX = re.compile("(\([a-zA-Z0-9-_:;.]+,[a-zA-Z0-9-_:;.]+\))")
# MUSCLE_CMD =["#MUSCLE_PATH", "-maxiters", "2", "-diags", "-sv", "-distance1", "kbit20_3", "-quiet"]
MUSCLE_CMD = ["/home/kamigiri/tools/muscle3.8.31_i86linux64", "-maxiters", "2", "-diags", "-sv", "-distance1", "kbit20_3", "-quiet"]  # kept for debugging
# FASTTREE_CMD = ["#FASTTREE_PATH", "-quiet"]
FASTTREE_CMD = ["/home/kamigiri/tools/FastTreeDouble", "-quiet", "-nosupport"]
OUTPUT_LOCK = RLock()


def usage():
	print "After consensus sequences are generated, concatenates them into one big file, NODE.pep, and creates a NODE_COMPLETE file."
	sys.exit(1)


def get_fasttree(stdin_data):
	process = subprocess.Popen(MUSCLE_CMD, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=DEVNULL)
	output1 = process.communicate(stdin_data)[0]
	process = subprocess.Popen(FASTTREE_CMD, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=DEVNULL)
	output2 = process.communicate(output1)[0]
	return (output1, output2)


def makeConsensus(tq, hamm_dist, consensus_pep):
	logger = logging.getLogger()
	while True:
		try:
			nos = tq.get(block=True, timeout=3)
			pep_data = nos[0]
			clusterID = nos[1]
			logger.info("Processing cluster " + clusterID)
			(mus_out, output) = get_fasttree("".join(pep_data))

			graph = nx.Graph()
			counter = 1
			leaves = []
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
						leaves.append(child[0])
					graph.add_edge(group, child[0], dist=float(child[1]))
				output = output[:l] + group + output[r + 1:]
			
			logger.info("Built graph for cluster " + clusterID)
			matrix_size = len(leaves) * (len(leaves) - 1) / 2
			while (matrix_size > 0):
				# nodes are headers/gene references
				leaves.sort()  # why?
				dist_matrix = numpy.empty(matrix_size, float)  # diagonal matrix stored as an array
				j = 0
				i = 1
				for n in leaves[1:]:
					for m in leaves[:i]:
						dist_matrix[j] = nx.shortest_path_length(graph, n, m, "dist")
						j += 1
					i += 1

				# n rows, 2 columns
				# 1st column = values, 2nd column = index
				# average_dist is actually a sum of dist
				# but since dividing all of its values by len(leaves) wouldn't change the median and we don't use these values elsewhere there is no need to do it
				average_dist = numpy.zeros((len(leaves), 2), float)
				# for i in xrange(len(nodes)):
				# 	average_dist[i][1] = i
				# 	for j in xrange(len(nodes)):  # still reads each cell twice, can try to see if more efficient to find to which 2 average_dist to add each cell
				# 		if i == j:
				# 			continue
				# 		elif i < j:
				# 			average_dist[i][0] += dist_matrix[(j * (j - 1) / 2) + i]
				# 		else:
				# 			average_dist[i][0] += dist_matrix[(i * (i - 1) / 2) + j]

				# read through the matrix only once to fill average_dist
				imax = 1
				pos = 0
				while(pos < matrix_size):
					average_dist[imax-1][1] = imax-1
					for i in xrange(imax):
						average_dist[i][0] += dist_matrix[pos]
						average_dist[imax][0] += dist_matrix[pos]
						pos += 1
					imax += 1

				# pick something close to the median (median of medians)
				# for_med = numpy.copy(average_dist)  # the median of medians search reorders elements
				index = int(average_dist[median_of_medians.for2DArray(average_dist)][1])
				representative_sequences.append(leaves[index])
				logger.info("Added representative to cluster " + clusterID)
				# n = math.floor((math.sqrt(8 * index + 1) - 1) / 2)
				# check which sequences are more distant to it than the threshold
				to_remove = numpy.full(len(leaves), -1, int)
				old_length = len(leaves)
				i = 0
				for j in xrange(old_length):
					if index == j:
						to_remove[i] = j
						leaves.remove(leaves[j-i])
						i += 1
					elif index < j:
						if dist_matrix[(j * (j - 1) / 2) + index] < DIST_THRESHOLD:
							to_remove[i] = j
							leaves.remove(leaves[j-i])
							i += 1
					else:
						if dist_matrix[(index * (index - 1) / 2) + j] < DIST_THRESHOLD:
							to_remove[i] = j
							leaves.remove(leaves[j-i])
							i += 1

				new_length = old_length - i
				new_dist_matrix = numpy.empty(((new_length - 1) * new_length) / 2, float)

				# -1 or -inf in old matrix where needed
				i = 0
				while(i < len(to_remove) and to_remove[i] != -1):
					t = to_remove[i] * (to_remove[i] - 1) / 2
					step = to_remove[i]
					for j in xrange(t, t + step):
						dist_matrix[j] = -float('Inf')
					pos = t + step - 1
					while (step < old_length - 1):
						step += 1
						pos += step
						dist_matrix[pos] = -float('Inf')
					i += 1

				i = 0
				for d in numpy.nditer(dist_matrix):
					if d != -float('Inf'):
						new_dist_matrix[i] = d
						i += 1

				dist_matrix = new_dist_matrix
				matrix_size = len(leaves) * (len(leaves) - 1) / 2
			# -float('Inf')
			# remove from matrix data that is not needed anymore
			# redo on remaining sequences using the matrix with only their sequences (prune the big matrix)

			mus_out = mus_out.split("\n")
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

			OUTPUT_LOCK.acquire()
			cons_out = open(consensus_pep, "a")
			for s in representative_sequences:
				cons_seq = "".join(mus_seqs[s])
				cons_seq = cons_seq.replace("-", "")
				cons_out.write(">" + clusterID + ";" + str(len(cons_seq)) + "\n")
				cons_out.write(cons_seq + "*\n")
			cons_out.close
			OUTPUT_LOCK.release()

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
	consensus_pep = node_dir + node + ".pep"
	os.system("rm " + consensus_pep)  # since the file is opened in append mode every time, we need to delete anything from a previous run

# 	with open("all_clusters_cons_pep_data.pkl") as f:
# 		all_pep_pkl = pickle.load(f)

	processes = [Process(target=makeConsensus, args=(notOKQ, h_dist, consensus_pep)) for i in range(numThreads)]

	for clusterID in pickleSeqs:
		notOKQ.put((pickleSeqs[clusterID], clusterID))

	for p in processes:
		p.start()
		logger.info("Starting %s" % (p.pid))
		# print "Starting",p.pid
	for p in processes:
		# print "Stopping",p.pid
		p.join()
		logger.info("Finished %s" % (p.pid))

	logger.info("All consensus sequences accounted for.")
	os.system("cat " + node_dir + "clusters/singletons.cons.pep >> " + consensus_pep)  # ">>" to append
	logger.info("Made pep file")

	node_done_file = node_dir + "NODE_COMPLETE"
	cr = open(node_done_file, 'w')
	cr.write("Way to go!\n")
	cr.close()
