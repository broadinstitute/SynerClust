#!/usr/bin/env python

import sys
import os
import cPickle as pickle
import logging
import multiprocessing
import subprocess
import networkx as nx
import numpy
import argparse

DEVNULL = open(os.devnull, 'w')
MUSCLE_CMD = ["#MUSCLE_PATH", "-maxiters", "2", "-diags", "-sv", "-distance1", "kbit20_3", "-quiet"]
FASTTREE_CMD = ["#FASTTREE_PATH", "-quiet", "-nosupport"]
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


class Selector(multiprocessing.Process):
	def __init__(self, cluster_queue, result_queue, lock, dist_threshold, cons_out):
		multiprocessing.Process.__init__(self)
		self.cluster_queue = cluster_queue
		self.result_queue = result_queue
		self.lock = lock
		self.dist_threshold = dist_threshold
		self.cons_out = cons_out
		self.logger = logging.getLogger()

	def run(self):
		cons_res = {}
		while True:
			next_task = self.cluster_queue.get()
			if next_task is None:
				self.cluster_queue.task_done()
				break
			# compute
			next_task(self.lock, cons_res, self.dist_threshold, self.cons_out, self.logger)
			# compute finished
			self.cluster_queue.task_done()
		self.result_queue.put(cons_res)  # put result here


class Select(object):
	def __init__(self, clusterID, pep_data):
		self.clusterID = clusterID
		self.pep_data = pep_data

	def __call__(self, lock, cons_res, dist_threshold, cons_out, logger):
		(mus_out, output) = get_fasttree("".join(self.pep_data))
		lengths = {}
		seqs = {}
		leaves = []
		for line in self.pep_data:
			if line[0] == ">":
				pep = line[1:].split("\n")
				seqs[pep[0]] = pep[1]  # could most likely use this dict instead of recreating one from mus_out
				lengths[pep[0]] = int(pep[0].split(";")[-1])
				leaves.append(pep[0])
		graph = nx.Graph()
		counter = 1
		representative_sequences = []
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
		# remove from matrix data that is not needed anymore
		# redo on remaining sequences using the matrix with only their sequences (prune the big matrix)

		# add single non represented sequence left, if any
		if leaves:
			representative_sequences.append(leaves.pop())

		mus_out = mus_out.split("\n")
		mus_str_seqs = {}
		for mo in mus_out:
			l = mo.rstrip()  # probably not needed since split was used before to create it
			if l.find('>') > -1:
				seqID = l[1:]
				mus_str_seqs[seqID] = ""
			else:
				mus_str_seqs[seqID] += l

		lock.acquire()
		cons_res[self.clusterID] = []
		for s in representative_sequences:
			out_buffer = ""
			cons_seq = mus_str_seqs[s]
			cons_seq = cons_seq.replace("-", "")
			out_buffer += ">" + self.clusterID + ";" + str(len(cons_seq) + 1) + "\n"
			out_buffer += cons_seq + "*\n"
			cons_res[self.clusterID].append(out_buffer)
			cons_out.write(out_buffer)
		cons_out.flush()
		lock.release()


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

	if "CLUSTERS_REFINED" not in os.listdir(args.node_dir):
		sys.exit("Error: Previous step (cluster refinement) has not been completed.\n")

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

	with open(args.node_dir + "pep_data.pkl", "r") as sdat:
		pickleSeqs = pickle.load(sdat)

	consensus_pep = args.node_dir + args.node + ".pep"
	os.system("rm " + consensus_pep)  # since the file is opened in append mode every time, we need to delete anything from a previous run

	cons_out = open(consensus_pep, "a")
	with open(args.node_dir + "pre_consensus_data.pkl", "r") as f:
		cons_pkl = pickle.load(f)

	lock = multiprocessing.Lock()
	cluster_queue = multiprocessing.JoinableQueue()
	result_queue = multiprocessing.Queue()

	selectors = [Selector(cluster_queue, result_queue, lock, args.dist, cons_out) for i in xrange(args.numThreads)]
	for w in selectors:
		w.start()

	for clusterID in pickleSeqs:
		cluster_queue.put(Select(clusterID, pickleSeqs[clusterID]))

	for i in xrange(args.numThreads):
		cluster_queue.put(None)

	cluster_queue.join()

	for i in xrange(args.numThreads):
		cons_pkl.update(result_queue.get())

	cons_out.close()

	with open(args.node_dir + "consensus_data.pkl", "w") as f:
		pickle.dump(cons_pkl, f)

	logger.info("All consensus sequences accounted for.")
	os.system("cat " + args.node_dir + "clusters/singletons.cons.pep >> " + consensus_pep)  # ">>" to append
	logger.info("Made pep file")

	node_done_file = args.node_dir + "NODE_COMPLETE"
	with open(node_done_file, 'w') as cr:
		cr.write("Way to go!\n")
