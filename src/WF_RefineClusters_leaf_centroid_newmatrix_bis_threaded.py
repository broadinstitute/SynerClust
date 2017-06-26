#!/usr/bin/env python

import sys
import os
# import NJ
import cPickle as pickle
import networkx as nx
import NetworkX_Extension as nxe
import numpy
import logging
# import subprocess
import argparse
import math
from operator import itemgetter
import multiprocessing

DEVNULL = open(os.devnull, 'w')
BEST_HIT_PROPORTION_THRESHOLD = 0.95
SYNTENY_THRESHOLD = 0.3
SYNTENY_DIFFERENCE_THRESHOLD = 0.2
synteny_data = {}
gene_to_rough_cluster = {}


def usage():
	print """From the rough cluster trees generated by WF_MakeRoughClusters, refine clusters so that all duplication events (paralogs) \
	occur after speciation from the most recent common ancestor or MRCA, which is [node]. [node_dir] is the directory that contains all \
	data about [node]. [flow_id] refers to a step in the workflow; [jobs_per_cmd] is the number of consensus sequence computations \
	distributed to a single node.

	WF_RefineClusters.py [node_dir] [node] [children...]
	"""
	sys.exit(1)


# class safeCounter(object):
# 	def __init__(self, initval=1):
# 		self.val = manager.Value("i", initval)
# 		self.lock = manager.Lock()

# 	def increment_after(self):
# 		with self.lock:
# 			res = self.val.value
# 			self.val.value += 1
# 		return self.val.value

# 	@property
# 	def value(self):
# 		with self.lock:
# 			return self.val.value

class Counter(object):
	def __init__(self, val=0):
		self.val = multiprocessing.Value('i', val)
		self.lock = multiprocessing.Lock()

	def safeIncrement(self):
		with self.lock:
			self.val.value += 1
		return self.val.value - 1


class Refinery(multiprocessing.Process):
	def __init__(self, cluster_queue, result_queue, mrca, cluster_counter):
		multiprocessing.Process.__init__(self)
		self.cluster_queue = cluster_queue
		self.result_queue = result_queue
		self.mrca = mrca
		self.cluster_counter = cluster_counter

	def run(self):
		genes_to_cluster = {}
		ok_trees = []
		potentials = {}
		identical_orphans_to_check = []
		identical_orphans_to_check_dict = {}
		identical_index = 0
		while True:
			next_task = self.cluster_queue.get()
			if next_task is None:
				self.cluster_queue.task_done()
				break
			# compute
			identical_index = next_task(self.mrca, genes_to_cluster, self.cluster_counter, ok_trees, identical_orphans_to_check, identical_orphans_to_check_dict, identical_index, potentials)
			# compute finished
			self.cluster_queue.task_done()
		# print "thread finished with " + str(len(identical_orphans_to_check))
		self.result_queue.put((genes_to_cluster, ok_trees, identical_orphans_to_check, identical_orphans_to_check_dict, potentials))  # put result here


class Refine(object):
	def __init__(self, cluster, graph):
		self.cluster = cluster
		self.graph = graph
		# self.mrca = mrca
		# self.identical_orphans_to_check = identical_orphans_to_check
		# self.identical_orphans_to_check_dict = identical_orphans_to_check_dict
		# self.identical_index = 0
		# self.identical_index_lock = identical_index_lock
		# self.cluster_counter = cluster_counter
		# self.cluster_counter_lock = cluster_counter_lock
		# self.ok_trees = ok_trees
		# self.genes_to_cluster = genes_to_cluster
		# self.gene_to_rough_cluster = gene_to_rough_cluster
		# self.potentials = potentials

	def __call__(self, mrca, genes_to_cluster, cluster_counter, ok_trees, identical_orphans_to_check, identical_orphans_to_check_dict, identical_index, potentials):
		leaves = self.graph.nodes()
		leaves.sort()
		syn = {}
		for n in leaves:  # genes
			leaf = "_".join(n.split("_")[:-1])
			syn[n] = [[], synteny_data[leaf][n]['count']]
			for m in synteny_data[leaf][n]['neighbors']:
				syn[n][0].append(gene_to_rough_cluster[m])
		syn_matrix = numpy.empty(len(leaves) * (len(leaves) - 1) / 2)
		i = 1
		pos = 0
		for m in leaves[1:]:
			syn_m = set(syn[m][0])
			count_m = float(syn[m][1])
			syn_m.discard(self.cluster)
			mSeqs = (len(syn[m][0]) - syn[m][0].count(self.cluster)) / count_m  # divide by count of leaves
			for n in leaves[:i]:
				syn_n = set(syn[n][0])
				count_n = float(syn[n][1])
				syn_n.discard(self.cluster)
				nSeqs = (len(syn[n][0]) - syn[n][0].count(self.cluster)) / count_n  # divide by count of leaves
				matches = 0
				if mSeqs == 0 or nSeqs == 0:
					syn_matrix[pos] = 1.0  # no neighbors in common if someone has no neighbors  # -= 0 ? does it change anything?
					pos += 1
					continue
				# all_neighbors = syn_m & set(syn[n])  # no need to .discard(cluster) since already did in syn_m, so won't be present in union
				all_neighbors = syn_m | syn_n
				for a in all_neighbors:
					t_m = max(syn[m][0].count(a), 0) / count_m  # divide by count so that merging a level N nodes with a leaf doesn't give a max synteny of 1/N only
					t_n = max(syn[n][0].count(a), 0) / count_n
					matches += min(t_m, t_n)
				# synFrac = float(matches) / float(max(mSeqs, nSeqs))
				synFrac = float(matches) / float(mSeqs + nSeqs - matches)
				# synFrac = float(matches) / float(max_neighbors_count)
				syn_matrix[pos] = 1.0 - synFrac
				pos += 1
			i += 1

		# logger.debug("Built matrices for " + cluster + " in " + str(time.time() - TIMESTAMP))
		# formatting matrix for output
		i = 0
		j = 1
		syn_buff = leaves[0] + "\n" + leaves[1] + "\t"
		for y in numpy.nditer(syn_matrix):
			syn_buff += str(y) + "\t"
			i += 1
			if i >= j:
				i = 0
				j += 1
				if j < len(leaves):
					syn_buff += "\n" + leaves[j] + "\t"
				else:
					syn_buff += "\n"

		new_graph = self.graph.copy()
		# check synteny matrix for cells lowest synteny (below a threshold, 0.2-0.5?)
		syntenic = []
		it = numpy.nditer(syn_matrix, flags=['f_index'])
		while not it.finished:
			if it[0] <= SYNTENY_THRESHOLD:
				position = it.index + 1  # formula works for indexes starting at 1, so need to offset
				j = int(round(math.ceil(math.sqrt(2 * position + 0.25) - 0.5)))  # formula is for lower triangular matrix, so need to offset distance matrix columns because we start at 0
				i = int(position - ((j - 1) * j / 2))  # no need for row offset because row 0 is empty in distance matrix
				syntenic.append((it[0], it.index, i - 1, j))
			it.iternext()
		syntenic.sort(key=itemgetter(0))  # key precised so that sort is only done on first element of lists and not on other ones for potential ties

		# for pair in syntenic:  # loop twice to allow pairs where best hit got pulled into another pair to be clustered on 2nd round?
		k = 0
		while k < len(syntenic):
			pair = syntenic[k]
			i = pair[2]
			j = pair[3]
			if not self.graph.has_edge(leaves[i], leaves[j]):  # if (i, j) is in the graph, (j, i) is also per construction/filtering of rough clusters; reverse is true too
				k += 1
				continue
			# check if this cell is the only low synteny for each member of the pair
			if len([p for p in syntenic if p[2] == i or p[3] == i or p[2] == j or p[3] == j]) == 1:
				# if yes, check if RBH hits, or among best hits
				if self.graph[leaves[i]][leaves[j]]['rank'] + self.graph[leaves[j]][leaves[i]]['rank'] <= 4:  # 2 = RBH, put 3 or 4 as limit?
					# if yes, cluster
					# merge leaves[i] and leaves[j]
					syn_dist = ":" + str(pair[0] / 2.0)
					new_node = "%s_%06d" % (mrca, cluster_counter.safeIncrement())
					# self.cluster_counter += 1
					ok_trees.append((new_node, (leaves[i], leaves[j]), ("(" + leaves[i] + ":" + str(self.graph[leaves[i]][leaves[j]]['rank']) + "," + leaves[j] + ":" + str(self.graph[leaves[j]][leaves[i]]['rank']) + ")", "(" + leaves[i] + syn_dist + "," + leaves[j] + syn_dist + ")")))
					nxe.merge(new_graph, self.graph, leaves[i], leaves[j], new_node)
					genes_to_cluster[leaves[i]] = (new_node, True)
					genes_to_cluster[leaves[j]] = (new_node, True)
					# remove other edges pointing to those nodes
					self.graph.remove_node(leaves[i])
					self.graph.remove_node(leaves[j])
			# if one of the genes has more than 1 syntenic gene, but still RBH
			else:
				# if best hit for ones for which it isn't have -m close, cluster
				# check first if best hit hasn't been clustered
				l = k + 1
				ties = [k]
				s = set(pair[2:])
				while l < len(syntenic) and pair[0] == syntenic[l][0]:  # syntenic tie
					if syntenic[l][2] in s or syntenic[l][3] in s:
						ties.append(l)
					l += 1
				ties_results = []
				for l in ties:
					good = 0
					sum_of_ranks = self.graph[leaves[i]][leaves[j]]['rank'] + self.graph[leaves[j]][leaves[i]]['rank']
					sum_of_m = self.graph[leaves[i]][leaves[j]]['m'] + self.graph[leaves[j]][leaves[i]]['m']
					if self.graph[leaves[i]][leaves[j]]['rank'] == 1:
						good += 1
					elif max([f['rank'] for f in self.graph[leaves[i]].values()]) == self.graph[leaves[i]][leaves[j]]['rank']:  # best hit has been clustered, this is the best remaining hit
						good += 1
					elif [f for f in self.graph[leaves[i]].values() if f['rank'] == 1]:  # if best hit hasn't been clustered
						# check if better hits are syntenic?
						if self.graph[leaves[i]][leaves[j]]['m'] >= BEST_HIT_PROPORTION_THRESHOLD:
							good += 1
					if self.graph[leaves[j]][leaves[i]]['rank'] == 1:
						good += 1
					elif max([f['rank'] for f in self.graph[leaves[j]].values()]) == self.graph[leaves[j]][leaves[i]]['rank']:
						good += 1
					elif [f for f in self.graph[leaves[j]].values() if f['rank'] == 1]:  # if best hit hasn't been clustered
						if self.graph[leaves[j]][leaves[i]]['m'] >= BEST_HIT_PROPORTION_THRESHOLD:
							good += 1
					ties_results.append((l, good, sum_of_ranks, sum_of_m))
				ties_results.sort(key=lambda tup: (-tup[1], tup[2], -tup[3]))
				if ties_results[0][1] == 2 and (len(ties_results) == 1 or (len(ties_results) > 1 and (ties_results[0][1] > ties_results[1][1] or ties_results[0][2] < ties_results[1][2] or ties_results[0][3] > ties_results[1][3]))):  # triple "or" because results are sorted, so if not better, equal
						pair = syntenic[ties_results[0][0]]
						i = pair[2]
						j = pair[3]
						# merge leaves[i] and leaves[j]
						syn_dist = ":" + str(pair[0] / 2.0)
						new_node = "%s_%06d" % (mrca, cluster_counter.safeIncrement())
						# self.cluster_counter += 1
						ok_trees.append((new_node, (leaves[i], leaves[j]), ("(" + leaves[i] + ":" + str(self.graph[leaves[i]][leaves[j]]['rank']) + "," + leaves[j] + ":" + str(self.graph[leaves[j]][leaves[i]]['rank']) + ")", "(" + leaves[i] + syn_dist + "," + leaves[j] + syn_dist + ")")))
						nxe.merge(new_graph, self.graph, leaves[i], leaves[j], new_node)
						genes_to_cluster[leaves[i]] = (new_node, True)
						genes_to_cluster[leaves[j]] = (new_node, True)
						# remove other edges pointing to those nodes
						self.graph.remove_node(leaves[i])
						self.graph.remove_node(leaves[j])
				ties_results.reverse()
				for tr in ties_results:
					del syntenic[tr[0]]
				continue
			k += 1

		# check for remaining RBH and cluster
		changed = True
		while changed:
			changed = False
			i = 0
			nodes_left = [n for n in self.graph.nodes() if self.graph[n]]  # nodes left that still have edges connecting them to other nodes
			while i < len(nodes_left):
				n1 = nodes_left[i]
				pair = None
				targets = [n2 for n2 in self.graph[n1] if self.graph[n1][n2]['rank'] == 1 and self.graph[n2][n1]['rank'] == 1 and n1[:32] != n2[:32]]  # RBH, and don't merge self-hits from leaves
				if len(targets) == 0:
					i += 1
					continue
				if len(targets) >= 1:
					pairs = []
					for n2 in targets:
						ii = leaves.index(n1)
						jj = leaves.index(n2)
						syn = 1.0
						if ii < jj:
							syn = syn_matrix[(jj * (jj - 1) / 2) + ii]
						else:
							syn = syn_matrix[(ii * (ii - 1) / 2) + jj]
						pairs.append([n2, syn])  # target, synteny
					pairs.sort(key=itemgetter(1))  # sort by ascending synteny distance
					# if pairs[0][1] < 1.0 and (len(pairs) == 1 or pairs[0][1] != pairs[1][1]):  # synteny evidance and no ex-aequo
					if pairs[0][1] < 1.0 and (len(pairs) == 1 or pairs[1][1] - pairs[0][1] > SYNTENY_DIFFERENCE_THRESHOLD):
						# CHECK IF THE 2ND NODE ALSO HAS THE LOWEST SYNTENY WITH THE CURRENT NODE
						# CHECK SYNTENY MATRIX AT 2ND NODE VALUES FIRST? REVERT INDEX TO CHECK IF EDGE EXISTS AND GOOD HIT?
						likely_pair = pairs[0][0]
						targets2 = [n2 for n2 in self.graph[likely_pair] if self.graph[likely_pair][n2]['rank'] == 1 and self.graph[n2][likely_pair]['rank'] == 1 and likely_pair[:32] != n2[:32]]
						if len(targets2) > 1:  # else it is the only hit so good (because reciprocity previously checked)
							pairs2 = []
							for n2 in targets2:
								ii = leaves.index(likely_pair)
								jj = leaves.index(n2)
								syn = 1.0
								if ii < jj:
									syn = syn_matrix[(jj * (jj - 1) / 2) + ii]
								else:
									syn = syn_matrix[(ii * (ii - 1) / 2) + jj]
								pairs2.append([n2, syn])  # target, synteny
							pairs2.sort(key=itemgetter(1))  # sort by ascending synteny distance
							if pairs2[0][0] == n1 and pairs2[1][1] - pairs2[0][1] > SYNTENY_DIFFERENCE_THRESHOLD:  # synteny evidance to 1st node and no ex-aequo
								pair = pairs[0][0]
							else:  # best node to merge from n1 side is not the best from n2 side
								i += 1
								continue
						else:
							pair = pairs[0][0]
					else:  # no evidance of which node is the good one to merge to
						i += 1
						continue
				else:
					pair = targets[0]
					# if graph[e[0]][e[1]]['rank'] == 1 and graph[e[1]][e[0]]['rank'] == 1:
					# MERGE e[0] and e[1]
				changed = True
				ii = leaves.index(n1)
				jj = leaves.index(pair)
				ma = max(ii, jj)
				mi = min(ii, jj)
				pos = (ma * (ma - 1) / 2) + mi
				syn_dist = ":" + str(syn_matrix[pos] / 2.0)
				new_node = "%s_%06d" % (mrca, cluster_counter.safeIncrement())
				# self.cluster_counter += 1
				ok_trees.append((new_node, (n1, pair), ("(" + n1 + ":" + str(self.graph[n1][pair]['rank']) + "," + pair + ":" + str(self.graph[pair][n1]['rank']) + ")", "(" + n1 + syn_dist + "," + pair + syn_dist + ")")))
				nxe.merge(new_graph, self.graph, n1, pair, new_node)
				genes_to_cluster[n1] = (new_node, True)
				genes_to_cluster[pair] = (new_node, True)
				# remove other edges pointing to those nodes
				self.graph.remove_node(n1)
				self.graph.remove_node(pair)
				if nodes_left.index(pair) < i:
					i -= 1
				nodes_left.remove(n1)
				nodes_left.remove(pair)

		# check remaining, is there any match in between them directly, any synteny?
		# remaining = all nodes in graph, check what are their edges in new_graph = potential inparalogs
			# if none, keep as potential inparalogs but don't cluster

		# current_cluster_identical_orphans = {}
		for node in new_graph.nodes():
			if mrca not in node:
				if node not in genes_to_cluster:
					new_orphan = "%s_%06d" % (mrca, cluster_counter.safeIncrement())
					# self.cluster_counter += 1
					genes_to_cluster[node] = (new_orphan, False)

		# identical_index = 0
		for node in new_graph.nodes():
			if mrca not in node:
				new_orphan = genes_to_cluster[node][0]
				(k, v) = min([[f, new_graph[node][f]['rank']] for f in new_graph[node]], key=itemgetter(1))  # add secondary ranking by synteny?
				if mrca in k:
					if new_orphan not in potentials:
						potentials[new_orphan] = [k]
					else:  # never?
						if k in potentials[new_orphan]:
							potentials[new_orphan].append(k)
				elif node[:32] == k[:32]:  # self hit, possible for full species tree leaves only
					new_orphan2 = genes_to_cluster[k][0]
					if new_orphan not in potentials:
						potentials[new_orphan] = [new_orphan2]
					else:  # never?
						if new_orphan2 in potentials[new_orphan]:
							potentials[new_orphan].append(new_orphan2)
				ok_trees.append((new_orphan, (node,), (node, node)))  # (node,) comma is required so its a tuple that can be looped on and not on the string itself
				# find all identical orphans to the current one, which should be the full set for this sequence since they would all be linked to each other
				identicals = [(f, genes_to_cluster[f][0]) for f in new_graph[node] if mrca not in f and new_graph[node][f]['identity'] == 1 and node in new_graph[f] and new_graph[f][node]['identity'] == 1]
				if identicals:
					identicals.append((node, new_orphan))
					if identicals[0][1] not in identical_orphans_to_check_dict:  # else this set has already been found by another member of the set
						identical_orphans_to_check.append(identicals)
						# with self.identical_index_lock:
						for identical in identicals:  # identical = tuple(child_id, parent_id)
							identical_orphans_to_check_dict[identical[1]] = identical_index
						identical_index += 1
		return identical_index


def main():
	usage = "usage: WF_RefineCluster_leaf_centroid_newmatrix.py [options]"
	parser = argparse.ArgumentParser(usage)
	parser.add_argument('-dir', dest="node_dir", required=True, help="Path to the \"nodes\" folder. (Required)")
	parser.add_argument('-node', dest="node", required=True, help="Current node name. (Required)")
	parser.add_argument('-t', type=int, dest="numThreads", default=4, help="Number of threads to use. (default = 4)")
	parser.add_argument('--no-synteny', dest="synteny", default=True, action='store_false', required=False, help="Disable use of synteny (required is information not available).")
	parser.add_argument('children', nargs=2, help="Children nodes. (Required)")
	args = parser.parse_args()

	repo_path = args.node_dir[:-6]
	mrca = args.node

	my_dir = args.node_dir + mrca + "/"

	if "TREES_FINISHED" not in os.listdir(my_dir):
		sys.exit("Error: Previous step (rough clustering) has not been completed.\n")

	if "CLUSTERS_REFINED" in os.listdir(my_dir):
		sys.exit(0)

	FORMAT = "%(asctime)-15s %(levelname)s %(module)s.%(name)s.%(funcName)s at %(lineno)d :\n\t%(message)s\n"
	logger = logging.getLogger()
	logging.basicConfig(filename=my_dir + 'RefineClusters_leaf_centroid.log', format=FORMAT, filemode='w', level=logging.DEBUG)
	# add a new Handler to print all INFO and above messages to stdout
	ch = logging.StreamHandler(sys.stdout)
	ch.setLevel(logging.INFO)
	logger.addHandler(ch)
	logger.info('Started')

	# read trees, resolve clusters
	# tree_dir = my_dir + "trees/"
	cluster_dir = my_dir + "clusters"
	if "clusters" in os.listdir(my_dir):
		if "old" not in os.listdir(my_dir):
			os.system("mkdir " + my_dir + "old")

		os.system("mv -f " + cluster_dir + "/ " + my_dir + "old/")
	os.system("mkdir " + cluster_dir)
	cluster_dir = cluster_dir + "/"

	# synteny_data = manager.dict(lock=False)
	# cluster_counter = safeCounter(1, manager)
	# ok_trees = manager.list()
	# genes_to_cluster = manager.dict()
	# potentials = manager.dict()  # potential inparalogs

	# identical_orphans_to_check = manager.list()
	# identical_orphans_to_check_dict = manager.dict()
	# identical_index = manager.Value('H', 0)
	# identical_index_lock = manager.Lock()

	global gene_to_rough_cluster
	with open(repo_path + "nodes/" + args.node + "/trees/gene_to_cluster.pkl", "r") as f:
		gene_to_rough_cluster = pickle.load(f)

	children_cons = {}
	pickleMaps = {}
	old_potentials = {}
	# load locus_mapping files from children
	for c in args.children:
		with open(args.node_dir + c + "/locus_mappings.pkl", 'r') as pklFile:
			pickleMaps[c] = pickle.load(pklFile)
		if args.synteny:
			with open(args.node_dir + c + "/synteny_data.pkl", 'r') as pklFile:
				synteny_data[c] = pickle.load(pklFile)
		if c[0] == "L":
			with open(args.node_dir + c + "/" + c + ".pkl", "r") as f:
				children_cons[c] = pickle.load(f)
		else:  # c[0] == "N"
			with open(args.node_dir + c + "/consensus_data.pkl", "r") as f:
				children_cons[c] = pickle.load(f)
			with open(args.node_dir + c + "/potential_inparalogs.pkl", "r") as f:
				old_potentials.update(pickle.load(f))

	blast_pep = {}
	for c in args.children:
		my_blast_pep = open(my_dir + c + ".blast.fa", 'r').readlines()
		curBlast = ""
		curPep = ""
		for m in my_blast_pep:
			m = m.rstrip()
			if len(m) < 1:
				continue
			if m.find(">") > -1:
				if len(curBlast) > 0:
					blast_pep[curBlast].append(curPep)
					curPep = ""
				line = m[1:]
				curBlast = line.split(";")[0]
				if curBlast not in blast_pep:
					blast_pep[curBlast] = []
			else:
				curPep += m
		if len(curPep) > 0:
			blast_pep[curBlast].append(curPep)

	newPickleMap = {}  # this will turn into the locus mappings for this node
	if args.synteny:
		newSyntenyMap = {}
	newNewickMap = {"children": [set(args.children)]}
	childToCluster = {}  # child gene/og --> og.id for this node

	graphs = {}
	with open(repo_path + "nodes/" + args.node + "/trees/cluster_graphs.dat", "r") as f:
		to_parse = []
		clusterID = None
		for line in f:
			if line == "//\n":
				if to_parse:
					graphs[clusterID] = nx.parse_edgelist(to_parse, create_using=nx.DiGraph(), nodetype=str, data=True)
				to_parse = []
				clusterID = None
			else:
				if clusterID:
					to_parse.append(line.rstrip())
				else:
					clusterID = line.rstrip()

	# manager = multiprocessing.Manager()
	# cluster_counter = 1
	cluster_counter = Counter(1)
	# cluster_counter_lock = multiprocessing.Lock()

	ok_trees = []
	genes_to_cluster = {}  # not to mistake with gene_to_rough_cluster that contains rough clustering for synteny calculation
	with open(repo_path + "nodes/" + args.node + "/trees/orphan_genes.txt", "r") as f:
		for line in f:
			node = line.rstrip()
			new_orphan = "%s_%06d" % (mrca, cluster_counter.safeIncrement())
			# cluster_counter += 1
			ok_trees.append((new_orphan, (node,), (node, node)))  # (node,) comma is required so its a tuple that can be looped on and not on the string itself
			genes_to_cluster[node] = (new_orphan, False)

	ok_trees_list = [ok_trees]
	identical_orphans_to_check_list = [[]]  # 1st element empty for orphans
	identical_orphans_to_check_dict_list = [{}]  # 1st element empty for orphans
	potentials = {}
	cluster_queue = multiprocessing.JoinableQueue()
	result_queue = multiprocessing.Queue()

	refiners = [Refinery(cluster_queue, result_queue, mrca, cluster_counter) for i in xrange(args.numThreads)]
	for w in refiners:
		w.start()

	for cluster in graphs:
		cluster_queue.put(Refine(cluster, graphs[cluster]))

	for i in xrange(args.numThreads):
		cluster_queue.put(None)

	cluster_queue.join()

	for i in xrange(args.numThreads):
		(genes_to_cluster_2, ok_trees, identical_orphans_to_check, identical_orphans_to_check_dict, potentials_2) = result_queue.get()
		genes_to_cluster.update(genes_to_cluster_2)
		ok_trees_list.append(ok_trees)
		identical_orphans_to_check_list.append(identical_orphans_to_check)
		identical_orphans_to_check_dict_list.append(identical_orphans_to_check_dict)
		potentials.update(potentials_2)

	# manager.shutdown()

	in_paralogs = {}
	for old in old_potentials:
		if genes_to_cluster[old][1] is False:  # gene is still in an orphan cluster
			if genes_to_cluster[old][0] not in in_paralogs:
				in_paralogs[genes_to_cluster[old][0]] = [genes_to_cluster[g][0] for g in old_potentials[old]]
			else:
				in_paralogs[genes_to_cluster[old][0]].extend([genes_to_cluster[g][0] for g in old_potentials[old]])

	for (k, v) in in_paralogs.items():
		if k not in potentials:
			potentials[k] = v
		else:
			i = 0
			for o in v:
				if o not in potentials[k]:
					potentials[k].insert(i, o)  # use of counter to add in_paralogs at the begining of the list in the same order, because 1st inparalog should be the closest to leaves
					i += 1

	pickleSeqs = {}
	pickleToCons = {}
	cons_pkl = {}
	singletons = cluster_dir + "singletons.cons.pep"
	combined_orphans_headers = open(my_dir + args.node + ".combined_orphans_headers.txt", "w")
	combined_orphans_translation = {}
	singles = open(singletons, 'w')
	sum_stats = my_dir + "summary_stats.txt"
	sstats = open(sum_stats, 'w')

	combined_orphans_list = {}
	new_orphan_index = 1

	for nt in xrange(args.numThreads + 1):  # +1 for orphans
		# print "thread " + str(nt)
		# print "number of trees " + str(len(ok_trees_list[nt]))
		identical_orphans_processed_indexes = set()
		for ok in ok_trees_list[nt]:
			clusterID = ok[0]
			identical_orphans_to_check = identical_orphans_to_check_list[nt]
			identical_orphans_to_check_dict = identical_orphans_to_check_dict_list[nt]
			newPickleMap[clusterID] = []

			if args.synteny:
				newSyntenyMap[clusterID] = {'count': 0, 'neighbors': [], 'children': []}

			# get list of leaf sequences to pull and organize in treeSeqs
			treeSeqs = {}
			tree_seq_count = 0
			# leafSeqs = {}
			child_leaves = {}
			taxa = set([])
			taxa_map = {}
			for g in ok[1]:
				child = "_".join(g.split("_")[:-1])
				if args.synteny:
					newSyntenyMap[clusterID]['children'].append(g)
				childToCluster[g] = clusterID
				leafKids = pickleMaps[child][g]
				if child not in child_leaves:
					child_leaves[child] = 0
				child_leaves[child] += len(leafKids)
	# TODO in this part, children pickle are opened every time, could probably only open each one
				for l in leafKids:
					newPickleMap[clusterID].append(l)
					lKid = "_".join(l.split("_")[:-1])
					taxa.add(lKid)
					if lKid not in taxa_map:
						taxa_map[lKid] = 0
					taxa_map[lKid] += 1
					if lKid not in pickleSeqs:
						seqFile = args.node_dir + lKid + "/" + lKid + ".pkl"
						pklFile = open(seqFile, 'rb')
						pickleSeqs[lKid] = pickle.load(pklFile)
						pklFile.close()
					seq = pickleSeqs[lKid][l]
					if seq not in treeSeqs:
						treeSeqs[seq] = []
					treeSeqs[seq].append(l)
					tree_seq_count += 1
					if args.synteny:
						newSyntenyMap[clusterID]['count'] += 1
			newNewickMap[clusterID] = list(ok[2])

			my_lengths = []
			min_taxa = len(taxa)
			max_taxa = 0
			for tm in taxa_map:
				tm_val = taxa_map[tm]
				if tm_val > max_taxa:
					max_taxa = tm_val
				if tm_val < min_taxa:
					min_taxa = tm_val
			for seq in treeSeqs:
				for me in treeSeqs[seq]:
					my_lengths.append(len(seq))
			avg = numpy.average(my_lengths)
			std = numpy.std(my_lengths)
			std_avg = std / avg
			out_dat = [clusterID, str(len(my_lengths)), str(len(taxa)), str(min_taxa), str(max_taxa), str(min(my_lengths)), str(max(my_lengths)), str(avg), str(std), str(std_avg)]

	# 		sys.exit()
			sstats.write("\t".join(out_dat) + "\n")

			if len(ok[1]) == 1:
				g = ok[1][0]
				skip = False
				unique = False

				orphan_index = None
				if clusterID in identical_orphans_to_check_dict:
					if identical_orphans_to_check_dict[clusterID] not in identical_orphans_processed_indexes:
						orphan_index = identical_orphans_to_check_dict[clusterID]
						# print "reading orphan index " + str(orphan_index) + " out of max index " + str(len(identical_orphans_to_check) - 1) + " at process " + str(nt)
						unique = True
						combined_orphans_list = []
						for identical_orphan_tuple in identical_orphans_to_check[orphan_index]:
							child = "_".join(identical_orphan_tuple[0].split("_")[:-1])
							if children_cons[child][identical_orphan_tuple[0]].count(">") > 1:	 # identical orphan is tuple (child_node, parent_node)
								unique = False
								# print "not unique"
								del combined_orphans_list
								break
							else:
								# print "appending " + str(identical_orphan_tuple)
								combined_orphans_list.append(identical_orphan_tuple[1])
					else:
						# print "duplicate"
						skip = True
					# print "unique " + str(unique) + "; skip " + str(skip)

				child = "_".join(ok[1][0].split("_")[:-1])
	# 			seq = children_cons[child][ok[0][0]]
				seqs = {}
				if child[0] == "L":
					seqs[g] = children_cons[child][g]
				else:
					for seq in children_cons[child][g]:
						i = 0
						identifier = None
						for s in seq.rstrip().split("\n"):
							if not i % 2:
								identifier = s[1:].split(";")[0]
							else:
								seqs[identifier] = s
							i += 1
				# else:
					# need to get list of sequences to write them all
				for k, s in seqs.iteritems():
					cons_pkl[clusterID] = [">" + clusterID + ";" + str(len(s)) + "\n" + s + "\n"]
					# singletons_pep[clusterID] = [">" + clusterID + ";" + str(len(s)) + "\n" + s + "\n"]
					if not skip:
						if orphan_index is not None and unique:
							# print 'combining orphan index ' + str(orphan_index) + " out of max index " + str(len(identical_orphans_to_check) - 1) + " at process " + str(nt)
							combined_clusterID = "combined_" + str(new_orphan_index)  # use of a new index that is common between all multiprocesses as not to overlap
							new_orphan_index += 1
							size_str = ";" + str(len(s)) + "\n"
							singles.write(">" + combined_clusterID + size_str + s + "\n")
							combined_orphans_headers.write((size_str).join(combined_orphans_list) + size_str)
							combined_orphans_translation[combined_clusterID + ";" + str(len(s))] = [c + ";" + str(len(s)) for c in combined_orphans_list]
							identical_orphans_processed_indexes.add(orphan_index)
							# break # unique sequence so only one loop
						else:
							singles.write(">" + clusterID + ";" + str(len(s)) + "\n" + s + "\n")
			else:
				pickleToCons[clusterID] = []
				# newNewickMap[clusterID].append([])
				for g in ok[1]:
					child = "_".join(g.split("_")[:-1])
					seqs = {}
					if child[0] == "L":
						seqs[g] = children_cons[child][g]
					else:
						for seq in children_cons[child][g]:
							# if seq[0] == ">":  # else its a leaf so only sequence is present
							i = 0
							identifier = None
							for s in seq.rstrip().split("\n"):
								if not i % 2:  # not so that 0 is True
									identifier = s[1:].split(";")[0]
								else:
									seqs[identifier] = s
								i += 1
	# 						else:
	# 							seqs = [seq]
					i = 0
					for k, s in seqs.iteritems():
						pickleToCons[clusterID].append(">" + k + ";" + str(i) + ";" + str(len(s)) + "\n" + s + "\n")  # str(i) is a unique part in name so that all different names for muscle/fasttree
						# newNewickMap[clusterID][2].append(">" + k + ";" + str(len(s)) + "\n" + s + "\n")
						i += 1
			# cluster_counter += 1
	singles.close()
	sstats.close()
	combined_orphans_headers.close()

	with open(my_dir + "combined_orphans_translation_table.pkl", 'w') as f:
		pickle.dump(combined_orphans_translation, f)

	with open(my_dir + "pep_data.pkl", 'w') as f:
		pickle.dump(pickleToCons, f)

	# update synteny data
	if args.synteny:
		for clust in newSyntenyMap:
			for child in newSyntenyMap[clust]['children']:
				lc = "_".join(child.split("_")[:-1])
				# logger.debug("%s splitted to %s" % (child, lc))
				for neigh in synteny_data[lc][child]['neighbors']:
					# logger.debug("newSyntenyMap[%s]['neighbors'].append(childToCluster[%s]" % (clust, neigh))
					newSyntenyMap[clust]['neighbors'].append(childToCluster[neigh])
		# pickle synteny data
		pklSyn = my_dir + "synteny_data.pkl"
		sdat = open(pklSyn, 'wb')
		pickle.dump(newSyntenyMap, sdat)
		sdat.close()

	# pickle the locus mappings
	pklMap = my_dir + "locus_mappings.pkl"
	sdat = open(pklMap, 'wb')
	pickle.dump(newPickleMap, sdat)
	sdat.close()

	with open(my_dir + "clusters_newick.pkl", "w") as f:
		pickle.dump(newNewickMap, f)

	with open(my_dir + "pre_consensus_data.pkl", "w") as f:
		pickle.dump(cons_pkl, f)

	with open(my_dir + "potential_inparalogs.pkl", "w") as f:
		pickle.dump(potentials, f)

	with open(my_dir + "current_inparalogs.pkl", "w") as f:
		pickle.dump(in_paralogs, f)

	# script complete call
	clusters_done_file = my_dir + "CLUSTERS_REFINED"
	cr = open(clusters_done_file, 'w')
	cr.write("Way to go!\n")
	cr.close()


if __name__ == "__main__":
	main()
