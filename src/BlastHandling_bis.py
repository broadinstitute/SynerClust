#!/usr/bin/env python

import operator
import logging
import networkx as nx
import pickle
import os
# import subprocess

DEVNULL = open(os.devnull, 'w')
JACCARD_THRESHOLD = 0.2


class BlastSegment:
	# logger = logging.getLogger("BlastSegment")

	def __init__(self, query, target, pID, align_length, qstart, qend, bitScore, evalue):
		self.query = query
		self.target = target
		# identity percentage
		self.pID = float(pID) / 100.0
		self.length = int(align_length)
		self.qLength = int(self.query.split(";")[1])
		self.qstart = int(qstart)
		self.qend = int(qend)
		self.tLength = int(self.target.split(";")[1])
		self.bitScore = float(bitScore)
		self.evalue = float(evalue)
		# percent = float(self.length) / float(min(self.qLength, self.tLength))  # mod for big BLAST
		# percent = float(self.length) / float(max(self.qLength, self.tLength))
		percent = float(self.length) / float(self.qLength)
		self.adjPID = self.pID * percent
		# self.score = (2.0 - self.adjPID) * 100000.0

	def getAdjPID(self):
		return self.adjPID

	# def getScore(self):
	# 	return self.score


class BlastParse:
	logger = logging.getLogger("BlastParse")
	node_path = None
	EVALUE_THRESHOLD = 1e-4
	CORE_HITS_COUNT_THRESHOLD = 5
	OVERLAP_PROPORTION_THRESHOLD = 0.5

	def __init__(self, max_size_diff, node_path):
		BlastParse.max_size_diff = max_size_diff
		BlastParse.node_path = node_path

	@staticmethod
	def getBestHits(q_hits, min_best_hit):
		bestAdjPID = 0.0
		lastAdjPID = 0.0
		best_evalue = 1.0
		current_rank = 1
		q_best = []
		for ts in q_hits:
			# ts_score = ts.getScore()
			if ts.pID < 0.5 or ts.length < 0.5 * ts.qLength:
				continue
			t = ts.target.split(";")[0]
			if ts.evalue < float(BlastParse.EVALUE_THRESHOLD):
				if best_evalue == 1.0:  # and ts.evalue < float(1e-3):  # TODO change hardcoded evalue threshold
					bestAdjPID = ts.getAdjPID()
					lastAdjPID = bestAdjPID
					best_evalue = ts.evalue
					# q_best.append((q, t, ts_score))
					q_best.append((t, current_rank, 1.0, ts))
					# current_rank += 1
				elif (ts.getAdjPID() > bestAdjPID * min_best_hit):  # and best_evalue < 1.0:
					# q_best.append((q, t, ts_score))
					if ts.getAdjPID() != lastAdjPID:  # if not a tie
						current_rank += 1
					q_best.append((t, current_rank, ts.getAdjPID() / bestAdjPID, ts))
					lastAdjPID = ts.getAdjPID()
		return q_best

	# hits are scored by cumulative percent identity
	# synData is unused
	@staticmethod
	def scoreHits(hits, headers, min_best_hit, minSynFrac):
		# bestHits = nx.Graph()
		# bestReciprocalHits = nx.Graph()
		bestReciprocalHits = nx.DiGraph()
		to_add = {}
		head = open(headers, 'r').readlines()
		myHead = set([])
		for h in head:
			h = h.rstrip()
			line = h.split(";")[0]  # mod for big BLAST
			myHead.add(line)  # mod for big BLAST
		head = None
		# bestHits.add_nodes_from(myHead)
		bestReciprocalHits.add_nodes_from(myHead)  # need to add nodes this way and not only through edges to have orphans

		for q in hits:
			q_hits = [hits[q][t] for t in hits[q]]
			q_hits.sort(key=operator.attrgetter('bitScore'), reverse=True)
			q_best = BlastParse.getBestHits(q_hits, min_best_hit)

			# if q_best:
			# 	BlastParse.logger.debug("len(q_best) = " + str(len(q_best)) + " for " + q_best[0][2].query)
			# checking for protein domain increasing number of hits
			# if len(q_best) >= BlastParse.CORE_HITS_COUNT_THRESHOLD:
			# 	identifying the potential domain
			# 	(overlap_start, overlap_end, overlap_count) = BlastParse.longest_maximal_overlap_interval(q_best)
			# 	if overlap_count >= BlastParse.CORE_HITS_COUNT_THRESHOLD and (overlap_end - overlap_start < q_best[0][2].qLength * BlastParse.OVERLAP_PROPORTION_THRESHOLD):
			# 		BlastParse.logger.debug("Masking a region for query q = " + q)
			# 		target_child = q_best[0][0][:q_best[0][0].rfind("_")]
			# 		query_child = q[:q.rfind("_")]

			# 		# get query sequence
			# 		cmd = ["cat", BlastParse.node_path + query_child + ".blast.fa"]
			# 		cmd2 = ["grep", "-A", "1", q_best[0][2].query]
			# 		process = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=DEVNULL)
			# 		fasta_file = process.communicate()[0]
			# 		process = subprocess.Popen(cmd2, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=DEVNULL)
			# 		query_seq = process.communicate(fasta_file)[0].split("\n")

			# 		# masking domain
			# 		new_query = query_seq[0] + "\n" + query_seq[1][:overlap_start] + query_seq[1][overlap_start:overlap_end].lower() + query_seq[1][overlap_end:] + "\n"

			# 		# run blastp
			# 		db = BlastParse.node_path + target_child + ".blast.fa"
			# 		cmd = ["#BLAST_PATHblastp", "-outfmt", "6", "-evalue", "1", "-lcase_masking", "-db", db]
			# 		process = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=DEVNULL)
			# 		output = process.communicate(new_query)[0]

			# 		# parse output
			# 		new_q_hits = BlastParse.readBlastM8(output.split("\n"))
			# 		new_q_best = BlastParse.getBestHits([new_q_hits[q][t] for t in new_q_hits[new_q_hits.keys()[0]]], min_best_hit)

			# 		# combine new hits with original ones
			# 		filtered_q_best = []
			# 		for new_h in new_q_best:
			# 			for h in q_best:
			# 				if new_h[0] == h[0]:
			# 					filtered_q_best.append(h)
			# 					break
			# 		BlastParse.logger.debug("Query q = " + q + "\nPre-masking: " + str(len(q_best)) + " hits; Post-masking: " + str(len(filtered_q_best)) + " hits; Masked length = " + str(overlap_end - overlap_start) + " on " + str(overlap_count) + "sequences.")
			# 		### Verify whether there is at least 1 hit left
			# 		# assign new result to be used in the graph
			# 		q_best = filtered_q_best

			# q_best = sorted(q_best, key=lambda tup: tup[2])
			q_best = sorted(q_best, key=lambda tup: tup[1])
			for hit in q_best:
				# if not bestHits.has_edge(hit[0], hit[1]):
				# 	bestHits.add_edge(hit[0], hit[1], weight=hit[2], query="_".join(q.split("_")[:-1]))
				# elif bestHits[hit[0]][hit[1]]["query"] != "_".join(q.split("_")[:-1]):  # not to add an edge in the reciprocal graph if there are simply multiple matches between same query and target
					# bestReciprocalHits.add_edge(hit[0], hit[1], weight=hit[2])
				if (hit[0], q) not in to_add:
					to_add[(q, hit[0])] = (hit[1], hit[2])
				else:
					bestReciprocalHits.add_edge(q, hit[0], rank=hit[1], m=hit[2])
					reciprocal_edge = to_add.pop((hit[0], q))
					bestReciprocalHits.add_edge(hit[0], q, rank=reciprocal_edge[0], m=reciprocal_edge[1])

				# if not bestHits.has_edge(q, hit[0]):
				# 	bestHits.add_edge(q, hit[0], rank=hit[1], m=hit[2], query="_".join(q.split("_")[:-1]))
				# elif bestHits[q][hit[0]]["query"] != "_".join(q.split("_")[:-1]):  # not to add an edge in the reciprocal graph if there are simply multiple matches between same query and target
				# 	bestReciprocalHits.add_edge(q, hit[0], rank=hit[1], m=hit[2])
		return bestReciprocalHits

	@staticmethod
	def makePutativeClusters(tree_dir, bestReciprocalHits):
		# numThreads = 4
		# MAX_HITS = numHits
		BlastParse.logger.info("len(best hits nodes) %d %d" % (len(bestReciprocalHits.nodes()), len(bestReciprocalHits.edges())))
		# subs = list(nx.connected_component_subgraphs(bestReciprocalHits))
		subs = list(nx.weakly_connected_component_subgraphs(bestReciprocalHits))  # only need weakly connected because the graph is built so that only reciprocal hits are part of the graph
		count = 1
		orphan_file = tree_dir + "orphan_genes.txt"
		orphans = open(orphan_file, 'w')
		BlastParse.logger.info("len(subs) = %s" % (len(subs)))
		# map child genes to rough clusters and vice versa
		geneToCluster = {}
		graphs = {}
		# clusterToGenes = {}
		# clusterToSub = {}
		gene_count = 0
		for s in subs:  # each subgraph is an initial cluster
			# if len(s.nodes()) == 1:
				# orphans.write(s.nodes()[0] + "\n")
				# locus = s.nodes()[0]
				# geneToCluster[locus] = clusterID
				# if clusterID not in clusterToGenes:
				# 	clusterToGenes[clusterID] = []
				# clusterToGenes[clusterID].append(locus)
			# else:
				# clusterToSub[clusterID] = s  # s is a graph object
			# if len(s.nodes()) > BlastParse.CORE_HITS_COUNT_THRESHOLD:
			# 	s2 = remove_weak_links(s)
			# else:
			s2 = [s]
			for s in s2:
				clusterID = "cluster_" + str(count)
				if len(s.nodes()) > 1:
					graphs[clusterID] = nx.generate_edgelist(s)  # data=True is default
				else:
					orphans.write(s.nodes()[0] + "\n")
				for locus in s.nodes():
					geneToCluster[locus] = clusterID
				# 	if clusterID not in clusterToGenes:
				# 		clusterToGenes[clusterID] = []
				# 	clusterToGenes[clusterID].append(locus)
				count += 1
				gene_count += len(s.nodes())
		orphans.close()
		BlastParse.logger.info("gene count = %d" % (gene_count))
		BlastParse.logger.info("count = %d" % (count))

		with open(tree_dir + "gene_to_cluster.pkl", "w") as f:
			pickle.dump(geneToCluster, f)
		with open(tree_dir + "cluster_graphs.dat", "w") as f:
			for clusterID in graphs:
				f.write(clusterID + "\n")
				for line in graphs[clusterID]:
					f.write(line + "\n")
				f.write("//\n")
		# with open(tree_dir + "cluster_to_genes.pkl", "w") as f:
		# 	pickle.dump(clusterToGenes, f)
		return 0

	@staticmethod
	def longest_maximal_overlap_interval(hits):
		starts = []
		ends = []
		for hit in hits:
			starts.append(hit[2].qstart)
			ends.append(hit[2].qend)
		starts.sort()
		ends.sort()
		n = len(hits)
		i = j = 0
		maxi = maxj = -1
		maximal_overlap = 0
		current_overlap = 0
		while i < n and j < n:
			if starts[i] < ends[j]:
				current_overlap += 1
				if current_overlap > maximal_overlap:
					maximal_overlap = current_overlap
					maxi = i
					maxj = j
				elif current_overlap == maximal_overlap:
					if ends[j] - starts[i] > ends[maxj] - starts[maxi]:  # keep longest
						maxi = i
						maxj = j
				i += 1
			else:
				current_overlap -= 1
				j += 1
		return (maxi, maxj, maximal_overlap)

	@staticmethod
	def overlap_intervals(hits):
		starts = []
		ends = []
		for hit in hits:
			starts.append((hit[2].qstart, hit[2].target))
			ends.append((hit[2].qend, hit[2].target))
		starts.sort(key=lambda tup: tup[0])
		ends.sort(key=lambda tup: tup[0])
		n = len(hits)
		i = j = 0
		# maxi = maxj = -1
		# maximal_overlap = 0
		# current_overlap = 0
		intervals = []
		current_targets = set()
		last_position = 0
		while i < n and j < n:
			if starts[i][0] < ends[j][0]:
				# current_overlap += 1
				current_targets.add(starts[i][1])
				# if current_overlap > maximal_overlap:
				# 	maximal_overlap = current_overlap
				# 	maxi = i
				# 	maxj = j
				# elif current_overlap == maximal_overlap:
				# 	if ends[j] - starts[i] > ends[maxj] - starts[maxi]:  # keep longest
				# 		maxi = i
				# 		maxj = j
				intervals.append(last_position, starts[i], current_targets.copy())
				last_position = starts[i]
				i += 1
			else:
				# current_overlap -= 1
				current_targets.remove(ends[j][1])
				last_position = ends[j]
				j += 1
		return intervals

	@staticmethod
	def readBlastM8FromFile(f):
		data = open(f, "r").readlines()
		return BlastParse.readBlastM8(data)

	# reads in the m8 file and returns hits, which is a dict of BlastSegments
	@staticmethod
	def readBlastM8(data):
		hits = {}
		for m in data:
			m = m.rstrip()
			line = m.split()
			if len(line) < 5:
				continue
			q = line[0]
			Q = q.split(";")[0]
			t = line[1]
			T = t.split(";")[0]
			if q == t:  # self hit
				continue
			elif int(q.split(";")[1]) > BlastParse.max_size_diff * int(t.split(";")[1]) or int(t.split(";")[1]) > BlastParse.max_size_diff * int(q.split(";")[1]):  # size difference too big
				continue
			mySeg = BlastSegment(q, t, line[2], line[3], line[6], line[7], line[11], line[10])  # query,target,pID,length,qstart,qend,bitScore,evalue
			if Q not in hits:
				hits[Q] = {}
			if T not in hits[Q]:
				hits[Q][T] = mySeg
			elif mySeg.bitScore > hits[Q][T].bitScore:  # Is this possible to find?
				hits[Q][T] = mySeg
		return hits


def remove_weak_links(graph):
	nodes = graph.nodes()
	to_remove = []
	for i in xrange(len(nodes)):
		for j in graph[nodes[i]].keys():
			if j in nodes[:i + 1]:
				continue
			if jaccard_similarity(graph, nodes[i], j) < JACCARD_THRESHOLD:
				to_remove.append((nodes[i], j))
	for tr in to_remove:
		graph.remove_edge(tr[0], tr[1])
	return list(nx.weakly_connected_component_subgraphs(graph))


def jaccard_similarity(graph, n1, n2):
	s1 = set(graph[n1].keys())
	s1.remove(n2)
	s2 = set(graph[n2].keys())
	s2.remove(n1)
	return len(s1 & s2) / float(len(s1 | s2))
