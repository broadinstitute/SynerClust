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

	def __init__(self, query, target, pID, align_length, bitScore, evalue):
		self.query = query
		self.target = target
		# identity percentage
		self.pID = float(pID) / 100.0
		self.length = int(align_length)
		self.qLength = int(query.split(";")[1])
		self.bitScore = float(bitScore)
		self.evalue = float(evalue)
		percent = float(self.length) / float(self.qLength)
		self.adjPID = self.pID * percent

	def getAdjPID(self):
		return self.adjPID


class BlastParse:
	logger = logging.getLogger("BlastParse")
	node_path = None
	translation_table = None
	EVALUE_THRESHOLD = 1e-4
	to_add = {}
	# CORE_HITS_COUNT_THRESHOLD = 5
	# OVERLAP_PROPORTION_THRESHOLD = 0.5

	def __init__(self, max_size_diff, node_path, translation_table):
		BlastParse.max_size_diff = max_size_diff
		BlastParse.node_path = node_path
		BlastParse.translation_table = translation_table

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
			if ts.length == ts.qLength and ts.length == int(ts.target.split(";")[1]) and ts.pID == 1.0:  # identical
				# q_best.append((t, 1, 1.0, ts, 1))
				q_best.append((t, 1, 1.0, 1))
			elif ts.evalue < float(BlastParse.EVALUE_THRESHOLD):
				if best_evalue == 1.0:  # and ts.evalue < float(1e-3):  # TODO change hardcoded evalue threshold
					bestAdjPID = ts.getAdjPID()
					lastAdjPID = bestAdjPID
					best_evalue = ts.evalue
					# q_best.append((q, t, ts_score))
					# q_best.append((t, current_rank, 1.0, ts, 0))
					q_best.append((t, current_rank, 1.0, 0))
					# current_rank += 1
				elif (ts.getAdjPID() > bestAdjPID * min_best_hit):  # and best_evalue < 1.0:
					if ts.getAdjPID() != lastAdjPID:  # if not a tie
						current_rank += 1
					# q_best.append((t, current_rank, ts.getAdjPID() / bestAdjPID, ts, 0))
					q_best.append((t, current_rank, ts.getAdjPID() / bestAdjPID, 0))
					lastAdjPID = ts.getAdjPID()
		return q_best

	@staticmethod
	def prepareDiGraph(headers):
		bestReciprocalHits = nx.DiGraph()
		myHead = set([])
		with open(headers, 'r') as head:
			for h in head:
				h = h.rstrip()
				line = h.split(";")[0]  # mod for big BLAST
				myHead.add(line)  # mod for big BLAST
		bestReciprocalHits.add_nodes_from(myHead)  # need to add nodes this way and not only through edges to have orphans
		return bestReciprocalHits

	# hits are scored by cumulative percent identity
	# synData is unused
	@staticmethod
	def scoreHits(hits, bestReciprocalHits, min_best_hit, minSynFrac):
		for q in hits:
			q_hits = [hits[q][t] for t in hits[q]]
			q_hits.sort(key=operator.attrgetter('bitScore', 'evalue', 'adjPID', 'target'), reverse=True)
			q_best = BlastParse.getBestHits(q_hits, min_best_hit)

			q_best = sorted(q_best, key=lambda tup: tup[1])
			for hit in q_best:
				if (hit[0], q) not in BlastParse.to_add:
					BlastParse.to_add[(q, hit[0])] = (hit[1], hit[2], hit[3])
				else:
					# check if combined node to uncombine now before inserting all in the graph
					if q[:9] == "combined_":
						qss = BlastParse.translation_table[q]
					else:
						qss = [q]
					if hit[0][:9] == "combined_":
						tss = BlastParse.translation_table[hit[0]]
					else:
						tss = [hit[0]]
					reciprocal_edge = BlastParse.to_add.pop((hit[0], q))
					for qs in qss:
						for ts in tss:
							bestReciprocalHits.add_edge(qs, ts, rank=hit[1], m=hit[2], identity=hit[3])
							bestReciprocalHits.add_edge(ts, qs, rank=reciprocal_edge[0], m=reciprocal_edge[1], identity=reciprocal_edge[2])
		return bestReciprocalHits

	@staticmethod
	def makePutativeClusters(tree_dir, bestReciprocalHits):
		BlastParse.logger.info("len(best hits nodes) %d %d" % (len(bestReciprocalHits.nodes()), len(bestReciprocalHits.edges())))
		subs = list(nx.weakly_connected_component_subgraphs(bestReciprocalHits))  # only need weakly connected because the graph is built so that only reciprocal hits are part of the graph
		count = 1
		orphan_file = tree_dir + "orphan_genes.txt"
		orphans = open(orphan_file, 'w')
		BlastParse.logger.info("len(subs) = %s" % (len(subs)))
		# map child genes to rough clusters and vice versa
		geneToCluster = {}
		graphs = {}
		gene_count = 0

		for s in subs:  # each subgraph is an initial cluster
			s2 = [s]
			for s in s2:
				clusterID = "cluster_" + str(count)
				if len(s.nodes()) > 1:
					graphs[clusterID] = nx.generate_edgelist(s)  # data=True is default
				else:
					orphans.write(s.nodes()[0] + "\n")
				for locus in s.nodes():
					geneToCluster[locus] = clusterID
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
		return 0

	@staticmethod
	def readBlastM8FromFile(f,):
		data = open(f, "r").readlines()
		return BlastParse.readBlastM8(data)

	# reads in the m8 file and returns hits, which is a dict of BlastSegments
	@staticmethod
	def readBlastM8(data):
		hits = {}
		for m in data:
			m = m.rstrip()
			line = m.split("\t")
			if len(line) < 5:
				continue
			q = line[0]
			t = line[1]
			if q == t:  # self hit
				continue
			Q = q.split(";")[0]
			T = t.split(";")[0]
			if line[2] < 50.0 or line[3] < 0.5 * int(q.split(";")[1]):  # filter less than 50% identity and less than 50% of sequence length matches
				continue
			elif int(q.split(";")[1]) > BlastParse.max_size_diff * int(t.split(";")[1]) or int(t.split(";")[1]) > BlastParse.max_size_diff * int(q.split(";")[1]):  # size difference too big
				continue
			mySeg = BlastSegment(q, t, line[2], line[3], line[11], line[10])  # query,target,pID,length,bitScore,evalue
			if Q not in hits:
				hits[Q] = {}
			if T not in hits[Q]:
				hits[Q][T] = mySeg
			elif mySeg.bitScore > hits[Q][T].bitScore:  # Is this possible to find?
				hits[Q][T] = mySeg
		return hits
