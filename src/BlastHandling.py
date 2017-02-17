#!/usr/bin/env python

import operator
import logging
import networkx as nx
import pickle


class BlastSegment:
	# logger = logging.getLogger("BlastSegment")

	def __init__(self, query, target, pID, align_length, bitScore, evalue):
		self.query = query
		self.target = target
		# identity percentage
		self.pID = float(pID) / 100.0
		self.length = int(align_length)
		self.qLength = int(self.query.split(";")[1])
		self.tLength = int(self.target.split(";")[1])
		self.bitScore = float(bitScore)
		self.evalue = float(evalue)
		# percent = float(self.length) / float(min(self.qLength, self.tLength))  # mod for big BLAST
		percent = float(self.length) / float(max(self.qLength, self.tLength))
		percent = float(self.length) / float(self.qLength)
		self.adjPID = self.pID * percent
		self.score = (2.0 - self.adjPID) * 100000.0

	def getAdjPID(self):
		return self.adjPID

	def getScore(self):
		return self.score


class BlastParse:
	logger = logging.getLogger("BlastParse")
	EVALUE_THRESHOLD = 1e-4

	def __init__(self, m8_file, max_size_diff):
		self.m8_file = m8_file
		BlastParse.max_size_diff = max_size_diff

	# hits are scored by cumulative percent identity
	# synData is unused
	@staticmethod
	def scoreHits(hits, headers, min_best_hit, synData, minSynFrac):
		bestHits = nx.Graph()
		bestReciprocalHits = nx.Graph()
		head = open(headers, 'r').readlines()
		myHead = set([])
		for h in head:
			h = h.rstrip()
			line = h.split(";")[0]  # mod for big BLAST
			myHead.add(line)  # mod for big BLAST
		head = None
		bestHits.add_nodes_from(myHead)
		bestReciprocalHits.add_nodes_from(myHead)

		for q in hits:
			bestAdjPID = 0.0
			best_evalue = 1.0
			q_best = []
			q_hits = [hits[q][t] for t in hits[q]]
			q_hits.sort(key=operator.attrgetter('bitScore'), reverse=True)
			for ts in q_hits:
				ts_score = ts.getScore()
				t = ts.target.split(";")[0]
				if ts.evalue < float(BlastParse.EVALUE_THRESHOLD):
					if best_evalue == 1.0:  # and ts.evalue < float(1e-3):  # TODO change hardcoded evalue threshold
						bestAdjPID = ts.getAdjPID()
						best_evalue = ts.evalue
						q_best.append((q, t, ts_score))
					elif (ts.getAdjPID() > bestAdjPID * min_best_hit):  # and best_evalue < 1.0:
						q_best.append((q, t, ts_score))

			q_best = sorted(q_best, key=lambda tup: tup[2])
			for hit in q_best:
				if not bestHits.has_edge(hit[0], hit[1]):
					bestHits.add_edge(hit[0], hit[1], weight=hit[2], query="_".join(q.split("_")[:-1]))
				elif bestHits[hit[0]][hit[1]]["query"] != "_".join(q.split("_")[:-1]):  # not to add an edge in the reciprocal graph if there are simply multiple matches between same query and target
					bestReciprocalHits.add_edge(hit[0], hit[1], weight=hit[2])
		return bestReciprocalHits

	@staticmethod
	def makePutativeClusters(tree_dir, synData, bestReciprocalHits):
		# numThreads = 4
		# MAX_HITS = numHits
		BlastParse.logger.info("len(best hits nodes) %d %d" % (len(bestReciprocalHits.nodes()), len(bestReciprocalHits.edges())))
		subs = list(nx.connected_component_subgraphs(bestReciprocalHits))
		count = 1
		orphan_file = tree_dir + "orphan_genes.txt"
		orphans = open(orphan_file, 'w')
		BlastParse.logger.info("len(subs) = %s" % (len(subs)))
		# map child genes to rough clusters and vice versa
		geneToCluster = {}
		clusterToGenes = {}
		clusterToSub = {}
		gene_count = 0
		for s in subs:  # each subgraph is an initial cluster
			clusterID = "cluster_" + str(count)
			if len(s.nodes()) == 1:
				locus = s.nodes()[0]
				orphans.write(locus + "\n")
				geneToCluster[locus] = clusterID
				if clusterID not in clusterToGenes:
					clusterToGenes[clusterID] = []
				clusterToGenes[clusterID].append(locus)
			else:
				clusterToSub[clusterID] = s  # s is a graph object
				for locus in s.nodes():
					geneToCluster[locus] = clusterID
					if clusterID not in clusterToGenes:
						clusterToGenes[clusterID] = []
					clusterToGenes[clusterID].append(locus)
					# BlastParse.logger()
			count += 1
			gene_count += len(s.nodes())
		orphans.close()
		BlastParse.logger.info("gene count = %d" % (gene_count))
		BlastParse.logger.info("count = %d" % (count))

		with open(tree_dir + "gene_to_cluster.pkl", "w") as f:
			pickle.dump(geneToCluster, f)
		with open(tree_dir + "cluster_to_genes.pkl", "w") as f:
			pickle.dump(clusterToGenes, f)
		return 0

	# reads in the m8 file and returns hits, which is a dict of BlastSegments
	def readBlastM8(self):
		m8 = open(self.m8_file, 'r').readlines()
		hits = {}
		for m in m8:
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
			mySeg = BlastSegment(q, t, line[2], line[3], line[11], line[10])  # query,target,pID,length,bitScore,evalue
			if Q not in hits:
				hits[Q] = {}
			if T not in hits[Q]:
				hits[Q][T] = mySeg
			elif mySeg.bitScore > hits[Q][T].bitScore:  # Is this possible to find?
				hits[Q][T] = mySeg
		return hits
