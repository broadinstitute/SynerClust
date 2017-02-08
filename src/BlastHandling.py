#!/usr/bin/env python

import operator
import logging
import networkx as nx
import numpy
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
		self.adjPID = 0.0
		self.score = 0.0

	# returns (pID/100)*qLength
	# adjusted identity percentage, 100 000 < adjPID < 200 000
	# lower means more identity
	def setAdjPID(self):
		# percent of the shortest sequence that the blast result represents
		# percent = float(self.length) / float(min(self.qLength, self.tLength))  # mod for big BLAST
		# changed to having an adjustement based on the size of the biggest sequence because ~100aa protein could have high match to ~2000aa protein, which is probably an artifact
		percent = float(self.length) / float(max(self.qLength, self.tLength))
		self.adjPID = self.pID * percent
		self.score = (2.0 - self.adjPID) * 100000.0

	def getAdjPID(self):
		return self.adjPID

	def getScore(self):
		return self.score


class BlastParse:
	logger = logging.getLogger("BlastParse")
	EVALUE_THRESHOLD = 1e-4

	def __init__(self, m8_file):
		self.m8_file = m8_file

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
						# 	if (ts.evalue < float(1.0e-150)) or (best_evalue > 0.0 and (math.log10(best_evalue) + 30.0 > math.log10(ts.evalue))):
						# 	qd_best.append((q, t, ts_score))
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
			if 'L_0000000_vpfkc0x3-cBWYIYFcmcPqA_019882' in s.nodes():
				debug_buffer = ''
				for e in s.edges():
					debug_buffer += e[0] + " " + e[1] + "\n"
				BlastParse.logger.debug("edges:\n" + debug_buffer + "\n\n")
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

	# creates a distance matrix based on blast hits, augment distances with syntenic fractions
	@staticmethod
	def makeDistanceMatrix(graph, bestDirHits, geneToCluster, clusterToGenes, synData, homScale, synScale):
		big_dist = 200000.0
		myHomDist = {}
		# calculate blast distances
		for n in graph.nodes():
			myHomDist[n] = {}
			for m in graph.nodes():
				myHomDist[n][m] = big_dist
				if n == m:
					myHomDist[n][m] = 0.0
					continue
				if m in bestDirHits[n]:
					myHomDist[n][m] = bestDirHits[n][m]['weight']
					# myDist[n][m] = myDist[n][m] - (bestDirHits[g][h]['weight']*homScale)
					# myHomDist[n][m] = myHomDist[n][m] - bestDirHits[n][m]['weight']
		# populate neighbor lists with rough cluster IDs
		# synData contains one entry for each child node with their synteny_data.pkl
		syn = {}
		for d in myHomDist:
			# print "d", d
			syn[d] = []
			node = "_".join(d.split("_")[:-1])
			for n in synData[node][d]['neighbors']:
				syn[d].append(geneToCluster[n])
		# pairwise compare for syntenic fraction

		mySynDist = {}
		pairs = set([])
		all_nodes = graph.nodes()
		all_nodes.sort()
		syn_matrix = numpy.empty(len(all_nodes) * (len(all_nodes) - 1) / 2)
		i = 1
		pos = 0
		for m in all_nodes[1:]:
			mySynDist[m] = {}
			syn_m = set(syn[m])
			mSeqs = len(syn[m])
			for n in all_nodes[:i]:
				mySynDist[m][n] = big_dist
				my_pair = (m, n)
				if my_pair in pairs:
					mySynDist[m][n] = mySynDist[n][m]
					continue
				else:
					pairs.add(my_pair)
					my_inv_pair = (n, m)
					pairs.add(my_inv_pair)
				if n == m:
					mySynDist[m][n] = 0.0
					continue  # self comparisons should have identical neighbor sets anyway
				nSeqs = len(syn[n])
				matches = 0
				if mSeqs == 0 or nSeqs == 0:
					mySynDist[m][n] -= 0.0  # no neighbors in common if someone has no neighbors  # -= 0 ? does it change anything?
					continue
				all_neighbors = syn_m & set(syn[n])
				for a in all_neighbors:
					t_m = max(syn[m].count(a), 0)
					t_n = max(syn[n].count(a), 0)
					matches += min(t_m, t_n)
				synFrac = float(matches) / float(min(mSeqs, nSeqs))  # why mSeqs and not len(syn_m) which is a set that removes duplicates?
				mySynDist[m][n] = ((2.0 - synFrac) * 100000.0)
				syn_matrix[pos] = 1.0 - synFrac
				pos += 1
			i += 1
		return (myHomDist, mySynDist)

	# reads in the m8 file and returns hits, which is a dict of BlastSegments
	def readBlastM8(self):
		m8 = open(self.m8_file, 'r').readlines()
		hits = {}
		# old_t = ""
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

			mySeg = BlastSegment(q, t, line[2], line[3], line[11], line[10])  # query,target,pID,length,bitScore,evalue
			mySeg.setAdjPID()
			if Q not in hits:
				hits[Q] = {}
			if T not in hits[Q]:
				hits[Q][T] = mySeg
			elif mySeg.bitScore > hits[Q][T].bitScore:  # Is this possible to find?
				hits[Q][T] = mySeg
		return hits
