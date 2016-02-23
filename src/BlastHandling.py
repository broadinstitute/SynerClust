#!/usr/bin/env python

import operator
import math
import logging
import networkx as nx


class BlastSegment:
	# logger = logging.getLogger("BlastSegment")

	def __init__(self, query, target, pID, align_length, bitScore, evalue):
		self.query = query
		self.target = target
		# identity percentage
		self.pID = float(pID) / 100.0
		self.length = int(align_length)
		# self.qStart = int(qStart)
		# self.tStart = int(tStart)
		# self.qStop = int(qStop)
		# self.tStop = int(tStop)
		# self.qLength = self.qStop-self.qStart +1
		# self.tLength = self.tStop-self.tStart +1
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
		percent = float(self.length) / float(min(self.qLength, self.tLength))  # mod for big BLAST
		self.adjPID = self.pID * percent
		self.score = (2.0 - self.adjPID) * 100000.0

	def getAdjPID(self):
		return self.adjPID

	def getScore(self):
		return self.score


class BlastParse:
	logger = logging.getLogger("BlastParse")

	def __init__(self, m8_file):
		self.m8_file = m8_file

	# hits are scored by cumulative percent identity
	# synData is unused
	@staticmethod
	def scoreHits(hits, headers, min_best_hit, synData, numHits, minSynFrac):
		# SYN_FRAC = minSynFrac
		# NUM_HITS = numHits
		bestHits = nx.Graph()
		bestDirHits = nx.DiGraph()
		head = open(headers, 'r').readlines()
		myHead = set([])
		for h in head:
			h = h.rstrip()
			line = h.split(";")[0]  # mod for big BLAST
			myHead.add(line)  # mod for big BLAST
		head = None
		bestHits.add_nodes_from(myHead)
		bestDirHits.add_nodes_from(myHead)

		for q in hits:
			q_node = "_".join(q.split("_")[:-1])
			q_hits_species = {}
			for t in hits[q]:
				t_spec = "_".join(t.split("_")[:-1])
				if t_spec not in q_hits_species:
					q_hits_species[t_spec] = []
				q_hits_species[t_spec].append(hits[q][t])
			for species in q_hits_species:
				q_hits = q_hits_species[species]
				q_best = []
				qd_best = []
				bestAdjPID = 0.0
				best_evalue = 1.0
				q_hits.sort(key=operator.attrgetter('bitScore'), reverse=True)
				# q_best_hash = {}
				for ts in q_hits:
					ts_score = ts.getScore()
					t = ts.target.split(";")[0]
					if best_evalue == 1.0 and ts.evalue < float(1e-3):  # TODO change hardcoded evalue threshold
						bestAdjPID = ts.getAdjPID()
						best_evalue = ts.evalue
						qd_best.append((q, t, ts_score))
						q_best.append((q, t, ts_score))
					elif (ts.getAdjPID() > bestAdjPID * min_best_hit) and best_evalue < 1.0:
						# TODO find out why using these hard coded values, and maybe change them?
						if (ts.evalue < float(1.0e-150)) or (best_evalue > 0.0 and (math.log10(best_evalue) + 30.0 > math.log10(ts.evalue))):
							qd_best.append((q, t, ts_score))
							q_best.append((q, t, ts_score))
				q_best = sorted(q_best, key=lambda tup: tup[2])
				qd_best = sorted(qd_best, key=lambda tup: tup[2])
				i = 0
				qbi = 1
				qdbi = 1
				while i < len(q_best) or i < len(qd_best):
					if (i < len(q_best)) and q_node != species:  # query species != target species, but why only in the not directed graph?
						if not bestHits.has_edge(q_best[i][0], q_best[i][1]):  # why only verify in the not directed graph?
							bestHits.add_edge(q_best[i][0], q_best[i][1], weight=q_best[i][2], rank=qbi)
							qbi += 1
					if i < len(qd_best):
						bestDirHits.add_edge(qd_best[i][0], qd_best[i][1], weight=qd_best[i][2], rank=qdbi)
						# print qd_best[i]
						qdbi += 1
					i += 1
		return (bestHits, bestDirHits)

	@staticmethod
	def makePutativeClusters(bestHits, tree_dir, synData, homScale, synScale, bestDirHits, numHits):
		# numThreads = 4
		# MAX_HITS = numHits
		BlastParse.logger.info("len(best hits nodes) %d %d" % (len(bestHits.nodes()), len(bestHits.edges())))
		subs = list(nx.connected_component_subgraphs(bestHits))
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
			# why are the subgraphs defined from the undirected graph and then used on the directed graph? Isn't RBH used?
			my_sub = bestDirHits.subgraph(s.nodes())
			bd_sub = my_sub.copy()
			# my_edges = [e for e in bd_sub.edges_iter()]
			my_edges = bd_sub.edges()
			edges_to_remove = []
			for med in my_edges:
				BlastParse.logger.debug("med[0][:3] == med[1][:3]\n\t\tmed[0] = %s\n\t\tmed[1] = %s" % (med[0], med[1]))
				# remove self hits
				if "_".join(med[0].split("_")[:-1]) == "_".join(med[1].split("_")[:-1]):
					edges_to_remove.append(med)
# 				elif not bd_sub[med[0]][med[1]]['rank'] <= MAX_HITS:
# 					edges_to_remove.append(med)
			bd_sub.remove_edges_from(edges_to_remove)

			while len(bd_sub.nodes()) > 0:
				clusterID = "cluster_" + str(count)
				sccs = list(nx.strongly_connected_component_subgraphs(bd_sub))
				scc = sccs[0]  # TODO verify : always only 1 possible, even when clustering more than 2 genomes?
				if len(scc.nodes()) == 1:
					locus = scc.nodes()[0]
					orphans.write(locus + "\n")
					geneToCluster[locus] = clusterID
					if clusterID not in clusterToGenes:
						clusterToGenes[clusterID] = []
					clusterToGenes[clusterID].append(locus)
				else:
					clusterToSub[clusterID] = scc  # scc is a graph object
					for locus in scc.nodes():
						geneToCluster[locus] = clusterID
						if clusterID not in clusterToGenes:
							clusterToGenes[clusterID] = []
						clusterToGenes[clusterID].append(locus)
				for NO in scc.nodes():
					bd_sub.remove_node(NO)
				count += 1
				gene_count += len(scc.nodes())
		orphans.close()
		BlastParse.logger.info("gene count = %d" % (gene_count))
		BlastParse.logger.info("count = %d" % (count))

		homout = open(tree_dir + "homology_matrices.dat", "w")
		synout = open(tree_dir + "synteny_matrices.dat", "w")
		for sub in clusterToSub:
			s = clusterToSub[sub]  # s is a graph object
			if len(s.nodes()) == 1:
				# is this possible???? If only 1 node it shoudn't be in clusterToSub
				continue
			(h_dist, s_dist) = BlastParse.makeDistanceMatrix(s, bestDirHits, geneToCluster, clusterToGenes, synData, homScale, synScale)
			# homology distance matrix
			h_string = ""
			for hd in h_dist:
				dists = []
				for he in h_dist[hd]:
					dists.append(str(h_dist[hd][he]))
				line = hd + "\t" + "\t".join(dists) + "\n"
				h_string = h_string + line
			h_string += "//\n"
# 			mat_file = tree_dir + sub + ".hom.dist"
# 			matout = open(mat_file, 'w')
			homout.write(h_string)
# 			matout.close()

			# synteny distance matrix
			s_string = ""
			for sd in s_dist:
				dists = []
				for se in s_dist[sd]:
					dists.append(str(s_dist[sd][se]))
				line = sd + "\t" + "\t".join(dists) + "\n"
				s_string = s_string + line
			s_string += "//\n"
# 			mat_file = tree_dir + sub + ".syn.dist"
# 			matout = open(mat_file, 'w')
			synout.write(s_string)
# 			matout.close()
		homout.close()
		synout.close()
		return 0

	# creates a distance matrix based on blast hits, augment distances with syntenic fractions
	@staticmethod
	def makeDistanceMatrix(graph, bestDirHits, geneToCluster, clusterToGenes, synData, homScale, synScale):

		# big_dist = homScale*200000.0+synScale*200000.0
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
		for m in graph.nodes():
			mySynDist[m] = {}
			syn_m = set(syn[m])
			mSeqs = len(syn[m])
			for n in graph.nodes():
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
					mySynDist[m][n] -= 0.0  # no neighbors in common if someone has no neighbors
					continue
				all_neighbors = syn_m & set(syn[n])
				for a in all_neighbors:
					t_m = max(syn[m].count(a), 0)
					t_n = max(syn[n].count(a), 0)
					matches += min(t_m, t_n)
				synFrac = float(matches) / float(min(mSeqs, nSeqs))
				mySynDist[m][n] = ((2.0 - synFrac) * 100000.0)
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
			# Q_len = q.split(";")[1]
			t = line[1]
			T = t.split(";")[0]
			# T_len = t.split(";")[1]
			if q == t:  # self hit
				continue

			mySeg = BlastSegment(q, t, line[2], line[3], line[11], line[10])  # query,target,pID,length,bitScore,evalue
			mySeg.setAdjPID()
			if Q not in hits:
				hits[Q] = {}
			if T not in hits[Q]:
				hits[Q][T] = mySeg
			elif mySeg.bitScore > hits[Q][T].bitScore:  # Is this possible to find?
				# print "better bits!", mySeg.bitScore, hits[Q][T].bitScore, Q, T
				hits[Q][T] = mySeg
		return hits
