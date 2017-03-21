#!/usr/bin/env python

import numpy
import math
import logging
from scipy.stats import poisson
import networkx as nx
import NetworkX_Extension as nxe
# from GSA import Edge


class NJTree:
	logger = logging.getLogger("NJTree")

	def __init__(self, mrca, alpha, beta, gamma, gain, loss, synteny):
		self.graph = nx.Graph()
		self.bigNode = ""
		self.alpha = float(alpha)
		self.beta = float(beta)
		self.gamma = float(gamma)
		self.gain = float(gain)
		self.loss = float(loss)
		self.synteny = synteny
		self.OK = "false"
		self.rootEdge = None
		self.mrca = mrca
		self.rootedTree = None
		self.hom_shortest_paths = None
		self.syn_shortest_paths = None
		self.paths = None
		# self.gl_map = {}  # node -> gain/loss tuple

	def readDistanceMatrix(self):
		return self.distance_matrix

	def buildGraphFromNewDistanceMatrix(self, hom_matrix, syn_matrix, leaves):
		for l in leaves:
			my_species = "_".join(l.split("_")[:-1])
			self.graph.add_node(l, species=my_species)
# TODO Verify order of the leaves names and of the data

		unadded_nodes = leaves
		unadded_count = len(unadded_nodes)
		while unadded_count > 2:
			uan = unadded_count
			uan_denom = float(uan) - 2.0
			matrix_size = unadded_count * (unadded_count - 1) / 2
			sum_of_hom = numpy.zeros(unadded_count, float)
			sum_of_syn = numpy.zeros(unadded_count, float)
			imax = 1
			pos = 0
			while(pos < matrix_size):
				for i in xrange(imax):
					sum_of_hom[i] += hom_matrix[pos]
					sum_of_syn[i] += syn_matrix[pos]
					sum_of_hom[imax] += hom_matrix[pos]
					sum_of_syn[imax] += syn_matrix[pos]
					pos += 1
				imax += 1
			min_nm = float("Inf")
			minp = []
			k = 0
			l = 1
			for i in xrange(matrix_size):
				nm = uan_denom * (hom_matrix[i] + syn_matrix[i]) - (sum_of_hom[k] + sum_of_syn[k]) - (sum_of_hom[l] + sum_of_syn[l])  # need to look back for these indices
				if nm < min_nm:
					min_nm = nm
					minp = [k, l]
				k += 1
				if not k < l:
					k = 0
					l += 1
# 			mp0_mp_dist = 0.5 * hom_matrix[minp[0] + ((minp[1] - 1) * minp[1] / 2)] + 0.5 * ((sum_of_hom[minp[0]] + sum_of_syn[minp[0]]) - (sum_of_hom[minp[1]] + sum_of_syn[minp[1]])) / uan_denom
			mp0_mp_dist = 0.5 * hom_matrix[minp[0] + ((minp[1] - 1) * minp[1] / 2)] + 0.5 * (sum_of_hom[minp[0]] - sum_of_hom[minp[1]]) / uan_denom
# 			syn0_mp_dist = 0.5 * syn_matrix[minp[0] + ((minp[1] - 1) * minp[1] / 2)] + 0.5 * ((sum_of_hom[minp[0]] + sum_of_syn[minp[0]]) - (sum_of_hom[minp[1]] + sum_of_syn[minp[1]])) / uan_denom
			syn0_mp_dist = 0.5 * syn_matrix[minp[0] + ((minp[1] - 1) * minp[1] / 2)] + 0.5 * (sum_of_syn[minp[0]] - sum_of_syn[minp[1]]) / uan_denom
# 			mp1_mp_dist = 0.5 * hom_matrix[minp[0] + ((minp[1] - 1) * minp[1] / 2)] + 0.5 * ((sum_of_hom[minp[1]] + sum_of_syn[minp[1]]) - (sum_of_hom[minp[0]] + sum_of_syn[minp[0]])) / uan_denom
			mp1_mp_dist = 0.5 * hom_matrix[minp[0] + ((minp[1] - 1) * minp[1] / 2)] + 0.5 * (sum_of_hom[minp[1]] - sum_of_hom[minp[0]]) / uan_denom
# 			syn1_mp_dist = 0.5 * syn_matrix[minp[0] + ((minp[1] - 1) * minp[1] / 2)] + 0.5 * ((sum_of_hom[minp[1]] + sum_of_syn[minp[1]]) - (sum_of_hom[minp[0]] + sum_of_syn[minp[0]])) / uan_denom
			syn1_mp_dist = 0.5 * syn_matrix[minp[0] + ((minp[1] - 1) * minp[1] / 2)] + 0.5 * (sum_of_syn[minp[1]] - sum_of_syn[minp[0]]) / uan_denom

			newNode = ";".join([unadded_nodes[minp[0]], unadded_nodes[minp[1]]])
			my_species = ""
			if self.graph.node[unadded_nodes[minp[0]]]['species'] == self.graph.node[unadded_nodes[minp[1]]]['species']:
				my_species = self.graph.node[unadded_nodes[minp[0]]]['species']
			else:
				my_species = self.mrca
			self.graph.add_node(newNode, species=my_species)
			# replace first merged leave by newNode then shift everything after the 2nd merged leave
			self.graph.add_edge(unadded_nodes[minp[0]], newNode, homology_dist=mp0_mp_dist, synteny_dist=syn0_mp_dist)
			self.graph.add_edge(unadded_nodes[minp[1]], newNode, homology_dist=mp1_mp_dist, synteny_dist=syn1_mp_dist)

			unadded_nodes[minp[0]] = newNode
			for i in xrange(minp[1], unadded_count - 1):
				unadded_nodes[i] = unadded_nodes[i + 1]
			unadded_count -= 1  # replaced 2 nodes with 1

		# replace the first line/column of the merging with the merging and shift values after second line/column
			k = 0
			l = 1
			offset = 0
			dfg_hom = hom_matrix[(minp[1] * (minp[1] - 1) / 2) + minp[0]]
			dfg_syn = syn_matrix[(minp[1] * (minp[1] - 1) / 2) + minp[0]]
			for pos in xrange(matrix_size):
				if k == minp[1] or l == minp[1]:
					offset += 1
				elif l == minp[0]:
					dfk_hom = hom_matrix[pos]
					dgk_hom = hom_matrix[(minp[1] * (minp[1] - 1) / 2) + k]
					dfk_syn = syn_matrix[pos]
					dgk_syn = syn_matrix[(minp[1] * (minp[1] - 1) / 2) + k]
					hom_matrix[pos] = 0.5 * (dfk_hom + dgk_hom - dfg_hom)
					syn_matrix[pos] = 0.5 * (dfk_syn + dgk_syn - dfg_syn)
				elif k == minp[0]:
					dfk_hom = hom_matrix[pos]
					dgk_hom = hom_matrix[pos + minp[1] - minp[0]]
					dfk_syn = syn_matrix[pos]
					dgk_syn = syn_matrix[pos + minp[1] - minp[0]]
					hom_matrix[pos - offset] = 0.5 * (dfk_hom + dgk_hom - dfg_hom)
					syn_matrix[pos - offset] = 0.5 * (dfk_syn + dgk_syn - dfg_syn)
				else:
					hom_matrix[pos - offset] = hom_matrix[pos]
					syn_matrix[pos - offset] = syn_matrix[pos]
				k += 1
				if not k < l:
					k = 0
					l += 1
		if unadded_count == 2:
			self.graph.add_edge(unadded_nodes[0], unadded_nodes[1], homology_dist=hom_matrix[0], synteny_dist=syn_matrix[0])  # check this
			unadded_nodes = [";".join(unadded_nodes[:2])]
		bigNode = unadded_nodes.pop()
		self.bigNode = bigNode
		return bigNode

	def getNewick(self):
		if self.rootedTree:
			processed = ['root']
			current_leaves = list(self.rootedTree['root'])
# 			nwk = "(" + ",".join(current_leaves) + ");"
# 			nwk = ",".join(current_leaves)
			nwk = "(" + current_leaves[0] + ":" + str(self.rootedTree['root'][current_leaves[0]]['homology_dist']) + ',' + current_leaves[1] + ":" + str(self.rootedTree['root'][current_leaves[1]]['homology_dist']) + ")"
			if self.synteny:
				nwk2 = "(" + current_leaves[0] + ":" + str(self.rootedTree['root'][current_leaves[0]]['synteny_dist']) + ',' + current_leaves[1] + ":" + str(self.rootedTree['root'][current_leaves[1]]['synteny_dist']) + ")"
			while current_leaves:
				n = current_leaves.pop()
				neighbors = list(self.rootedTree[n])
				if len(neighbors) > 1:  # if not a leaf
					for neighbor in neighbors:
						if neighbor in processed:
							neighbors.remove(neighbor)
							break
					processed.append(n)
# 					new_nwk = ",".join(neighbors)
					new_nwk = neighbors[0] + ":" + str(self.rootedTree[n][neighbors[0]]['homology_dist']) + ',' + neighbors[1] + ":" + str(self.rootedTree[n][neighbors[1]]['homology_dist'])
					nwk = nwk.replace(n, "(" + new_nwk + ")")
					if self.synteny:
						new_nwk2 = neighbors[0] + ":" + str(self.rootedTree[n][neighbors[0]]['synteny_dist']) + ',' + neighbors[1] + ":" + str(self.rootedTree[n][neighbors[1]]['synteny_dist'])
						nwk2 = nwk2.replace(n, "(" + new_nwk2 + ")")
					current_leaves.extend(neighbors)
			if self.synteny:
				return [nwk, nwk2]
			else:
				return [nwk]
		else:
			NJTree.logger.critical("Tried to get Newick from a tree that has no rootTree: %s" % (self.bigNode))

	@staticmethod
	def toNewick(graph):
		up = []  # unprocessed
		leaf = []
		for n in graph.nodes():
			if len(graph[n]) > 1:
				up.append(n)
			else:
				leaf.append((n, graph.node[n]['species']))
		curNode = None
		last_string = ""
		if len(graph.nodes()) == 2:
			ew = str(graph[leaf[0][0]][leaf[1][0]]['homology_dist'])
			last_string = "(" + leaf[0][0] + ":" + ew + "," + leaf[1][0] + ":" + ew + ")"
		while len(up) > 0:
			(curNode, e_count) = NJTree.calcMostEdgesToLeaves(up, leaf, graph)
			leaves = []
			for e in graph[curNode]:
				for l in leaf:
					if l[0] == e:
						e_i = leaf.index(l)
						e_text = e
						if 'child_newick' in graph.node[e]:
							if e_count > 2 and len(up) > 1:
								continue
							e_text = graph.node[e]['child_newick']
						leaf.pop(e_i)
						ew = graph[curNode][e]['homology_dist']
						text = e_text + ":" + str(ew)
						leaves.append(text)
			# add newick text to curNode
			node_text = "(" + ",".join(leaves) + ")"
			last_string = node_text
			graph.node[curNode]['child_newick'] = node_text
			# change curNode to leaf
			cn_i = up.index(curNode)
			up.pop(cn_i)
			leaf.append((curNode, graph.node[curNode]['species']))
		if len(leaf) == 2 and len(up) == 0 and len(graph.nodes()) > 2:
			ew = str(graph[leaf[0][0]][leaf[1][0]]['homology_dist'])
			last_string = "(" + graph.node[leaf[0][0]]['child_newick'] + ":" + ew + "," + graph.node[leaf[1][0]]['child_newick'] + ":" + ew + ")"
		last_string = last_string.replace("(", "(\n")
		last_string = last_string.replace(",", ",\n")
		last_string = last_string.replace(")", ")\n")
		last_string = last_string.rstrip()
		return last_string + ";"

	def rootTree(self):
		"""Return Score, root edge, number of losses
		"""
		# for each edge in 'tree' graph, score the tree
		roots = []
		min_gl = len(self.graph.nodes()) * 2

		if self.rootEdge is not None:  # check speed of this part in case of big cluster that's already been split once and is still big
			# self.hom_shortest_paths = nxe.all_pairs_path_length(self.graph, 'homology_dist')  # should already be present from splitNewTree if rooted because it was splitted from larger tree
			self.paths = nx.shortest_path(self.graph, None, None)
			# self.syn_shortest_paths = nxe.all_pairs_path_length(self.graph, 'synteny_dist')
			(score, tree, gl_sum, loss) = self.scoreEdge(self.rootEdge, min_gl)
			self.rootedTree = tree
			return (score, self.rootEdge, loss)

		if self.synteny:
			([self.hom_shortest_paths, self.syn_shortest_paths], self.paths) = nxe.all_pairs_path_length(self.graph, ['homology_dist', 'synteny_dist'])
		else:
			([self.hom_shortest_paths], self.paths) = nxe.all_pairs_path_length(self.graph, ['homology_dist'])
		# self.syn_shortest_paths = nxe.all_pairs_path_length(self.graph, 'synteny_dist')[0]
		# store shortest path matrix - it is the same for everyone
		if len(self.graph.nodes()) > 100:
			limit = len(self.graph.nodes()) / 2
			right_stack = [self.graph.nodes()[0]]
			left_stack = []
			degrees = {right_stack[0]: 1}
			while True:
				if right_stack:
					current_node = right_stack.pop()
				else:
					(current_node, to_degree) = left_stack.pop()
					neighbors = self.graph[to_degree].keys().remove(current_node)
					degrees[to_degree] = degrees[neighbors[0]] + degrees[neighbors[1]]
					if degrees[to_degree] >= limit:
						pair = neighbors[0] if degrees[neighbors[0]] >= degrees[neighbors[1]] else neighbors[1]
						return (-1.0, [to_degree, pair], len(self.graph.nodes()))
				right = False
				neighbors = self.graph[current_node].keys()
				if len(neighbors) == 1:
					degrees[current_node] = 1
				else:
					for neighbor in neighbors:
						if neighbor in degrees:
							continue
						if not right:
							right_stack.append(neighbor)
							right = True
						else:
							left_stack.append([neighbor, current_node])

			# # big_e = 0.0
			# big_combo = 0.0
			# e_pair = None
			# for e in self.graph.edges():
			# 	if e[0].find(";") > -1 and e[1].find(";") > -1:
			# 		my_semi = min(e[0].count(";"), e[1].count(";"))
			# 		my_big_e = self.graph[e[0]][e[1]]['homology_dist']
			# 		my_big_combo = float(my_semi) * my_big_e
			# 		if my_big_combo > big_combo:
			# 			e_pair = e
			# 			big_combo = my_big_combo
			# return (-1.0, e_pair, len(self.graph.nodes()))

		else:
			# self.paths = nx.shortest_path(self.graph, None, None)
			for e in self.graph.edges():
				(score, tree, gl_sum, loss) = self.scoreEdge(e, min_gl)
				if gl_sum < min_gl:
					# ~ print "gl_sum", gl_sum
					min_gl = gl_sum
				self.graph[e[0]][e[1]]['root_score'] = score
				roots.append((score, e, tree, loss))
		roots = sorted(roots, key=lambda tup: tup[0], reverse=True)
		if len(roots) == 1 or roots[0][0] == roots[1][0]:  # how does the second condition make the tree correct?
			self.OK = "true"
		self.rootedTree = roots[0][2]
		return (roots[0][0], roots[0][1], roots[0][3])

	def scoreEdge(self, e, min_gl):
		# get homology distances, calculate variance
		# h_dists = self.getHomologyDistances(e)
		h_dists = self.getDistances(e, 'homology_dist')
		h_var = numpy.var(h_dists)
		# s_dists = self.getSyntenyDistances(e)
		if self.synteny:
			s_dists = self.getDistances(e, 'synteny_dist')
			s_var = numpy.var(s_dists)
		else:
			s_var = 0  # syn_factor becomes a constant at -2 so doesn't change ranking
		# get gain/loss count
		my_gl = 0
		gain = 0
		loss = 0
		if len(self.graph.nodes()) > 100:
			(gain, loss, tree) = self.getGainLossCount(e, -1)
			# why not keep the results from getGainLossCount?
			gain = min_gl / 2
			loss = min_gl / 2
		else:
			(gain, loss, tree) = self.getGainLossCount(e, min_gl)
			my_gl = gain + loss
		# score root
		my_poisson = 0.0
# 		my_gain_poisson = 0.0
		my_gain_poisson = poisson.pmf((gain), self.gain)
		if my_gain_poisson > 0:
			my_poisson += math.log10(my_gain_poisson)
# 		my_loss_poisson = 0.0
		my_loss_poisson = poisson.pmf((loss), self.loss)
		if my_loss_poisson > 0:
			my_poisson += math.log10(my_loss_poisson)
		gl_factor = self.gamma * my_poisson
# 		dist_factor = (-self.beta) * h_var
		if h_var > 0:
			dist_factor = (-self.alpha) * math.log10(h_var)
		else:  # if var = 0, it means we're right in the middle
			dist_factor = (-self.alpha) * -2.0
# 		syn_factor = (-self.alpha) * s_var
		if s_var > 0:
			syn_factor = (-self.beta) * math.log10(s_var)
		else:
			syn_factor = (-self.beta) * -2.0
		score = math.exp(gl_factor) * math.exp(dist_factor) * math.exp(syn_factor)
		score2 = gl_factor + dist_factor + syn_factor
		# score = math.exp(gl_factor + dist_factor + syn_factor) should be the same
		# and since exp is a strictly increasing function, the ranking would remain the same without applying it
		return (score, tree, my_gl, loss)

	# returns a list of distances from the root edge midpoint to all leaf nodes
	def getDistances(self, e, attr):
		# bigNode is a concatenation of all leaf nodes in this tree, separated by a ;
		nodes = self.bigNode.split(";")
		dists = []
		for n in nodes:
			dists.append(self.getInterNodeDist(n, e, attr))
		return dists

	# def getHomologyDistances(self, e):
	# 	# bigNode is a concatenation of all leaf nodes in this tree, separated by a ;
	# 	nodes = self.bigNode.split(";")
	# 	dists = []
	# 	for n in nodes:
	# 		dists.append(self.getHomInterNodeDist(n, e, 'homology_dist'))
	# 	return dists

	# def getSyntenyDistances(self, e):
	# 	nodes = self.bigNode.split(";")
	# 	dists = []
	# 	for n in nodes:
	# 		dists.append(self.getSynInterNodeDist(n, e, 'synteny_dist'))
	# 	return dists

	# this calculates the distance from leaf to root-edge midpoint
	def getInterNodeDist(self, n, e, attr):
		dist = None
		if len(self.paths[e[0]][e[1]]) < len(self.paths[e[1]][e[0]]):
			dist = self.hom_shortest_paths[n][e[0]]
		else:
			dist = self.hom_shortest_paths[n][e[1]]
		if self.graph[e[0]][e[1]][attr] >= 0.0:
			return dist + (self.graph[e[0]][e[1]][attr] / 2)
		else:
			return dist

	# def getHomInterNodeDist(self, n, e, attr):
	# 	# find distance to farthest node on potential root edge
	# 	raw_dist = min(self.hom_shortest_paths[n][e[0]], self.hom_shortest_paths[n][e[1]])
	# 	near_node = e[0]
	# 	if self.hom_shortest_paths[n][e[1]] == raw_dist:
	# 		near_node = e[1]
	# 	raw_len = 0.0
	# 	raw_path = self.paths[n][near_node]
	# 	for i in range(len(raw_path) - 1):
	# 		raw_len += self.graph[raw_path[i]][raw_path[i + 1]]['homology_dist']
	# 	# subract half of root edge length from raw_dist to get the distance from the root (edge mid-point) to this node, n
	# 	edge_length = self.graph[e[0]][e[1]][attr]
	# 	mid_edge = edge_length / 2.0
	# 	dist = raw_len + mid_edge
	# 	return dist

	# # this calculates the distance from leaf to root-edge midpoint
	# def getSynInterNodeDist(self, n, e, attr):
	# 	# find distance to farthest node on potential root edge
	# 	raw_dist = min(self.syn_shortest_paths[n][e[0]], self.syn_shortest_paths[n][e[1]])
	# 	near_node = e[0]
	# 	if self.syn_shortest_paths[n][e[1]] == raw_dist:
	# 		near_node = e[1]
	# 	raw_len = 0.0
	# 	raw_path = self.paths[n][near_node]
	# 	for i in range(len(raw_path) - 1):
	# 		raw_len += self.graph[raw_path[i]][raw_path[i + 1]]['synteny_dist']

	# 	# subract half of root edge length from raw_dist to get the distance from the root (edge mid-point) to this node, n
	# 	edge_length = self.graph[e[0]][e[1]][attr]
	# 	mid_edge = edge_length / 2.0
	# 	dist = raw_len + mid_edge
	# 	return dist

	def getGainLossCount(self, e, min_gl):
		# TODO verify how the species of a node is set to MRCA
		# returns 0 gain when there needs to be a duplication event for the tree to exist (the loss after is found)
		gain = 0
		loss = 0
		gl_total = 0
		tGraph = self.graph.copy()
		newWeight = tGraph[e[0]][e[1]]['homology_dist'] / 2.0
		if self.synteny:
			newWeight2 = tGraph[e[0]][e[1]]['synteny_dist'] / 2.0
# 		newSpecies = ""
		newID = "root"
		tGraph.remove_edge(e[0], e[1])
		tGraph.add_node(newID, species=self.mrca)
		if self.synteny:
			tGraph.add_edge(e[0], newID, homology_dist=newWeight, synteny_dist=newWeight2)
			tGraph.add_edge(e[1], newID, homology_dist=newWeight, synteny_dist=newWeight2)
		else:
			tGraph.add_edge(e[0], newID, homology_dist=newWeight)
			tGraph.add_edge(e[1], newID, homology_dist=newWeight)
		# TODO no modification to synteny weights??
		up = []  # up = unprocessed, length is number of edges
		leaf = []  # leaf nodes
		for n in tGraph.nodes():
			if len(tGraph[n]) > 1:
				up.append(n)
			else:
				leaf.append((n, tGraph.node[n]['species']))
		# get node with most edges to leaves to process
		# note: species of a node can change during this process based on rooting, the e_leaf situation is for handling that
		while len(up) > 0:
			curNode = (NJTree.calcMostEdgesToLeaves(up, leaf, tGraph))[0]  # curNode = AA node instead of root?
			curNodeSpecies = ""
			# if curNode in self.gl_map:
			# 	gain += self.gl_map[curNode]['gain']
			# 	loss += self.gl_map[curNode]['loss']
			# 	curNodeSpecies = self.gl_map[curNode]['species']
			# 	gl_total = gain + loss
		# else:
			childSpecies = set([])
			child = []
			for e in tGraph[curNode]:
				e_leaf = None
				for l in leaf:
					if l[0] == e:
						e_leaf = l
				if e_leaf:
					e_i = leaf.index(e_leaf)
					leaf.pop(e_i)
					childSpecies.add(e_leaf[1])
					child.append(e_leaf)
			if len(childSpecies) == 1:
				if self.mrca in childSpecies:
					curNodeSpecies = self.mrca
					pass  # useless?
				else:
					curNodeSpecies = child[0][1]
					gain += 1
					gl_total += 1
			else:
				curNodeSpecies = self.mrca
				# 2 child species
				if self.mrca in childSpecies:
					# shouldn't there be a gain somewhere too in this case?
					loss += 1
					gl_total += 1
				# ~ self.gl_map[curNode] = {'gain':gain,'loss':loss, 'species':curNodeSpecies}
			cn_i = up.index(curNode)
			up.pop(cn_i)
			leaf.append((curNode, curNodeSpecies))
			tGraph.node[curNode]['species'] = curNodeSpecies
		return (gain, loss, tGraph)

	@staticmethod
	def calcMostEdgesToLeaves(unprocN, leaf, TG):
		"""
		unprocN = list of nodes that are not leaves?
		leaf = list of (node, species_name)
		TG = graph
		"""
		mostLeaves = -1
		retNode = None
		l_zero = []
		for l in leaf:
			l_zero.append(l[0])
		for n in unprocN:
			e_count = 0
			for e in TG[n]:
				if e in l_zero:
					e_count += 1
			if e_count > mostLeaves:
				mostLeaves = e_count
				retNode = n
			if e_count == 2:
				return (retNode, mostLeaves)
		return (retNode, mostLeaves)

	# checks if tree needs to be split
	def checkTree(self, root):
		""" Returns "true" if they are 2 species in the tree and the mrca is not present.
		Returns "false" if they are 2 species in the tree but one of them is the mrca, or if more than 100 nodes in the graph.
		Returns "orphan" if they is only 1 species in the tree
		mrca = current node
		"""
		if len(self.graph.nodes()) > 100:  # TODO why this arbitrary limit on the size of the tree?
			NJTree.logger.debug("Too big graph detected and arbitrarly split")
			self.OK = "false"
			return self.OK
		# check species of root node neighbors
		edge = root[1]
# 		loss = root[2]
		set_species = set([])
		species = []
# 		edge_weight = self.graph[edge[0]][edge[1]]['homology_weight']
		# TODO what is the point of having twice the results?
		set_species.add(self.rootedTree.node[edge[0]]['species'])
		set_species.add(self.rootedTree.node[edge[1]]['species'])
		species.append(self.rootedTree.node[edge[0]]['species'])
		species.append(self.rootedTree.node[edge[1]]['species'])
		if len(species) == 1 or len(set_species) == 1 and (self.mrca not in set_species):
			# all leaf nodes get added to orphan pile
			# NJTree.logger.debug("Orphan in checkTree still happens")
			self.OK = "orphan"
			return self.OK
		if len(species) == 2:
			if self.mrca in species:
				self.OK = "false"
				return self.OK
			else:
				self.OK = "true"
				return self.OK
		# should add a verification that there are no more than 2 species

	def splitNewTree(self, root):
		new_trees = []
		new_root_edges = []
		futur_roots = [root[1][0], root[1][1]]
		self.graph.remove_edge(root[1][0], root[1][1])
		new_graphs = nx.connected_component_subgraphs(self.graph)
		for n in new_graphs:
			for futur_root in futur_roots:  # loop here, but next if selects only 1 iteration
				if futur_root in n.nodes():
					new_tree = NJTree(self.mrca, self.alpha, self.beta, self.gamma, self.gain, self.loss, self.synteny)
					if len(n.nodes()) == 1:  # only 1 leaf on this half
						new_tree.bigNode = futur_root
						new_tree.rootEdge = (futur_root, futur_root)
						new_root_edges.append(new_tree.rootEdge)
						new_trees.append(new_tree)
						break
					new_hom_weight = 0.0
					new_syn_weight = 0.0
					children = []
					to_remove = []
					hom_attributes = nx.get_edge_attributes(n, 'homology_dist')
					if self.synteny:
						syn_attributes = nx.get_edge_attributes(n, 'synteny_dist')
					for child in nx.all_neighbors(n, futur_root):
						if child in futur_roots:
							continue  # other half of the graph
						children.append(child)
						try:
							new_hom_weight += hom_attributes[(child, futur_root)]
							if self.synteny:
								new_syn_weight += syn_attributes[(child, futur_root)]
						except KeyError:
							new_hom_weight += hom_attributes[(futur_root, child)]
							if self.synteny:
								new_syn_weight += syn_attributes[(futur_root, child)]
						to_remove.append([child, futur_root])
					for pair in to_remove:
						n.remove_edge(pair[0], pair[1])
					n.remove_node(futur_root)
					if self.synteny:
						n.add_edge(children[0], children[1], homology_dist=new_hom_weight, synteny_dist=new_syn_weight)
					else:
						n.add_edge(children[0], children[1], homology_dist=new_hom_weight)
					new_tree.rootEdge = (children[0], children[1])
					new_root_edges.append(new_tree.rootEdge)
					new_BigNode = ";".join([f for f in n.nodes() if len(n[f]) == 1])
					new_tree.bigNode = new_BigNode
					new_tree.graph = n
					new_tree.hom_shortest_paths = self.hom_shortest_paths
					if self.synteny:
						new_tree.syn_shortest_paths = self.syn_shortest_paths
					new_trees.append(new_tree)
					break
		return (new_trees, new_root_edges)
