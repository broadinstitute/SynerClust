#!/usr/bin/env python

import sys, numpy, math, logging
from scipy.stats import poisson
import networkx as nx
import NetworkX_Extension as nxe

class NJTree:
	logger = logging.getLogger("NJTree")

	def __init__(self, newick_file, syn_mat_file, mrca, alpha, beta, gamma, gain, loss):
		#self.newick_file = newick_file
		self.distance_file = newick_file
		self.dMatrix = ""
		self.newick = ""
		self.graph = nx.Graph()
		self.syntenyGraph = nx.Graph()
		self.syntenyMatrix = {}
		self.syntenyIndex = None
		self.synteny_file = syn_mat_file
		self.mrca = mrca
		self.bigNode = ""
		self.alpha = 10.0 # TODO take into account the actual value given
		#~ self.alpha = float(alpha)
		self.beta = 0.01 # TODO take into account the actual value given
		self.gamma = float(gamma)
		self.gain = float(gain)
		self.loss = float(loss)
		self.OK = "false"
		self.rootedTree = None
		self.hom_shortest_paths = None
		self.syn_shortest_paths = None
		self.paths = None
		self.centroid = ""
		self.gl_map = {} #node -> gain/loss tuple
	
	def readSyntenyMatrix(self, valid_nodes):
		matrix_data = open(self.synteny_file,'r').readlines()
		index = {}
		rev_index = {}
# 		count = 0
		valid_index = set([])
		tcount = 0
		for m in matrix_data:
			gene = m.split()[0]
			if gene in valid_nodes:
				valid_index.add(tcount)
				index[gene]=tcount
				rev_index[tcount]=gene
			tcount+=1
		for m in matrix_data:
			m=m.rstrip()
			dat = m.split()
			gene = dat[0]
			if not gene in valid_nodes:
				continue
			self.syntenyMatrix[gene] = {}
			dists = dat[1:]
			for vi in valid_index:
				self.syntenyMatrix[gene][rev_index[vi]] = float(dists[vi])
		self.syntenyIndex = index
			
	def readDistanceMatrix(self):
		matrix_data=open(self.distance_file,'r').readlines()
		return matrix_data
	
	def buildGraphFromDistanceMatrix(self, matrix_data):
		matrix={}
		unadded_nodes = set([])
		initial_indices = []
		for m in matrix_data:
			g = m.split()[0]
			initial_indices.append(g)
		for m in matrix_data:
			m = m.rstrip()
			dat = m.split()
			gene = dat[0]
			unadded_nodes.add(gene)
			my_species = "_".join(gene.split("_")[:-1])
			self.graph.add_node(gene,species=my_species)
			dists = dat[1:]
			matrix[gene] = {}
			for i in range(len(matrix_data)):
				matrix[gene][initial_indices[i]]=float(dists[i])
		self.readSyntenyMatrix(unadded_nodes)
		while len(unadded_nodes)>2:
			Udists = {}
			synUdists = {}
			uan = len(unadded_nodes)
			uan_denom = float(uan)-2.0
			for n in unadded_nodes:
				u = 0.0
				synu = 0.0
				for m in unadded_nodes:
					if n == m:
						continue
					u += matrix[n][m]#/uan_denom
					synu += self.syntenyMatrix[n][m]#/uan_denom
				Udists[n] = u
				synUdists[n] = synu
			min_nm = 1000000000.0
			minp = []
			for n in unadded_nodes:
				for m in unadded_nodes:
					if n==m:
						continue
					nm = uan_denom*matrix[n][m] - Udists[n] - Udists[m]
					if nm<min_nm:
						min_nm = nm
						minp = [n,m]
			minp.sort()
			mp0_mp_dist = 0.5*matrix[minp[0]][minp[1]] + 0.5*(Udists[minp[0]]-Udists[minp[1]])/uan_denom
			syn_mp0_mp_dist = 0.5*self.syntenyMatrix[minp[0]][minp[1]] + 0.5*(synUdists[minp[0]]-synUdists[minp[1]])/uan_denom
			mp1_mp_dist = 0.5*matrix[minp[0]][minp[1]] + 0.5*(Udists[minp[1]]-Udists[minp[0]])/uan_denom
			syn_mp1_mp_dist = 0.5*self.syntenyMatrix[minp[0]][minp[1]] + 0.5*(synUdists[minp[1]]-synUdists[minp[0]])/uan_denom

						
			newNode = ";".join(minp)
			my_species = ""
			if self.graph.node[minp[0]]['species'] == self.graph.node[minp[1]]['species']:
				my_species = self.graph.node[minp[0]]['species']
			else:
				my_species = self.mrca
			self.graph.add_node(newNode,species=my_species)
			for m in minp:
				unadded_nodes.remove(m)
# 			new_hrow = []
# 			new_srow = []
			matrix[newNode] = {}
			self.syntenyMatrix[newNode] = {}
			for k in unadded_nodes:
				dik = matrix[minp[0]][k]
				djk = matrix[minp[1]][k]
				dij = matrix[minp[0]][minp[1]]
				new_dist = (dik+djk-dij)/2.0
				matrix[newNode][k]=new_dist
				matrix[k][newNode]=new_dist
				
				dik = self.syntenyMatrix[minp[0]][k]
				djk = self.syntenyMatrix[minp[1]][k]
				dij = self.syntenyMatrix[minp[0]][minp[1]]
				new_dist = (dik+djk-dij)/2.0
				self.syntenyMatrix[newNode][k] = new_dist
				self.syntenyMatrix[k][newNode] = new_dist
			unadded_nodes.add(newNode)
			self.graph.add_edge(minp[0], newNode, homology_weight = mp0_mp_dist, synteny_weight = syn_mp0_mp_dist)
			self.graph.add_edge(minp[1], newNode, homology_weight = mp1_mp_dist, synteny_weight = syn_mp1_mp_dist)
		# TODO this if should be out of the while loop
		if len(unadded_nodes)==2:
			minp = []
			for ua in unadded_nodes:
				minp.append(ua)
			self.graph.add_edge(minp[0], minp[1], homology_weight = matrix[minp[0]][minp[1]], synteny_weight = self.syntenyMatrix[minp[0]][minp[1]])
			unadded_nodes = set([])
			minp.sort()
			big_md = ";".join(minp)
			unadded_nodes.add(big_md)
		bigNode = unadded_nodes.pop()
		self.bigNode = bigNode
		return bigNode
	
	def parseNewick(self, listt):
		nodes = []
		extinct_nodes = []
		ln_count = 0
		p_count = 0
		tp_count = 0
		for l in listt:
			ln = l.rstrip()
			self.newick+=ln
			if ln.find("(") > -1:
				p_count += 1
				tp_count += 1
			if ln.find(")") > -1:
				p_count -= 1
				tp_count += 1
			if len(ln) > 1:
				if ln.find(";") > -1:
					ln = ln[0:-1]
				if ln[0].find(":") == -1:
					c_paren = ln[-1].find(")") > -1
					nodes.append( (ln_count,ln[0:-1],p_count,tp_count, c_paren))
					identifier = ln.split(":")[0]
					species = identifier[0:identifier.find("_")]
					self.graph.add_node(identifier, row=ln_count, species=species)
				else:
					c_paren = ln[-1].find(")") > -1
					extinct_nodes.append( (ln_count, ln[0:-1], p_count, tp_count,c_paren))
			ln_count +=1
		nodes = sorted(nodes, key=lambda tup: tup[0])
		extinct = sorted(extinct_nodes, key=lambda tup: tup[0])
		return (nodes, extinct)
	
	@staticmethod
	def toNewick(graph):
		up = [] #unprocessed
		leaf = []
		for n in graph.nodes():
			if len(graph[n]) > 1:
				up.append(n)
			else:
				leaf.append((n,graph.node[n]['species']))
		curNode = None
		last_string = ""
		if len(graph.nodes()) == 2:
			ew = str(graph[leaf[0][0]][leaf[1][0]]['homology_weight'])
			last_string = "("+leaf[0][0]+":"+ew+","+leaf[1][0]+":"+ew+")"
		while len(up) > 0:
			(curNode,e_count) = NJTree.calcMostEdgesToLeaves(up, leaf, graph)
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
						ew = graph[curNode][e]['homology_weight']
						text = e_text+":"+str(ew)
						leaves.append(text)
			#add newick text to curNode
			node_text = "("+",".join(leaves)+")"
			last_string = node_text
			graph.node[curNode]['child_newick'] = node_text
			#change curNode to leaf
			cn_i = up.index(curNode)
			up.pop(cn_i)
			leaf.append((curNode, graph.node[curNode]['species']))
		if len(leaf) == 2 and len(up)==0 and len(graph.nodes()) > 2:
			ew = str(graph[leaf[0][0]][leaf[1][0]]['homology_weight'])
			last_string = "("+graph.node[leaf[0][0]]['child_newick']+":"+ew+","+graph.node[leaf[1][0]]['child_newick']+":"+ew+")"
		last_string = last_string.replace("(","(\n")
		last_string = last_string.replace(",",",\n")
		last_string = last_string.replace(")",")\n")
		last_string = last_string.rstrip()
		return last_string+";"
		
	def toNewickFromRooted(self):
		up = [] #unprocessed
		leaf = []
		for n in self.rootedTree.nodes():
			if len(self.rootedTree[n]) > 1:
				up.append(n)
			else:
				leaf.append((n,self.rootedTree.node[n]['species']))
		curNode = None
		while len(up) > 0:
			curNode = NJTree.calcMostEdgesToLeaves(up,leaf,self.rootedTree)[0]
			leaves = []
			for e in self.rootedTree[curNode]:
				for l in leaf:
					if l[0] == e:
						e_i = leaf.index(l)
						leaf.pop(e_i)
						e_text = l[0]
						if 'child_newick' in self.rootedTree.node[e]:
							e_text = self.rootedTree.node[l[0]]['child_newick']
						ew = self.rootedTree[curNode][l[0]]['homology_weight']
						text = e_text+":"+str(ew)
						leaves.append(text)
			#add newick text to curNode
			node_text = "("+",".join(leaves)+")"
			self.rootedTree.node[curNode]['child_newick'] = node_text
			#change curNode to leaf
			cn_i = up.index(curNode)
			up.pop(cn_i)
			leaf.append((curNode,self.rootedTree.node[curNode]['species']))
		return self.rootedTree.node[curNode]['child_newick']+";"
	
	def buildGraph(self, nodes, extinct):
		while len(nodes) > 3:
			nodes = sorted(nodes, key=lambda tup: tup[0])
			md = NJTree.getMinDistancePair(nodes)
			#remove nodes from current list
			m1_i = nodes.index(md[0])
			nodes.pop(m1_i)
			m2_i = nodes.index(md[1])
			nodes.pop(m2_i)
			#find next extinct node
			row_avg = float(md[0][0]+md[1][0])/2.0
			max_row = float(max(md[0][0],md[1][0]))
		
			m0_id = md[0][1][0:md[0][1].find(":")]
			m1_id = md[1][1][0:md[1][1].find(":")]
			m_id = [m0_id, m1_id]
			m_id.sort()
			identifier = ";".join(m_id)
			e_node = None
			min_pcount = min(md[0][2],md[1][2])
			for e in extinct:
				if float(e[0]) > max_row and e[2] <= min_pcount:
					e_node = e
					break
			#remove appropriate extinct node from extinct list; create new node based on merged children IDs
			e_i = extinct.index(e_node)
			extinct.pop(e_i)
			nodes.append((row_avg,identifier+e[1],e[2],e[3],e[4]))
			sp = ""
			if self.graph.node[m0_id]['species'] == self.graph.node[m1_id]['species']:
				sp = self.graph.node[m0_id]['species']
			else:
				sp = self.mrca
			self.graph.add_node(identifier,row=row_avg,species=sp)
			#add edges from new node to children
			for n in md:
				edge_weight = n[1].split(":")[1]
				ew = float(edge_weight[0:edge_weight.find(".")+5])
				self.graph.add_edge(n[1].split(":")[0], identifier, homology_weight=ew)
		if len(nodes) == 2:
			edge_weight = nodes[0][1].split(":")[1]
			ew = float(edge_weight[0:edge_weight.find(".")+5])
			self.graph.add_edge(nodes[0][1].split(":")[0], nodes[1][1].split(":")[0], homology_weight=ew)
			self.bigNode = ";".join([nodes[0][1][0:nodes[0][1].find(":")],nodes[1][1][0:nodes[1][1].find(":")]])
			return self.bigNode
		
		elif len(nodes) < 4:
			#merge all nodes into one extinct node
			semicolons = []
			rows = []
			node_ids = []
			for n in nodes:
				myID = n[1][0:n[1].find(":")]
				mySC = myID.count(";")
				semicolons.append((myID,mySC,n[0]))
				rows.append(n[0])
				node_ids.append(myID)
			sc = sorted(semicolons, key=lambda tup: tup[1])
			node_ids.sort()
			newID = ";".join(node_ids)
			row_sum = numpy.mean(rows)
			newSpecies = self.mrca
			if self.graph.node[sc[0][0]]['species'] == self.graph.node[sc[1][0]]['species']:
				newSpecies = self.graph.node[sc[0][0]]['species']
			#add new node to graph
			self.graph.add_node(newID,row=row_sum,species=newSpecies)
			#add edges to new node
			for n in nodes:
				edge_weight = n[1].split(":")[1]
				ew = float(edge_weight[0:edge_weight.find(".")+5])
				self.graph.add_edge(n[1].split(":")[0], newID, homology_weight=ew)
			self.bigNode = newID
			return self.bigNode

	@staticmethod
	def getMinDistancePair(nodes):
		min_dist = nodes[-1][0] -nodes[0][0]
		i = len(nodes) -1
		while i >=1:
			j=i-1
			if j >= 0:
				if i == j:
					j -= 1
				min_row = min(i,j)
				if nodes[min_row][4] == True:
					j-=1
				elif nodes[i][4] == False and nodes[j][4] == False:
					j-=1
				else:
					if abs(nodes[i][0] - nodes[j][0]) < min_dist:
						min_dist = abs(nodes[i][0] - nodes[j][0])
					j -= 1
					if min_dist == 1:
						j=-1
						i=-1
			i -= 1
		md_pairs = []
		while len(md_pairs) ==0:
			max_rparen = 0
			for n in nodes:
				#~ for m in nodes:
				m_index = nodes.index(n)+1
				if m_index<len(nodes):
					m = nodes[m_index]
					if m[0] < n[0]:
						continue
					if abs(n[0]-m[0]) <= min_dist:
						if n[4] == True:
							continue
						if n[4] == False and m[4] == False:
							continue
						if m[2] == 0 and len(nodes) > 2:
							continue
						my_mrp = max(n[2], m[2])
						if my_mrp > max_rparen:
							max_rparen = my_mrp
						md_pairs.append((n,m))
			if len(md_pairs) == 0:
				min_dist+=1
		#get pairs that have the maximum rparen counts
		max_rparen_pairs = []
		for md in md_pairs:
			if md[0][2] == max_rparen or md[1][2] == max_rparen:
				max_rparen_pairs.append(md)
		max_rparen_pairs = sorted(max_rparen_pairs, key=lambda pair: pair[1][3], reverse=True)
		return max_rparen_pairs[0]
		
	def rootTree(self):
		"""Return Score, root edges, number of losses
		"""
		#for each edge in 'tree' graph, score the tree
		roots = []
		min_gl = len(self.graph.nodes())*2
		#store shortest path matrix - it is the same for everyone
		if len(self.graph.nodes())>100:
# 			big_e = 0.0
			big_combo = 0.0
			e_pair = None
			for e in self.graph.edges():
				if e[0].find(";")>-1 and e[1].find(";")>-1:
					my_semi = min(e[0].count(";"),e[1].count(";"))
					my_big_e = self.graph[e[0]][e[1]]['homology_weight']
					my_big_combo = float(my_semi)*my_big_e
					if my_big_combo > big_combo:
						e_pair = e
						big_combo = my_big_combo
			return (-1.0, e_pair,len(self.graph.nodes()))
				
		else:
			self.hom_shortest_paths = nxe.all_pairs_path_length(self.graph, 'homology_weight')
			self.paths = nx.shortest_path(self.graph,None,None)
			self.syn_shortest_paths = nxe.all_pairs_path_length(self.graph, 'synteny_weight')
			for e in self.graph.edges():
				(score, tree, gl_sum, loss) = self.scoreEdge(e, min_gl)
				if gl_sum<min_gl:
					#~ print "gl_sum", gl_sum
					min_gl = gl_sum
				self.graph[e[0]][e[1]]['root_score'] = score
				roots.append((score, e, tree, loss))
		roots = sorted(roots, key=lambda tup: tup[0], reverse=True)
		if len(roots) == 1 or roots[0][0] == roots[1][0]:
			self.OK = "true"
		self.rootedTree = roots[0][2]
		return (roots[0][0], roots[0][1], roots[0][3])
			
	def findMidpointEdge(self):
		sp = self.shortest_paths
		longest = 0.0
		long_pair = None
		pairs = set([])
		for m in sp:
			for n in sp:
				my_pair = (m,n)
				if my_pair in pairs:
					continue
				else:
					pairs.add(my_pair)
					my_inv_pair = (n,m)
					pairs.add(my_inv_pair)
				if m==n:
					continue
				if sp[m][n] > longest:
					longest = sp[m][n]
					long_pair = (m,n)
				elif long_pair == None:
					long_pair = (m,n)
		path = nx.shortest_path(self.graph,long_pair[0],long_pair[1])
		cur_node = path.pop()
		cur_length = 0.0
		root_edge = None
		mid = longest/2.0
		while path:
			nextNode = path.pop()
			cur_length+= self.graph[cur_node][nextNode]['homology_weight']
			if cur_length > mid:
				root_edge = (cur_node,nextNode)
				break
			cur_node = nextNode
		if root_edge == None:
			#this basically means all of the sequences are identical
			#~ print "no root edge, setting to 'longest' pair"
			root_edge = long_pair
		return root_edge
	
	def scoreEdge(self, e, min_gl):
		#get homology distances, calculate variance
		h_dists = self.getHomologyDistances(e)
		h_var = numpy.var(h_dists)
		s_dists = self.getSyntenyDistances(e)
		s_var = numpy.var(s_dists)
		#get gain/loss count
		my_gl = 0
		gain = 0
		loss = 0
		if len(self.graph.nodes())>100:
			(gain, loss, tree) = self.getGainLossCount(e, -1)
			# why not keep the results from getGainLossCount?
			gain = min_gl/2
			loss = min_gl/2
		else:
			(gain, loss, tree) = self.getGainLossCount(e, min_gl)
			my_gl = gain + loss
		#score root
		my_poisson = 0.0
# 		my_gain_poisson = 0.0
		my_gain_poisson = poisson.pmf((gain), self.gain)
		if my_gain_poisson>0:
			my_poisson += math.log10(my_gain_poisson)
# 		my_loss_poisson = 0.0
		my_loss_poisson = poisson.pmf((loss), self.loss)
		if my_loss_poisson > 0:
			my_poisson += math.log10(my_loss_poisson)
		gl_factor = self.gamma * my_poisson
		dist_factor = (-self.beta) * h_var
		if h_var > 0:
			dist_factor = (-self.beta) * math.log10(h_var)
		syn_factor = (-self.alpha) * s_var
		if s_var > 0:
			syn_factor = (-self.alpha) * math.log10(s_var)
		score = math.exp(gl_factor) * math.exp(dist_factor) * math.exp(syn_factor)
		# score = math.exp(gl_factor + dist_factor + syn_factor) should be the same
		# and since exp is a strictly increasing function, the ranking would remain the same without applying it
		return (score,tree, my_gl,loss)
				
	#returns a list of distances from the root edge midpoint to all leaf nodes
	def getHomologyDistances(self, e):
		#bigNode is a concatenation of all leaf nodes in this tree, separated by a ;
		nodes = self.bigNode.split(";")
		dists = []
		for n in nodes:
			dists.append(self.getHomInterNodeDist(n, e, 'homology_weight'))
		return dists
		
	def getSyntenyDistances(self, e):
		nodes = self.bigNode.split(";")
		dists = []
		for n in nodes:
			dists.append(self.getSynInterNodeDist(n, e, 'synteny_weight'))
		return dists
		
	#this calculates the distance from leaf to root-edge midpoint
	def getHomInterNodeDist(self,n,e,attr):
		#find distance to farthest node on potential root edge
		raw_dist = min(self.hom_shortest_paths[n][e[0]], self.hom_shortest_paths[n][e[1]])
		near_node = e[0]
		if self.hom_shortest_paths[n][e[1]] == raw_dist:
			near_node = e[1]
		raw_len = 0.0
		raw_path = self.paths[n][near_node]
		for i in range(len(raw_path)-1):
			raw_len+=self.graph[raw_path[i]][raw_path[i+1]]['homology_weight']
		#subract half of root edge length from raw_dist to get the distance from the root (edge mid-point) to this node, n
		edge_length = self.graph[e[0]][e[1]][attr]
		mid_edge = edge_length/2.0
		dist = raw_len + mid_edge
		return dist	
		
	#this calculates the distance from leaf to root-edge midpoint
	def getSynInterNodeDist(self,n,e,attr):
		#find distance to farthest node on potential root edge
		raw_dist = min(self.syn_shortest_paths[n][e[0]], self.syn_shortest_paths[n][e[1]])
		near_node = e[0]
		if self.syn_shortest_paths[n][e[1]]==raw_dist:
			near_node = e[1]
		raw_len = 0.0
		raw_path = self.paths[n][near_node]
		for i in range(len(raw_path)-1):
			raw_len+=self.graph[raw_path[i]][raw_path[i+1]]['synteny_weight']
		
		#subract half of root edge length from raw_dist to get the distance from the root (edge mid-point) to this node, n
		edge_length = self.graph[e[0]][e[1]][attr]
		mid_edge = edge_length/2.0
		dist = raw_len + mid_edge
		return dist
		
	def getGainLossCount(self, e, min_gl):
		# TODO verify how the species of a node is set to MRCA
		# returns 0 gain when there needs to be a duplication event for the tree to exist (the loss after is found)
		gain = 0
		loss = 0
		gl_total = 0
		tGraph = self.graph.copy()
		newWeight = tGraph[e[0]][e[1]]['homology_weight']/2.0
# 		newSpecies = ""
		newID = "root"
		tGraph.remove_edge(e[0], e[1])
		tGraph.add_node(newID, species=self.mrca)
		tGraph.add_edge(e[0], newID, homology_weight=newWeight)
		tGraph.add_edge(e[1], newID, homology_weight=newWeight)
		# TODO no modification to synteny weights??
		up = [] #up = unprocessed, length is number of edges
		leaf = [] #leaf nodes
		for n in tGraph.nodes():
			if len(tGraph[n]) > 1:
				up.append(n)
			else:
				leaf.append((n,tGraph.node[n]['species']))
		#get node with most edges to leaves to process
		#note: species of a node can change during this process based on rooting, the e_leaf situation is for handling that
		while len(up) > 0:
			curNode = (NJTree.calcMostEdgesToLeaves(up, leaf, tGraph))[0] # curNode = AA node instead of root?
			curNodeSpecies = ""
			if curNode in self.gl_map:
				gain += self.gl_map[curNode]['gain']
				loss += self.gl_map[curNode]['loss']
				curNodeSpecies = self.gl_map[curNode]['species']
				gl_total = gain+loss
			else:
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
						pass # useless?
					else:
						curNodeSpecies = child[0][1]
						gain+=1
						gl_total+=1
				else:
					curNodeSpecies = self.mrca
					#2 child species
					if self.mrca in childSpecies:
						loss+=1
						gl_total+=1
				#~ self.gl_map[curNode] = {'gain':gain,'loss':loss, 'species':curNodeSpecies}
			cn_i = up.index(curNode)
			up.pop(cn_i)
			leaf.append((curNode,curNodeSpecies))
			tGraph.node[curNode]['species'] = curNodeSpecies
		return (gain, loss,tGraph)
	
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
			if e_count ==2:
				return (retNode,mostLeaves)
		return (retNode,mostLeaves)

	#checks if tree needs to be split
	def checkTree(self, root, NO_BREAK_EW):
		""" Returns "true" if they are 2 species in the tree and the mrca is not present.
		Returns "false" if they are 2 species in the tree but one of them is the mrca, or if more than 100 nodes in the graph.
		Returns "orphan" if they is only 1 species in the tree
		"""
		if len(self.graph.nodes())>100: # TODO why this arbitrary limit on the size of the tree?
			self.OK = "false"
			return self.OK
		#check species of root node neighbors
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
		if len(species) == 1 or len(set_species) == 1 and (not self.mrca in set_species):
			#all leaf nodes get added to orphan pile
			self.OK = "orphan"
			return self.OK
		if len(species) == 2:
			if self.mrca in species:
				self.OK = "false"
				return self.OK
			else:
				self.OK = "true"
				return self.OK
		# How can the MRCA be one of the species??
		# should add a verification that there are no more than 2 species

	#Centroid is based on distance. 
	def calcCentroid(self):
		if self.shortest_paths == None:
			self.shortest_paths = nx.shortest_path_length(self.graph,None,None,'homology_weight')
		nodes = self.bigNode.split(";")
		minDist = float(len(nodes)*2)
		minRoot = ""
		for r in nodes:
			dist = 0.0
			for n in nodes:
				dist+=self.shortest_paths[r][n]
			if dist < minDist:
				minDist = dist
				minRoot = r
		self.centroid = minRoot
	
	def splitTree(self, root):
		self.graph.remove_edge(root[1][0],root[1][1])
		new_graphs = nx.connected_component_subgraphs(self.graph)
		newicks = []
		matrices = []
		for n in new_graphs:
			#remove "root" node, put an edge in its place... these trees shouldn't be rooted! 
			#then get the newick format of the unrooted tree
			myRoot = None
			if n.has_node(root[1][1]):
				myRoot = root[1][1]
			else:
				myRoot = root[1][0]
			children = []
			if not len(n[myRoot]) ==2:
				if len(n.nodes()) == 1:
					newicks.append(("("+myRoot+");"))
					matrices.append((myRoot+"\t0.0"))
				else:
					sys.exit("splitTree: not 2 edges from temp root!:  "+str(len(n[myRoot])))
			else:
				for e in n[myRoot]:
					ew = n[myRoot][e]['homology_weight']
					children.append((e,ew))
				n.remove_node(myRoot)
				n.add_edge(children[0][0],children[1][0],homology_weight=(children[0][1]+children[1][1]))
				matrices.append(NJTree.makeDistanceMatrix(myRoot, n))
				newicks.append(NJTree.toNewick(n))
		return (newicks,matrices)
	
	@staticmethod
	def makeDistanceMatrix(root, N):
		d = nx.shortest_path(N)
		leaf = []
		for nod in N.nodes():
			if nod.count(";")==0:
				leaf.append(nod)
		mat_string = ""
		for m in leaf:
			if len(mat_string)>0:
				mat_string+="\n"
			m_dists = []
			m_dists.append(m)
			for n in leaf:
				mn_dist = 0.0
				if not m==n:
					path = d[m][n]
					while len(path)>1:
						mn_dist+= N[path[0]][path[1]]['homology_weight']
						path.pop(0)
				m_dists.append(str(mn_dist))
			mat_string += "\t".join(m_dists)
		return mat_string

	def splitBigTree(self):
		g_edges = self.graph.edges(data=True)
		biggest_edge_weight = 0.0
		biggest_edge = (None, None)
		for ge in g_edges:
			gw = ge[2]['homology_weight']
			if gw > biggest_edge_weight:
				biggest_edge_weight = gw
				biggest_edge = (ge[0],ge[1])
		self.graph.remove_edge(biggest_edge[0],biggest_edge[1])
		unchecked_new_graphs = nx.connected_component_subgraphs(self.graph)
		too_big = 1
		new_graphs = []
		while too_big==1:
			ok = 1
			while len(unchecked_new_graphs) > 0:
				ng = unchecked_new_graphs.pop()
				if len(ng.nodes())>50:
					ok=0
					g_edges = ng.edges(data=True)
					biggest_edge_weight = 0.0
					biggest_edge = (None, None)
					for ge in g_edges:
						gw = ge[2]['homology_weight']
						if gw > biggest_edge_weight:
							biggest_edge_weight = gw
							biggest_edge = (ge[0],ge[1])
					ng.remove_edge(biggest_edge[0],biggest_edge[1])
					my_subs = nx.connected_component_subgraphs(ng)
					for ms in my_subs:
						unchecked_new_graphs.append(ms)
				else:
					new_graphs.append(ng)

			if ok ==1:
				too_big=0
		newicks = []
		for n in new_graphs:
			#remove "root" node, put an edge in its place... these trees shouldn't be rooted! 
			#then get the newick format of the unrooted tree
			myRoot = None
			if n.has_node(biggest_edge[0]):
				myRoot = biggest_edge[0]
			else:
				myRoot = biggest_edge[1]
# 			children = []
			if not len(n[myRoot]) ==2:
				if len(n.nodes()) == 1:
					newicks.append(("("+myRoot+");"))
				else:
					print n[myRoot]
					sys.exit("splitTree: not 2 edges from temp root!:  "+str(len(n[myRoot])))
			else:
				newicks.append(NJTree.toNewick(n))
		return newicks
