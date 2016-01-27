#!/usr/bin/env python

import sys, math, getopt, random, string, os, re, logging, hashlib, base64
import networkx as nx
import traceback

class Tree:
	logger = logging.getLogger("Tree")

	def __init__(self, tree_file, genomeToLocusFile):
		self.tree_file = tree_file
		self.genomeToLocusFile = genomeToLocusFile
		self.genomeToLocus = {}
		self.locusToGenome = {}
		self.tree_string = ""
		self.tree = nx.Graph()
		self.ancestral = []
		self.rooted_tree = nx.DiGraph()
		self.genomes = []
		
	def readTree(self):
		# tree = open(self.tree_file,'r').read()
		with open(self.tree_file) as t:
			tree = t.read()
			#~ print tree
			trees = re.split('\)\d+[\.\d+]*:',tree)
			jstring = "):"
			tree = jstring.join(trees)
			tree = tree.rstrip()
			print tree
			#~ sys.exit()
			self.tree_string = tree
		
	def readGenomeToLocusFile(self):
		tags = open(self.genomeToLocusFile,'r').readlines()
		for t in tags:
			t = t.rstrip()
			l = t.split()
			if not len(l) > 1:
				continue
			self.genomeToLocus[l[0]] = l[1]
			self.locusToGenome[l[1]] = l[0]
		
	def codeGenomeID(self, genome):
		Tree.logger.debug("".join(traceback.format_stack()))
		# Tree.logger.debug("\t"*len(traceback.format_stack()) + "codeGenomeID")
		tag = ''
		if genome in self.genomeToLocus:
			tag = self.genomeToLocus[genome]
			self.locusToGenome[tag] = genome
		else:
			children = genome.split(";")
			# max node degree between children nodes
			degree = max(int(children[0].split("_")[1]), int(children[1].split("_")[1])) + 1
			tag = "N_%07d_%s" %(degree, base64.urlsafe_b64encode(hashlib.md5(genome).digest())[:-2])
			# tag = ''.join(random.choice(string.ascii_uppercase) for x in range(3))
			# while tag in self.locusToGenome:
			# 	tag = ''.join(random.choice(string.ascii_uppercase) for x in range(3))
			self.genomeToLocus[genome] = tag
			self.locusToGenome[tag] = genome
			# print genome, tag
		return tag
		
	def reregisterGenomeID(self, id, newChildren):
		oldGenome = self.locusToGenome[id]
		newGenome = ";".join(newChildren)
		del self.genomeToLocus[oldGenome]
		self.genomeToLocus[newGenome] = id
		self.locusToGenome[id] = newGenome
		
	def parseTree(self):
		myTreeString = self.tree_string[0:-1]
		myExcisableRegions = self.indexParens(myTreeString)
		while len(myExcisableRegions) > 0:
			myExcisableRegions = sorted(myExcisableRegions, key=lambda reg: reg[4], reverse=True)
			newSpeciesLocus = ""
			for mer in myExcisableRegions:
				region = mer[2][1:-1]

				genomes = region.split(",")
				child_nodes = []
				#~ print "genomes", len(genomes), len(myTreeString), len(mer[2])
				#~ print genomes
				for g in genomes:
					#~ print g
					nome = g.split(":")[0]
					dist = float(g.split(":")[1])
					locus = ""
					if nome in self.locusToGenome:
						locus = nome
					else:
						locus = self.codeGenomeID(nome)

					child_nodes.append((locus,dist))
				same_length = 0
				if len(child_nodes) ==1:
					print "child nodes ==1"
					sys.exit()
					
				has_dist = 0
				if len(mer[2]) == len(myTreeString) and len(genomes)==2:
					same_length =1
					tnodes = []
					for cn in child_nodes:
						self.tree.add_node(cn[0])
						tnodes.append(cn[0])
					self.tree.add_edge(tnodes[0],tnodes[1],weight=child_nodes[0][1])
					myTreeString = region
					break
				while len(child_nodes) > 1 and same_length == 0:
					#~ print "WHILE", len(child_nodes)
					if len(child_nodes)>2 and child_nodes[0][1]>0.0:
						has_dist = 1
						#~ print "has_dist", has_dist
						break
					tkids = []
					#~ print len(child_nodes)
					tkids.append(child_nodes.pop(0))
					tkids.append(child_nodes.pop(0))
					#~ print len(child_nodes)
					tnodes = []
					for tk in tkids:
						self.tree.add_node(tk[0])
						tnodes.append(tk[0])
					tnodes.sort()
					newSpeciesLocus = self.codeGenomeID(";".join(tnodes))
					#~ print tkids,";".join(tnodes), newSpeciesLocus
					self.tree.add_node(newSpeciesLocus)
					weight = 0.0
					for tk in tkids:
						weight = tk[1]
						self.tree.add_edge(newSpeciesLocus,tk[0],weight=tk[1])
						#~ print tk[0],self.tree.edge[tk[0]]
					#~ print newSpeciesLocus, self.tree.edge[newSpeciesLocus]
					child_nodes.append((newSpeciesLocus,weight))
				if has_dist ==1:
					min_pair = (child_nodes[0], child_nodes[1])
					min_dist = child_nodes[0][1]+child_nodes[1][1]
					#~ print min_pair, min_dist
					subPaths = {}
					for c in child_nodes:
						#~ print c
						if self.tree_string.find(self.locusToGenome[c[0]])>-1:
							self.tree.add_node(c[0])
						for q in child_nodes:
							if q[0]==c[0]:
								continue
							weight_sum = q[1]+c[1]
							if weight_sum<min_dist:
								min_dist = weight_sum
								min_pair = (c, q)
					tkids = []
					#~ print len(child_nodes), child_nodes
					tkids.append(child_nodes.pop(child_nodes.index(min_pair[0])))
					tkids.append(child_nodes.pop(child_nodes.index(min_pair[1])))
					#~ print len(child_nodes), child_nodes
					tnodes = []
					for tk in tkids:
						tnodes.append(tk[0])
					tnodes.sort()
					newSpeciesLocus = self.codeGenomeID(";".join(tnodes))
					#~ print tkids,";".join(tnodes), newSpeciesLocus
					self.tree.add_node(newSpeciesLocus)
					weight = child_nodes[0][1]
					for tk in tkids:
						self.tree.add_edge(newSpeciesLocus,tk[0],weight=tk[1])
					self.tree.add_edge(newSpeciesLocus,child_nodes[0][0],weight=child_nodes[0][1])
					my_dist = ":"+str(child_nodes[0][1])
					newSpeciesLocus = child_nodes[0][0]+my_dist+","+newSpeciesLocus+my_dist
					#~ newSpeciesLocus = "("+",".join(forReplace)+")"
				myTreeString = myTreeString.replace(mer[2],newSpeciesLocus)
				
				#~ print myTreeString
				child_nodes = []
			myExcisableRegions = self.indexParens(myTreeString)
		return myTreeString
			
	# def parseTree2(self):
	# 	myTreeString = self.tree_string[0:-1]
	# 	myExcisableRegions = self.indexParens(myTreeString)
	# 	while len(myExcisableRegions) > 0:
	# 		myExcisableRegions = sorted(myExcisableRegions, key=lambda reg: reg[3], reverse=True)
	# 		#~ print "excisable", len(myExcisableRegions)
	# 		newSpeciesLocus = ""
	# 		for mer in myExcisableRegions:
	# 			region = mer[2][1:-1]
	# 			genomes = region.split(",")
	# 			child_nodes = []
	# 			#~ print "genomes", len(genomes)
	# 			#~ print genomes
	# 			for g in genomes:
	# 				#~ print g
	# 				nome = g.split(":")[0]
	# 				dist = float(g.split(":")[1])
	# 				locus = ""
	# 				if nome in self.locusToGenome:
	# 					locus = nome
	# 				else:
	# 					locus = self.codeGenomeID(nome)

	# 				child_nodes.append((locus,dist))
	# 			if len(child_nodes) ==1:
	# 				print "child nodes ==1"
	# 				sys.exit()
	# 			if len(mer[2]) == len(myTreeString) and len(genomes)==2:
	# 				same_length =1
	# 				tnodes = []
	# 				for cn in child_nodes:
	# 					self.tree.add_node(cn[0])
	# 					tnodes.append(cn[0])
	# 				self.tree.add_edge(tnodes[0],tnodes[1],weight=child_nodes[0][1])
	# 				myTreeString = region
	# 				break
					
	# 			has_dist = 0
	# 			while len(child_nodes) > 1:
	# 				#~ print "WHILE", len(child_nodes)
	# 				if len(child_nodes)>2 and child_nodes[0][1]>0.0:
	# 					has_dist = 1
	# 					#~ print "has_dist", has_dist
	# 					break
	# 				tkids = []
	# 				#~ print len(child_nodes)
	# 				tkids.append(child_nodes.pop(0))
	# 				tkids.append(child_nodes.pop(0))
	# 				#~ print len(child_nodes)
	# 				tnodes = []
	# 				for tk in tkids:
	# 					if self.tree_string.find(self.locusToGenome[tk[0]])>-1:
	# 						self.genomes.append(self.locusToGenome[tk[0]])
	# 					self.tree.add_node(tk[0])
	# 					tnodes.append(tk[0])
	# 				tnodes.sort()
	# 				newSpeciesLocus = self.codeGenomeID(";".join(tnodes))
	# 				#~ print tkids,";".join(tnodes), newSpeciesLocus
	# 				self.tree.add_node(newSpeciesLocus)
	# 				weight = 0.0
	# 				for tk in tkids:
	# 					weight = tk[1]
	# 					self.tree.add_edge(newSpeciesLocus,tk[0],weight=tk[1])
	# 					#~ print tk[0],self.tree.edge[tk[0]]
	# 				#~ print newSpeciesLocus, self.tree.edge[newSpeciesLocus]
	# 				child_nodes.append((newSpeciesLocus,weight))
	# 			if has_dist ==1:
	# 				#~ print myTreeString
	# 				min_pair = (child_nodes[0], child_nodes[1])
	# 				min_dist = child_nodes[0][1]+child_nodes[1][1]
	# 				#~ print min_pair, min_dist
	# 				subPaths = {}
	# 				for c in child_nodes:
	# 					#~ print c
	# 					if self.tree_string.find(self.locusToGenome[c[0]])>-1:
	# 						self.genomes.append(self.locusToGenome[c[0]])
	# 						self.tree.add_node(c[0])
	# 					paths = nx.shortest_path_length(self.tree,c[0],None,'weight')
	# 					longest_path = 0.0
	# 					#~ print paths
	# 					for p in paths:
	# 						if p == c[0]:
	# 							continue
	# 						if paths[p] > longest_path:
	# 							longest_path = paths[p]
	# 							#~ print longest_path
	# 					longest_path += c[1]
	# 					subPaths[c[0]] = longest_path
	# 				for c in child_nodes:
	# 					for n in child_nodes:
	# 						if c==n:
	# 							continue
	# 						myDist = subPaths[c[0]]+subPaths[n[0]]
	# 						#~ print c[0], subPaths[c[0]], n[0], subPaths[n[0]], myDist
	# 						if myDist < min_dist:
	# 							min_dist = myDist
	# 							min_pair = (c, n)
	# 							#~ print min_pair, min_dist
	# 				tkids = []
	# 				#~ print len(child_nodes), child_nodes
	# 				tkids.append(child_nodes.pop(child_nodes.index(min_pair[0])))
	# 				tkids.append(child_nodes.pop(child_nodes.index(min_pair[1])))
	# 				#~ print len(child_nodes), child_nodes
	# 				tnodes = []
	# 				for tk in tkids:
	# 					tnodes.append(tk[0])
	# 				tnodes.sort()
	# 				newSpeciesLocus = self.codeGenomeID(";".join(tnodes))
	# 				#~ print tkids,";".join(tnodes), newSpeciesLocus
	# 				self.tree.add_node(newSpeciesLocus)
	# 				weight = child_nodes[0][1]
	# 				for tk in tkids:
	# 					self.tree.add_edge(newSpeciesLocus,tk[0],weight=tk[1])
	# 				self.tree.add_edge(newSpeciesLocus,child_nodes[0][0],weight=child_nodes[0][1])
					
	# 				newSpeciesLocus = child_nodes[0][0]
	# 				#~ newSpeciesLocus = "("+",".join(forReplace)+")"

	# 			myTreeString = myTreeString.replace(mer[2],newSpeciesLocus)
	# 			#~ print myTreeString
	# 			child_nodes = []
	# 		myExcisableRegions = self.indexParens(myTreeString)
	# 	print myTreeString
	# 	if len(self.genomes) == 1:
	# 		return "child nodes == 1"
	# 	else:
	# 		print "GENOMES", len(self.genomes)
	# 		return "Good tree"
			
	def findCentroid(self):
		minDist = -1.0
		minRoot = ""
		minRootLen = 0
		for rg in self.genomes:
			print minDist, minRootLen
			r = self.genomeToLocus[rg]
			dist = 0.0
			for ng in self.genomes:
				n = self.genomeToLocus[ng]
				dist+=nx.shortest_path_length(self.tree,r,n,'weight')
			print dist
			if dist < minDist or minDist == -1.0:
				minDist = dist
				minRoot = r
				minRootLen = int(rg.split(";")[1])
			elif dist == minDist:
				root_len = int(rg.split(";")[1])
				if root_len > minRootLen:
					minRoot = r
					minRootLen = root_len
		return self.locusToGenome[minRoot]
			
	def trimTaxa(self, taxaToKeep):
		ancestral = set([])
		for n in self.tree.nodes():
			n_nome = self.locusToGenome[n]
			if n_nome.find(";")>-1:
				ancestral.add(n)
				continue
			if not n_nome in taxaToKeep:
				#~ print "removed", n_nome
				self.tree.remove_node(n)
		for n in self.tree.nodes():
			n_nome = self.locusToGenome[n]
			if n_nome.find(";") > -1 and len(self.tree[n])<2:
				self.tree.remove_node(n)
		badAncestral = 1
		while (badAncestral):
			badAncestral = 0
			for n in self.tree.nodes():
				if n in ancestral and len(self.tree[n])==2:
					badAncestral+=1
					myedges = []
					myweights = []
					for e in self.tree[n]:
						myedges.append(e)
						myweights.append(self.tree[n][e]['weight'])
					newweight = myweights[0]+myweights[1]
					#~ print myedges,myedges[0],myedges[1], newweight
					self.tree.add_edge(myedges[0],myedges[1],weight=newweight)
					self.tree.remove_node(n)
	
	#Roots tree by midpoint based on branch length.
	def rootByMidpoint(self):
		#get all shortest paths between all nodes
		paths = nx.shortest_path_length(self.tree,None,None,'weight')
		#identify the pair of nodes with the longes shortest path. ancestral nodes are also evaluated, but can never be a member of the longest shortest path.
		longest_path = 0.0
		longest_pair = None
		path_dist = {}
		for p in paths:
			path_dist[p] = {}
			for r in paths[p]:
				#~ print r, paths[p][r]
				if r==p:
					path_dist[p][r] = 0.0
					#~ continue
				if r in path_dist[p]:
					if paths[p][r] < path_dist[p][r]:
						path_dist[p][r] = paths[p][r]
				else:
					path_dist[p][r] = paths[p][r]
				if paths[p][r] > longest_path:
					longest_pair = (p, r)
					longest_path = paths[p][r]
				elif longest_pair == None:
					longest_pair = (p,r)
		#trace the longest shortest path to identify the edge containing the midpoint
		path = nx.shortest_path(self.tree, longest_pair[0],longest_pair[1])
		#~ print "longest pair",longest_pair, longest_path
		mid = longest_path/2.0
		cur_node = path.pop()
		cur_length = 0.0
		root_edge = None
		while path:
			nextNode = path.pop()
			cur_length+= self.tree[cur_node][nextNode]['weight']
			if cur_length > mid:
				root_edge = (cur_node,nextNode)
				break
			cur_node = nextNode
		if root_edge == None:
			#this basically means all of the sequences are identical
			print "no root edge, setting to 'longest' pair"
			root_edge = longest_pair
		
		#root tree by adding a root node on midpoint edge
		#~ for n in self.tree.nodes():
			#~ print n, len(self.tree.edge[n]), self.tree.edge[n]
		print "root edge", root_edge
		self.rootTree(root_edge)



	#make a directed graph
	def rootTree(self, root_edge):
		re_weight = (self.tree.edge[root_edge[0]][root_edge[1]]['weight'])/2.0
		self.tree.remove_edge(root_edge[0], root_edge[1])
		self.tree.add_node("root")

		mod_root_edge = []
		for re in root_edge:
			print re, len(self.tree.edges(re)), self.tree.edges(re)
			if not len(self.tree.edges(re))==1:
				print "edge from root to", re
				self.tree.add_edge("root",re,weight=re_weight)
				mod_root_edge.append(re)
			else:
				print "single", re, self.tree.edge[re]
				new_rooter = ""
				for e in self.tree.edge[re]:
					if not e == re:
						new_rooter = e
				myWeight = re_weight+self.tree.edge[new_rooter][re]['weight']
				self.tree.remove_edge(new_rooter, re)
				self.tree.remove_node(re)
				self.tree.add_edge("root",new_rooter,weight=myWeight)
				mod_root_edge.append(new_rooter)
				re = new_rooter
			print re, len(self.tree.edges(re)), self.tree.edge[re]
		self.rooted_tree = self.tree.to_directed()
		paths = nx.shortest_path(self.rooted_tree,source="root")
		for e in self.tree.edges():

			#the node closer to the root becomes the ancestral node of the farther node
			e0 = len(paths[e[0]])
			e1 = len(paths[e[1]])
			if e0<e1:
				self.rooted_tree.remove_edge(e[1],e[0])
			else:
				self.rooted_tree.remove_edge(e[0],e[1])
		ttree = self.rooted_tree.copy()
		up_tree = nx.DiGraph()
		while len(ttree.nodes()) > 1:
			remove_these = set([])
			new_map = {}
			for n in ttree.nodes():
				if n in remove_these:
					pass
				elif len(ttree.out_edges(n)) == 0:
					#~ print "0 out edges",n,self.locusToGenome[n]
					#~ print "ttree.edges(n)",ttree.in_edges(n)
					#~ up_tree.add_node(n)
					parent = self.rooted_tree.in_edges(n)[0][0]
					kids = []
					for oe in ttree.out_edges(parent):
						if len(ttree.out_edges(oe[1])) == 0:
							kids.append(oe[1])
					if len(kids) ==2:
						kids.sort()
						#~ print kids
						up_parent = self.codeGenomeID(";".join(kids))
						if parent == "root":
							up_parent = "root"
						new_map[up_parent] = kids
						
			for nm in new_map:
				
				if not nm in up_tree.nodes():
					up_tree.add_node(nm)
				for ki in new_map[nm]:
					if not ki in up_tree.nodes():
						up_tree.add_node(ki)
					up_tree.add_edge(nm,ki)
					remove_these.add(ki)
			#~ print "removal length", len(remove_these), remove_these
			for rt in remove_these:
				if rt in ttree.nodes():
					ttree.remove_node(rt)
			#~ print "ttree nodes",len(ttree.nodes()), ttree.nodes()
			
		

		
		for n in self.rooted_tree.nodes():
			#~ print n
			if len(self.rooted_tree.out_edges(n)) > 0:
				self.ancestral.append(n)
		
	def indexParens(self, tree_string):
		left_parens = []
		right_parens = []
		parens = []
		left = tree_string.find("(")

		while left < len(tree_string)-1 and left > -1:
			left_parens.append(left)
			parens.append((left,"("))
			left = tree_string.find("(",left+1)

		right = tree_string.find(")",0)
		while right < len(tree_string) and right > -1:
			right_parens.append(right)
			parens.append((right,")"))
			right = tree_string.find(")",right+1)
		parens = sorted(parens, key=lambda tup: tup[0])
		myRegions = []
		i = 1
		while i<len(parens):
			p = parens[i]
			# might want to add a verification that the tree is correctly formated first
			if p[1] == ")" and parens[i-1][1] == "(":
				l = parens[i-1][0]
				r = p[0]+1
				length = r-l
				#~ print tree_string[l:r]
				later = tree_string[r:]
				myweight = re.match(':\d+\.\d+',later)
				region_weight = "None"
				if not myweight ==None:
					region_weight = later[myweight.start()+1:myweight.end()]
				#~ print "looking for", region_weight
				myRegions.append((l,p[0],tree_string[l:r], length,region_weight))
			i+=1
		return myRegions
		
	def writeLocusTagFile(self):
		tag_out = open(self.genomeToLocusFile, 'w')
		for g in self.genomeToLocus:
			line = "\t".join([g, self.genomeToLocus[g]])+"\n"
			tag_out.write(line)
		tag_out.close()
		print "Wrote locus tags to locus_tag_file.txt"
		
	def toNewickLabels(self):
		graph = self.rooted_tree.copy()
		up = [] #unprocessed
		leaf = []
		for n in graph.nodes():
			if len(graph.out_edges(n)) > 0:
				up.append(n)
			elif len(graph.out_edges(n)) ==0:
				leaf.append(n)
				
		curNode = None
		last_string = ""
		up.sort()
		leaf.sort()
		if len(graph.nodes()) == 2:
			last_string = "("+",".join(leaf)+")"
		while len(up) > 0:
			(curNode,e_count) = self.calcMostEdgesToLeaves(up,leaf,graph)
			leaves = []
			for e in graph[curNode]:
				for l in leaf:
					if l == e:
						e_i = leaf.index(l)
						e_text = e
						if 'child_newick' in graph.node[e]:
							if e_count > 2 and len(up) > 1:
								continue
							e_text = graph.node[e]['child_newick']+e
						leaf.pop(e_i)
						leaves.append(e_text)
			#add newick text to curNode
			node_text = "("+",".join(leaves)+")"
			last_string = node_text
			graph.node[curNode]['child_newick'] = node_text
			#change curNode to leaf
			cn_i = up.index(curNode)
			up.pop(cn_i)
			leaf.append(curNode)
		if len(leaf) == 2 and len(up)==0 and len(graph.nodes()) > 2:
			last_string = "("+graph.node[leaf[0]]['child_newick']+","+graph.node[leaf[1]]['child_newick']+")"
		last_string = last_string.rstrip()
		return last_string+";"
		
	def toNewickNoLabels(self):
		graph = self.rooted_tree.copy()
		up = [] #unprocessed
		leaf = []
		for n in graph.nodes():
			if len(graph.out_edges(n)) > 0:
				up.append(n)
			elif len(graph.out_edges(n)) ==0:
				leaf.append(n)
				
		curNode = None
		last_string = ""
		up.sort()
		leaf.sort()
		if len(graph.nodes()) == 2:
			last_string = "("+",".join(leaf)+")"
		while len(up) > 0:
			(curNode,e_count) = self.calcMostEdgesToLeaves(up,leaf,graph)
			leaves = []
			for e in graph[curNode]:
				for l in leaf:
					if l == e:
						e_i = leaf.index(l)
						e_text = e
						if 'child_newick' in graph.node[e]:
							if e_count > 2 and len(up) > 1:
								continue
							e_text = graph.node[e]['child_newick']
						leaf.pop(e_i)
						leaves.append(e_text)
			#add newick text to curNode
			node_text = "("+",".join(leaves)+")"
			last_string = node_text
			graph.node[curNode]['child_newick'] = node_text
			#change curNode to leaf
			cn_i = up.index(curNode)
			up.pop(cn_i)
			leaf.append(curNode)
		if len(leaf) == 2 and len(up)==0 and len(graph.nodes()) > 2:
			last_string = "("+graph.node[leaf[0]]['child_newick']+","+graph.node[leaf[1]]['child_newick']+")"
		last_string = last_string.rstrip()
		return last_string+";"


	def calcMostEdgesToLeaves(self,unprocN,leaf,TG):
		mostLeaves = -1
		retNode = None
		for n in unprocN:
			e_count = 0
			for e in TG[n]:
				for l in leaf:
					if e == l:
						e_count += 1
			if e_count > mostLeaves:
				mostLeaves = e_count
				retNode = n
		return (retNode,mostLeaves)
		
	def makeClassSplits(self):
		ancestor_to_children = {}
		for a in self.ancestral:
			ancestor_to_children[a] = set([])
			leaves = []
			kids = self.getLeaves(a, leaves)
			for k in kids:
				ancestor_to_children[a].add(self.locusToGenome[k])
		return ancestor_to_children
		
	def getLeaves(self, ancestral_node, leaves):
		ret_list = []
		for c in self.rooted_tree.successors(ancestral_node):
			if len(self.rooted_tree.out_edges(c)) == 0:
				leaves.append(c)
			else:
				self.getLeaves(c, leaves)
		return leaves
			
def usage():
	print """
	python tree_poke.py [options]
	-t, --tree [file]
		where [file] is a species tree in newick format without bootstrap values
	-s, --sequence [file]
		where [file] is a fasta file where the headers correspond with the leaf nodes in the tree
	-n, --node
		if flagged, ancestral nodes will be labelled.  Mappings will be stored in the locus file
	-l, --locus [file]
		where [file] is the location of the locus mappings, stored in locus_tag_file.txt by default
	-r, --root
		if flagged, tree will be rooted by the midpoint based on branch length distances
		CURRENTLY ALL TREES WILL BE ROOTED. OOPS?
	-k, --keep [file]
		where [file] is a one taxa per line file of taxa to retain in the tree.  The topology of the existing tree will be retained.
		No labels will be retained or reported.
	-h, --help
		prints this and exits
	"""
def main(argv):
	tree_file = ""
	seq_file = ""
	out_file = ""
	locus_tag = "locus_tag_file.txt"
	root = 1 #always on... need to develop
	labels = 0
	taxa_to_keep = ""
	
	try:
		opts, args = getopt.getopt(argv, "t:l:k:s:o:rnh",["tree=","locus=","keep=","sequence=","output=","root","node","help"])

	except getopt.GetoptError:
		usage()
		sys.exit(2)
	for opt, arg in opts:
		if opt in ("-t","--tree"):
			tree_file = arg
		elif opt in ("-s","--sequence"):
			seq_file = arg
		elif opt in ("-o","--output"):
			out_file = arg
		elif opt in ("-l","--locus"):
			locus_tag = arg
		elif opt in ("-r","--root"):
			root = 1
		elif opt in ("-n","--node"):
			labels = 1
		elif opt in ("-k","--keep"):
			taxa_to_keep = arg
		elif opt in ("-h","--help"):
			sys.exit(usage())
	
	myTree = Tree(tree_file, locus_tag)
	clusterID = (out_file.split("/")[-1]).split(".")[0]
	mySeqs = open(seq_file,'r').readlines()
	seqs = {}
	curseq = ""
	centroidSeq = ""

	for s in mySeqs:
		s = s.rstrip()
		if s.find(">")>-1:
			curseq = s[1:]
			seqs[curseq] = ""
		else:
			seqs[curseq] = seqs[curseq]+s
	if len(seqs)==1:
		for s in seqs:
			centroidSeq=seqs[s]
	else:

		myTree.readTree()
		success = myTree.parseTree()
		myCentroid = myTree.findCentroid()
		centroidSeq = seqs[myCentroid]
	conPep = open(out_file,'w')
	conPep.write(">"+clusterID+";"+str(len(centroidSeq))+"\n")
	conPep.write(centroidSeq+"\n")
	conPep.close()

	
if __name__ == "__main__":
	if len(sys.argv) == 1:
		sys.exit(usage())
	else:
		main(sys.argv[1:])