#!/usr/bin/env python

import sys
import getopt
import re
import logging
import hashlib
import base64
import networkx as nx


class Tree:
	logger = logging.getLogger("Tree")

	def __init__(self, tree_file, genomeToLocusFile):
		self.tree_file = tree_file
		self.genomeToLocusFile = genomeToLocusFile
		self.genomeToLocus = {}
		self.locusToGenome = {}
		self.nodeChildrenCount = {}
		self.tree_string = ""
		self.tree = nx.Graph()
		self.root = None
		self.ancestral = []
		self.rooted_tree = nx.DiGraph()
		self.genomes = []

	def readTree(self):
		with open(self.tree_file) as t:
			tree = t.read()
			trees = re.split('\)\d+[\.\d+]*:', tree)
			jstring = "):"
			tree = jstring.join(trees)
			tree = tree.rstrip()
			print tree
			self.tree_string = tree

	def readGenomeToLocusFile(self):
		tags = open(self.genomeToLocusFile, 'r').readlines()
		for t in tags:
			t = t.rstrip()
			l = t.split()
			if not len(l) > 1:
				continue
			self.genomeToLocus[l[0]] = l[1]
			self.locusToGenome[l[1]] = l[0]
			self.nodeChildrenCount[l[1]] = int(l[2])

	def codeGenomeID(self, genome):
		tag = ''
		if genome in self.genomeToLocus:
			tag = self.genomeToLocus[genome]
			self.locusToGenome[tag] = genome
		else:
			if genome.count(";") > 0:
				children = genome.split(";")
				children.sort()
				# max node degree between children nodes
				degree = max(int(children[0].split("_")[1]), int(children[1].split("_")[1])) + 1
				n = self.nodeChildrenCount[children[0]] + self.nodeChildrenCount[children[1]]
			else:  # should never happen since otherwise its a leaf
				exit("Error: Genome " + genome + " found in the tree but not in the repo_spec.")
				degree = 0
			tag = "N_%07d_%s" % (degree, base64.urlsafe_b64encode(hashlib.md5(genome).digest())[:-2])
			Tree.logger.debug("Created new tag %s for %s" % (tag, genome))
			self.genomeToLocus[genome] = tag
			self.locusToGenome[tag] = genome
			self.nodeChildrenCount[tag] = n
		return tag

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
				for g in genomes:
					if g[0] == "'":
						nome = g[1:g.rfind("'")]
					elif g[0] == '"':
						nome = g[1:g.rfind('"')]
					else:
						nome = ":".join(g.split(":")[:-1])
					dist = float(g.split(":")[-1])
					locus = ""
					if nome in self.locusToGenome:
						locus = nome
					else:  # will error as no ";" because leaf
						locus = self.codeGenomeID(nome)

					child_nodes.append((locus, dist))
				same_length = 0
				if len(child_nodes) == 1:
					print "child nodes ==1"
					sys.exit()

				has_dist = 0
				if len(mer[2]) == len(myTreeString) and len(genomes) == 2:
					same_length = 1
					tnodes = []
					for cn in child_nodes:
						self.tree.add_node(cn[0])
						tnodes.append(cn[0])
					self.tree.add_edge(tnodes[0], tnodes[1], weight=child_nodes[0][1])
					myTreeString = region
					break
				while len(child_nodes) > 1 and same_length == 0:
					if len(child_nodes) > 2 and child_nodes[0][1] > 0.0:
						has_dist = 1
						break
					if len(child_nodes) == 2 and len(mer[2]) == len(myTreeString):
						newSpeciesLocus = child_nodes[0][0] + ":" + str(child_nodes[0][1]) + "," + child_nodes[1][0] + ":" + str(child_nodes[1][1])
						self.tree.add_edge(child_nodes[0][0], child_nodes[1][0], weight=child_nodes[0][1])
						break
					tkids = []
					tkids.append(child_nodes.pop(0))
					tkids.append(child_nodes.pop(0))
					tnodes = []
					for tk in tkids:
						self.tree.add_node(tk[0])
						tnodes.append(tk[0])
					tnodes.sort()
					newSpeciesLocus = self.codeGenomeID(";".join(tnodes))
					self.tree.add_node(newSpeciesLocus)
					weight = 0.0
					for tk in tkids:
						weight = tk[1]
						self.tree.add_edge(newSpeciesLocus, tk[0], weight=tk[1])
					child_nodes.append((newSpeciesLocus, weight))
				if has_dist == 1:
					min_pair = (child_nodes[0], child_nodes[1])
					min_dist = child_nodes[0][1] + child_nodes[1][1]
					for c in child_nodes:
						if self.tree_string.find(self.locusToGenome[c[0]]) > -1:
							self.tree.add_node(c[0])
						for q in child_nodes:
							if q[0] == c[0]:
								continue
							weight_sum = q[1] + c[1]
							if weight_sum < min_dist:
								min_dist = weight_sum
								min_pair = (c, q)
					tkids = []
					tkids.append(child_nodes.pop(child_nodes.index(min_pair[0])))
					tkids.append(child_nodes.pop(child_nodes.index(min_pair[1])))
					tnodes = []
					for tk in tkids:
						tnodes.append(tk[0])
					tnodes.sort()
					newSpeciesLocus = self.codeGenomeID(";".join(tnodes))
					self.tree.add_node(newSpeciesLocus)
					weight = child_nodes[0][1]
					for tk in tkids:
						self.tree.add_edge(newSpeciesLocus, tk[0], weight=tk[1])
					self.tree.add_edge(newSpeciesLocus, child_nodes[0][0], weight=child_nodes[0][1])
					my_dist = ":" + str(child_nodes[0][1])
					newSpeciesLocus = child_nodes[0][0] + my_dist + "," + newSpeciesLocus + my_dist
				myTreeString = myTreeString.replace(mer[2], newSpeciesLocus)
				child_nodes = []
			myExcisableRegions = self.indexParens(myTreeString)
		return myTreeString

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
				dist += nx.shortest_path_length(self.tree, r, n, 'weight')
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

	# make a directed graph
	def rootTree(self, root_edge):
		re_weight = (self.tree.edge[root_edge[0]][root_edge[1]]['weight']) / 2.0
		self.tree.remove_edge(root_edge[0], root_edge[1])
		self.root = self.codeGenomeID(";".join(root_edge))
		self.tree.add_node(self.root)

		mod_root_edge = []
		for redge in root_edge:
			print redge, len(self.tree.edges(redge)), self.tree.edges(redge)
			if not len(self.tree.edges(redge)) == 1:
				print "edge from root to", redge
				self.tree.add_edge(self.root, redge, weight=re_weight)
				mod_root_edge.append(redge)
			else:
				print "single", redge, self.tree.edge[redge]
				new_rooter = ""
				for e in self.tree.edge[redge]:
					if not e == redge:
						new_rooter = e
				myWeight = re_weight + self.tree.edge[new_rooter][redge]['weight']
				self.tree.remove_edge(new_rooter, redge)
				self.tree.remove_node(redge)
				self.tree.add_edge(self.root, new_rooter, weight=myWeight)
				mod_root_edge.append(new_rooter)
				redge = new_rooter
			print redge, len(self.tree.edges(redge)), self.tree.edge[redge]
		self.rooted_tree = self.tree.to_directed()
		paths = nx.shortest_path(self.rooted_tree, source=self.root)
		for e in self.tree.edges():

			# the node closer to the root becomes the ancestral node of the farther node
			e0 = len(paths[e[0]])
			e1 = len(paths[e[1]])
			if e0 < e1:
				self.rooted_tree.remove_edge(e[1], e[0])
			else:
				self.rooted_tree.remove_edge(e[0], e[1])

		for n in self.rooted_tree.nodes():
			# print n
			if len(self.rooted_tree.out_edges(n)) > 0:
				self.ancestral.append(n)

	def indexParens(self, tree_string):
		left_parens = []
		right_parens = []
		parens = []
		left = tree_string.find("(")

		while left < len(tree_string) - 1 and left > -1:
			left_parens.append(left)
			parens.append((left, "("))
			left = tree_string.find("(", left + 1)

		right = tree_string.find(")", 0)
		while right < len(tree_string) and right > -1:
			right_parens.append(right)
			parens.append((right, ")"))
			right = tree_string.find(")", right + 1)
		parens = sorted(parens, key=lambda tup: tup[0])
		myRegions = []
		i = 1
		while i < len(parens):
			p = parens[i]
			# might want to add a verification that the tree is correctly formated first
			if p[1] == ")" and parens[i - 1][1] == "(":
				l = parens[i - 1][0]
				r = p[0] + 1
				length = r - l
				# print tree_string[l:r]
				later = tree_string[r:]
				myweight = re.match(':\d+\.\d+', later)
				region_weight = "None"
				if myweight is not None:
					region_weight = later[myweight.start() + 1:myweight.end()]
				# print "looking for", region_weight
				myRegions.append((l, p[0], tree_string[l:r], length, region_weight))
			i += 1
		return myRegions


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
	# root = 1  # always on... need to develop
	# labels = 0
	# taxa_to_keep = ""

	try:
		opts, args = getopt.getopt(argv, "t:l:k:s:o:rnh", ["tree=", "locus=", "keep=", "sequence=", "output=", "root", "node", "help"])

	except getopt.GetoptError:
		usage()
		sys.exit(2)
	for opt, arg in opts:
		if opt in ("-t", "--tree"):
			tree_file = arg
		elif opt in ("-s", "--sequence"):
			seq_file = arg
		elif opt in ("-o", "--output"):
			out_file = arg
		elif opt in ("-l", "--locus"):
			locus_tag = arg
		elif opt in ("-h", "--help"):
			sys.exit(usage())

	myTree = Tree(tree_file, locus_tag)
	clusterID = (out_file.split("/")[-1]).split(".")[0]
	mySeqs = open(seq_file, 'r').readlines()
	seqs = {}
	curseq = ""
	centroidSeq = ""

	for s in mySeqs:
		s = s.rstrip()
		if s.find(">") > -1:
			curseq = s[1:]
			seqs[curseq] = ""
		else:
			seqs[curseq] = seqs[curseq] + s
	if len(seqs) == 1:
		for s in seqs:
			centroidSeq = seqs[s]
	else:

		myTree.readTree()
		# success = myTree.parseTree()
		myTree.parseTree()
		myCentroid = myTree.findCentroid()
		centroidSeq = seqs[myCentroid]
	conPep = open(out_file, 'w')
	conPep.write(">" + clusterID + ";" + str(len(centroidSeq)) + "\n")
	conPep.write(centroidSeq + "\n")
	conPep.close()


if __name__ == "__main__":
	if len(sys.argv) == 1:
		sys.exit(usage())
	else:
		main(sys.argv[1:])
