#!/usr/bin/env python

# Gets the filesystem set prepped to run the actual algorithm and writes the commands to do so
import os
import logging
import collections


class Tree:
	logger = logging.getLogger("Tree")

	def __init__(self, tree_obj, blast_eval, num_cores, min_best_hit, syn_dist, minSynFrac, mutrate, synteny):
		self.tree_obj = tree_obj
		self.genomeToLocusFile = tree_obj.genomeToLocusFile
		self.genomeToLocus = tree_obj.genomeToLocus
		self.locusToGenome = tree_obj.locusToGenome
		self.nodeChildrenCount = tree_obj.nodeChildrenCount
		self.tree_string = ""
		self.tree = tree_obj.tree
		self.rooted_tree = tree_obj.rooted_tree
		self.root = None
		self.blast_eval = blast_eval
		self.num_cores = num_cores
		self.min_best_hit = min_best_hit
		self.syn_dist = int(syn_dist)
		self.min_syn_frac = minSynFrac
		self.mutrate = mutrate
		self.synteny = synteny
		self.syn2_path = "#SYNERCLUST_PATH"
		Tree.logger.debug("Tree initialized")

	def codeGenomeID(self, genome):
		return self.genomeToLocus[genome]

	def reregisterGenomeID(self, identifier, newChildren):  # actually usefull?
		oldGenome = self.locusToGenome[identifier]
		newGenome = ";".join(newChildren)
		del self.genomeToLocus[oldGenome]
		self.genomeToLocus[newGenome] = identifier
		self.locusToGenome[identifier] = newGenome
		self.nodeChildrenCount[identifier] = self.nodeChildrenCount[oldGenome]
		del self.nodeChildrenCount[oldGenome]

	def writeLocusTagFile(self):
		tag_out = open(self.genomeToLocusFile, 'w')
		for g in self.genomeToLocus:
			line = "\t".join([g, self.genomeToLocus[g], str(self.nodeChildrenCount[self.genomeToLocus[g]])]) + "\n"
			tag_out.write(line)
		tag_out.close()
		Tree.logger.info("Wrote locus tags to locus_tag_file.txt")

	def writeCodedNewick(self, coded_nwk_out):
		to_replace = [self.root]
		nwk_str = self.root + ":" + str(self.nodeChildrenCount[self.root])
		while to_replace:
			n = to_replace.pop()
			children_tag = self.rooted_tree.out_edges(n)  # should always be 2 if internal node and 0 if leaf
			if not children_tag:
				continue
			children_tag = [children_tag[0][1], children_tag[1][1]]
			insert_str = ""
			pos = nwk_str.find(n)
			if children_tag:
				insert_str = "(" + children_tag[0] + ":" + str(self.nodeChildrenCount[children_tag[0]]) + "," + children_tag[1] + ":" + str(self.nodeChildrenCount[children_tag[1]]) + ")"
			nwk_str = nwk_str[:pos] + insert_str + nwk_str[pos:]
			to_replace.extend([children_tag[0], children_tag[1]])
		nwk_str += ";\n"
		with open(coded_nwk_out, 'w') as f:
			f.write(nwk_str)

	def calculateMinTotalSubTiers(self, subnodes, node):
		minTotal = 1000000000
		minPair = None
		# print "node",node,"subtiers",subnodes
		for r in subnodes:
			for s in subnodes:
				if r == s:
					continue
				# print subnodes[r]+subnodes[s], r, s
				if (subnodes[r] + subnodes[s]) < minTotal:
					# print "min!", r, s
					Tree.logger.info("min! %s %s" % (str(r), str(s)))
					minTotal = subnodes[r] + subnodes[s]
					minPair = (r, s)
		minList = [minPair[0], minPair[1]]
		minList.sort()
		minPair = (minList[0], minList[1])
		return minPair

	def calculateNodeDependencies(self, working_dir):
		node_dir = working_dir + "nodes"
		if "nodes" not in os.listdir(working_dir):
			os.system("mkdir " + node_dir)
		gene_nodes = []
		gene_nodes_num_tiers = {}
		noGene_nodes = {}  # node maps to children
		ready_nodes = []
		n_list = os.listdir(node_dir + "/")
		pick_count = 0
		for n in self.rooted_tree.nodes():
			# print "rooted_tree node", n, self.rooted_tree[n]
			noGene_nodes[n] = []
			for e in self.rooted_tree.edges(n):
				# print e
				noGene_nodes[n].append(e[1])
			if n not in n_list:
				os.system("mkdir " + node_dir + "/" + n)
			if len(noGene_nodes[n]) == 0:
				pick_count += 1

		cmd_count = 0  # also tiers below
		nodeTier = {}
		childToParent = {}
		while len(noGene_nodes.keys()) > 0:
			new_nodes = []
			for n in noGene_nodes:
				ready = 1
				for d in noGene_nodes[n]:
					if d not in gene_nodes:
						ready = 0
				if ready == 1:
					# this section handles 3-way merges which, while legitimate, cause a ton of headaches for the rest of this code
					if len(noGene_nodes[n]) > 2:
						# print "THREE-WAY MERGE", noGene_nodes[n]
						Tree.logger.info("THREE-WAY MERGE %s" % (noGene_nodes[n]))
						# calc number of sub tiers below each child node
						subs = {}
						for q in noGene_nodes[n]:
							subs[q] = gene_nodes_num_tiers[q]
						# find pair of children with minimum # of sub tiers
						ready_pair = self.calculateMinTotalSubTiers(subs, n)
						excluded = ""
						for q in noGene_nodes[n]:
							if q not in ready_pair:
								excluded = q
						# make a new node using ready_pair
						intermediate_node = self.codeGenomeID(";".join([ready_pair[0], ready_pair[1]]))
						Tree.logger.info("new intermediate node %s %s" % (intermediate_node, ready_pair))
						# print "new intermediate node", intermediate_node, ready_pair
						os.system("mkdir " + node_dir + "/" + intermediate_node)
						new_nodes.append((intermediate_node, ready_pair[0], ready_pair[1]))

						# reset children of n
						noGene_nodes[n] = []
						noGene_nodes[n].append(intermediate_node)
						noGene_nodes[n].append(excluded)
						# modify tree
						self.rooted_tree.add_node(intermediate_node)
						self.rooted_tree.remove_edge(ready_pair[0], n)
						self.rooted_tree.remove_edge(ready_pair[1], n)
						self.rooted_tree.add_edge(ready_pair[0], intermediate_node)
						self.rooted_tree.add_edge(ready_pair[1], intermediate_node)
						self.rooted_tree.add_edge(intermediate_node, n)

						# reregister locus tag of n
						noGene_nodes[n].sort()
						my_id = n
						if ";".join(noGene_nodes[n]) in self.genomeToLocus:
							my_id = self.genomeToLocus[";".join(noGene_nodes[n])]

						self.reregisterGenomeID(my_id, noGene_nodes[n])

						ready_nodes.append(intermediate_node)
					else:
						ready_nodes.append(n)
			if len(new_nodes) > 0:
				for n in new_nodes:
					noGene_nodes[n[0]] = []
					noGene_nodes[n[0]].append(n[1])
					noGene_nodes[n[0]].append(n[2])
				new_nodes = []

			# print "ready_nodes", ready_nodes
			nodeTier[cmd_count] = []
			for r in ready_nodes:
				nodeTier[cmd_count].append(r)
				# print "noGene_nodes[r]",r, noGene_nodes[r]
				for d in noGene_nodes[r]:
					childToParent[d] = r
				del noGene_nodes[r]
				gene_nodes.append(r)
				gene_nodes_num_tiers[r] = cmd_count
			ready_nodes = []
			cmd_count += 1
		# print "c2p\n",childToParent
		self.makeNodeFlowWorkflowControl(nodeTier, childToParent, working_dir)

	def makeNodeFlowWorkflowControl(self, nodeTier, childToParent, working_dir):
		self.root = nodeTier[max(nodeTier)][0]
		stack = []
		queue = collections.deque([[self.root, 0, []]])  # current_node, current_id, [child1_id, child2_id]
		count = 1
		while len(queue) > 0:
			current = queue.pop()
			stack.append(current)
			for e in self.rooted_tree.edges(current[0]):  # get children
				queue.append([e[1], count, []])  # add child to the queue (~breadth first search)
				current[2].append(count)  # add child to parent dependency
				count += 1
		with open(working_dir + "uger_jobs.sh", "w") as uge_out:
			with open(working_dir + "jobs.sh", "w") as out:
				uge_out.write("#! /bin/bash\n\nTIME=$(date +%s)\n")
				while len(stack) > 0:
					current = stack.pop()
					if current[0][0] != "L":  # not a leaf
						uge_out.write("qsub -N j${TIME}" + str(current[1]))
						if len(current[2]) != 0:
							uge_out.write(" -hold_jid j${TIME}" + str(current[2][0]) + ",j${TIME}" + str(current[2][1]))
						uge_out.write(" " + working_dir + "nodes/" + str(current[0]) + "/" + str(current[0]) + "_uge.sh\n")
						out.write(working_dir + "nodes/" + str(current[0]) + "/" + str(current[0]) + ".sh\n")
		os.chmod(working_dir + "uger_jobs.sh", 0775)
		os.chmod(working_dir + "jobs.sh", 0775)
		all_proc_nodes = []
		for n in nodeTier[1]:  # tier 1 is the first tier above the leaf nodes
			curNode = n
			while curNode not in all_proc_nodes:
				kids = []
				for e in self.rooted_tree.edges(curNode):
					kids.append(e[1])
				self.makeSingleNodeFlow(working_dir, curNode, kids)
				all_proc_nodes.append(curNode)
				if curNode not in childToParent:
					break  # curNode is root
				curNode = childToParent[curNode]

	def makeSingleNodeFlow(self, working_dir, curNode, kids):
		my_dir = working_dir + "nodes/" + curNode + "/"
		child1 = kids[0]
		child2 = kids[1]

		syn2_path = self.syn2_path

		sh_uge_file = syn2_path + "NewNodeShTemplate.sh"
		sh_file = syn2_path + "NewNodeTemplate.sh"

		my_sh_uge_file = my_dir + curNode + "_uge.sh"
		my_sh_file = my_dir + curNode + ".sh"

		with open(sh_uge_file, 'r') as f:
			s_file = f.read()
		s_file = s_file.replace('#SYNERCLUST_PATH', syn2_path)
		s_file = s_file.replace('#WORKING_DIR', working_dir)
		s_file = s_file.replace('#CHILD1', child1)
		s_file = s_file.replace('#CHILD2', child2)
		s_file = s_file.replace('#NODE', curNode)
		s_file = s_file.replace('#BLAST_EVAL', str(self.blast_eval))
		s_file = s_file.replace('#NUM_CORES', str(self.num_cores))
		s_file = s_file.replace('#MUTRATE', str(self.mutrate))
		s_file = s_file.replace('#MIN_BEST_HIT', str(self.min_best_hit))
		s_file = s_file.replace('#MIN_SYNTENIC_FRACTION', str(self.min_syn_frac))
		s_file = s_file.replace('#WORKING_DIR', working_dir)
		if self.synteny:
			s_file = s_file.replace('#NOSYNTENY', "")
		else:
			s_file = s_file.replace('#NOSYNTENY', "--no-synteny")
		with open(my_sh_uge_file, 'w') as my_sh:
			my_sh.write(s_file)
		os.chmod(my_sh_uge_file, 0775)

		with open(sh_file, 'r') as f:
			s_file = f.read()
		s_file = s_file.replace('#SYNERCLUST_PATH', syn2_path)
		s_file = s_file.replace('#WORKING_DIR', working_dir)
		s_file = s_file.replace('#CHILD1', child1)
		s_file = s_file.replace('#CHILD2', child2)
		s_file = s_file.replace('#NODE', curNode)
		s_file = s_file.replace('#BLAST_EVAL', str(self.blast_eval))
		s_file = s_file.replace('#NUM_CORES', str(self.num_cores))
		s_file = s_file.replace('#MUTRATE', str(self.mutrate))
		s_file = s_file.replace('#MIN_BEST_HIT', str(self.min_best_hit))
		s_file = s_file.replace('#MIN_SYNTENIC_FRACTION', str(self.min_syn_frac))
		s_file = s_file.replace('#WORKING_DIR', working_dir)
		if self.synteny:
			s_file = s_file.replace('#NOSYNTENY', "")
		else:
			s_file = s_file.replace('#NOSYNTENY', "--no-synteny")
		with open(my_sh_file, 'w') as my_sh:
			my_sh.write(s_file)
		os.chmod(my_sh_file, 0775)

	def getGenomeToLocus(self):
		return self.genomeToLocus

	def getLocusToGenome(self):
		return self.locusToGenome
