#!/usr/bin/env python

# Gets the filesystem set prepped to run the actual algorithm and writes the commands to do so
import os
import pickle
import logging
import collections


class Tree:
	logger = logging.getLogger("Tree")

	def __init__(self, tree_obj, blast_eval, num_cores, alpha, beta, gamma, gain, loss, min_best_hit, syn_dist, minSynFrac, hamming):
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
		self.alpha = alpha
		self.beta = beta
		self.gamma = gamma
		self.gain = gain
		self.loss = loss
		self.min_best_hit = min_best_hit
		self.syn_dist = int(syn_dist)
		self.min_syn_frac = minSynFrac
		self.hamming = hamming
		self.syn2_path = "#SYNERGY2_PATH"
		Tree.logger.debug("Tree initialized")

	def codeGenomeID(self, genome):
		return self.genomeToLocus[genome]
		# tag = ''
		# if genome in self.genomeToLocus:
		# 	tag = self.genomeToLocus[genome]
		# else:
		# 	# TODO CHECK IF THIS IS EVER CALLED
		# 	Tree.logger.error("If this is called, check when and why to modify the tag generated.")
		# 	import pdb
		# 	pdb.set_trace()
		# 	tag = ''.join(random.choice(string.ascii_uppercase) for x in range(3))
		# 	while tag in self.locusToGenome:
		# 		tag = ''.join(random.choice(string.ascii_uppercase) for x in range(3))
		# self.genomeToLocus[genome] = tag
		# self.locusToGenome[tag] = genome
		# return tag

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

	def makePicklesForSingleGenome(self, working_dir, genome, node):
		gdat = open(working_dir + "genomes/" + genome + "/annotation.txt", 'r').readlines()
		# this file might not be necessary in the long-term, in fact these genome directories may not be necessary at all
		# everything should be pulled out and stored in pickles after the repo files are parsed.  Whatever.
		ndat = open(working_dir + "nodes/" + node + "/" + node + ".pep", 'w')
		# print node
		x = gdat[1].split()[1]
		y = "_".join(x.split("_")[:-1])
		if not y == node:
			Tree.logger.warning("%s %s", (node, y))
			Tree.logger.warning("%s" % (gdat[1]))

		# these hashes will be pickles
		genes = {}
		gene_map = {}
		locusToTID = {}
		# this dict is for synteny information, keys are scaffold IDs
		neighbors = {}
		for g in gdat[1:]:
			g = g.rstrip()
			l = g.split()
			if len(l) < 8:
				Tree.logger.warning(g)
				continue
			if not l[2] in neighbors:
				neighbors[l[2]] = []
			gene_tup = (l[1], int(l[3]), int(l[4]), l[5], int(l[6]))  # scaffold -> locus,start,stop,strand,length
			locusToTID[l[1]] = l[0]
			neighbors[l[2]].append(gene_tup)
			genes[l[1]] = l[7]
			gene_map[l[1]] = [l[1]]
			line = ">" + l[1] + ";" + l[6] + "\n" + l[7] + "\n"
			ndat.write(line)
		ndat.close()

		# make synteny pickle
		self.makeSyntenyPickle(working_dir, genome, node, neighbors)
		# dump pickles!
		pdat = open(working_dir + "nodes/" + node + "/" + node + ".pkl", 'wb')
		pickle.dump(genes, pdat)
		pdat.close()
		ldat = open(working_dir + "nodes/" + node + "/locus_mappings.pkl", 'wb')
		pickle.dump(gene_map, ldat)
		ldat.close()
		tdat = open(working_dir + "genomes/" + genome + "/locusToTranscript.pkl", 'wb')
		pickle.dump(locusToTID, tdat)
		tdat.close()
		pcomp = open(working_dir + "nodes/" + node + "/PICKLES_COMPLETE", 'w')
		pcomp.write("Way to go!")
		pcomp.close()
		# print "Pickles complete for ", genome, node

	def makeSyntenyPickle(self, working_dir, genome, node, neighbors):
		# make a pickle!
		gsyn = {}
		nsyn = {}
		MAX_DIST = self.syn_dist
		for n in neighbors:
			# n is a scaffold ID
			genes = neighbors[n]
			genes = sorted(genes, key=lambda tup: tup[1])
			for g in genes:
				# g_i = genes.index(g)
				locus = g[0]
				lend = g[1]
				rend = g[2]
				strand = g[3]
				# length = g[4]
				gsyn[locus] = []
				nsyn[locus] = {}
				nsyn[locus]['neighbors'] = []
				nsyn[locus]['count'] = 1
				gmid = (rend - lend + 1) / 2 + lend
				minmid = lend - MAX_DIST
				maxmid = rend + MAX_DIST
				for h in genes:
					mid = (h[2] - h[1] + 1) / 2 + h[1]
					if h[0] == locus:
						continue
					if mid >= minmid and mid <= maxmid:
						# here, direction AND strand yields stream direction, it's bitwise 'and' or something
						stream = 0  # -1 means upstream, 1 means downstream
						dist = -1
						direction = None
						if h[1] > lend:
							dist = mid - gmid + 1
							direction = "+"
						else:
							dist = gmid - mid + 1
							direction = "-"
						if direction == strand:
							stream = -1
						else:
							stream = 1
						genome_tup = (h[0], h[1], h[2], dist, stream)
						gsyn[locus].append(genome_tup)
						nsyn[locus]['neighbors'].append(h[0])
					if mid > maxmid:
						break
# TODO inspect here
# List of all genes. For each gene, list of all other genes in the form (gene, left end, right end, distance, stream(?))
		gdat = open(working_dir + "genomes/" + genome + "/synteny_data.pkl", 'wb')
		pickle.dump(gsyn, gdat)
		gdat.close()
		ndat = open(working_dir + "nodes/" + node + "/synteny_data.pkl", 'wb')
		pickle.dump(nsyn, ndat)
		ndat.close()

	def makeCommandForNodeFlow(self, node, children, node_dir, cmd):
		c_text = " ".join(children)
		myCmd = cmd + node_dir + "/ " + node + " " + c_text + "\n"
		return myCmd

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
				if "PICKLES_COMPLETE" not in os.listdir(node_dir + "/" + n):
					self.makePicklesForSingleGenome(working_dir, self.locusToGenome[n], n)
					Tree.logger.info("%s %s %s" % (pick_count, n, self.locusToGenome[n]))
					# print pick_count, n, self.locusToGenome[n]

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
		all_proc_nodes = []
		serial_sets = {}
		sets = 1
		for n in nodeTier[1]:  # tier 1 is the first tier above the leaf nodes
			Tree.logger.info("set number %s %s" % (sets, n))
			# print "set number ", sets, n
			local_proc_nodes = []
			curNode = n
			set_cmd_count = 1
			next_cmd_id = "1." + str(sets) + "." + str(set_cmd_count)
			while curNode not in all_proc_nodes:
				kids = []
				for e in self.rooted_tree.edges(curNode):
					# print curNode, e
					kids.append(e[1])
				self.makeSingleNodeFlow(working_dir, curNode, next_cmd_id, kids)
				all_proc_nodes.append(curNode)
				local_proc_nodes.append((next_cmd_id, curNode))
				set_cmd_count += 1
				next_cmd_id = "1." + str(sets) + "." + str(set_cmd_count)
				if curNode not in childToParent:
					break  # curNode is root
				curNode = childToParent[curNode]
			serial_sets[sets] = local_proc_nodes
			sets += 1
			set_cmd_count = 1
		Tree.logger.info(serial_sets)
		# print serial_sets
		# self.makeNodeFlowLauncher(working_dir, serial_sets)

	def makeSingleNodeFlow(self, working_dir, curNode, cmd_id, kids):
		my_dir = working_dir + "nodes/" + curNode + "/"
		child1 = kids[0]
		child2 = kids[1]

		syn2_path = self.syn2_path

		sh_uge_file = syn2_path + "NewNodeShTemplate.sh"
		sh_file = syn2_path + "NewNodeTemplate.sh"

		my_sh_uge_file = my_dir + curNode + "_uge.sh"
		my_sh_file = my_dir + curNode + ".sh"

		s_file = open(sh_uge_file, 'r').read()
		s_file = s_file.replace('#SYNERGY2_PATH', syn2_path)
		s_file = s_file.replace('#WORKING_DIR', working_dir)
		s_file = s_file.replace('#CHILD1', child1)
		s_file = s_file.replace('#CHILD2', child2)
		s_file = s_file.replace('#NODE', curNode)
		s_file = s_file.replace('#ID', cmd_id)
		s_file = s_file.replace('#BLAST_EVAL', str(self.blast_eval))
		s_file = s_file.replace('#NUM_CORES', str(self.num_cores))
		s_file = s_file.replace('#HAMMING', str(self.hamming))
		s_file = s_file.replace('#ALPHA', str(self.alpha))
		s_file = s_file.replace('#BETA', str(self.beta))
		s_file = s_file.replace('#GAMMA', str(self.gamma))
		s_file = s_file.replace('#GAIN', str(self.gain))
		s_file = s_file.replace('#LOSS', str(self.loss))
		s_file = s_file.replace('#MIN_BEST_HIT', str(self.min_best_hit))
		s_file = s_file.replace('#MIN_SYNTENIC_FRACTION', str(self.min_syn_frac))
# 		s_file = s_file.replace('#HOMOLOGY_SCALE', str(self.homScale))
# 		s_file = s_file.replace('#SYNTENY_SCALE', str(self.synScale))
		s_file = s_file.replace('#WORKING_DIR', working_dir)

		my_sh = open(my_sh_uge_file, 'w')
		my_sh.write(s_file)
		my_sh.close()

		s_file = open(sh_file, 'r').read()
		s_file = s_file.replace('#SYNERGY2_PATH', syn2_path)
		s_file = s_file.replace('#WORKING_DIR', working_dir)
		s_file = s_file.replace('#CHILD1', child1)
		s_file = s_file.replace('#CHILD2', child2)
		s_file = s_file.replace('#NODE', curNode)
		s_file = s_file.replace('#ID', cmd_id)
		s_file = s_file.replace('#BLAST_EVAL', str(self.blast_eval))
		s_file = s_file.replace('#NUM_CORES', str(self.num_cores))
		s_file = s_file.replace('#HAMMING', str(self.hamming))
		s_file = s_file.replace('#ALPHA', str(self.alpha))
		s_file = s_file.replace('#BETA', str(self.beta))
		s_file = s_file.replace('#GAMMA', str(self.gamma))
		s_file = s_file.replace('#GAIN', str(self.gain))
		s_file = s_file.replace('#LOSS', str(self.loss))
		s_file = s_file.replace('#MIN_BEST_HIT', str(self.min_best_hit))
		s_file = s_file.replace('#MIN_SYNTENIC_FRACTION', str(self.min_syn_frac))
# 		s_file = s_file.replace('#HOMOLOGY_SCALE', str(self.homScale))
# 		s_file = s_file.replace('#SYNTENY_SCALE', str(self.synScale))
		s_file = s_file.replace('#WORKING_DIR', working_dir)

		my_sh = open(my_sh_file, 'w')
		my_sh.write(s_file)
		my_sh.close()

	# def toNewick(self, graph):
	# 	up = []  # unprocessed
	# 	leaf = []
	# 	for n in graph.nodes():
	# 		if len(graph[n]) > 1:
	# 			up.append(n)
	# 		else:
	# 			leaf.append((n, n[0]))
	# 	curNode = None
	# 	last_string = ""
	# 	if len(graph.nodes()) == 2:
	# 		ew = str(graph[leaf[0][0]][leaf[1][0]]['weight'])
	# 		last_string = "(" + leaf[0][0] + ":" + ew + "," + leaf[1][0] + ":" + ew + ")"
	# 	while len(up) > 0:
	# 		(curNode, e_count) = self.calcMostEdgesToLeaves(up, leaf, graph)
	# 		leaves = []
	# 		for e in graph[curNode]:
	# 			for l in leaf:
	# 				if l[0] == e:
	# 					e_i = leaf.index(l)
	# 					e_text = e
	# 					if 'child_newick' in graph.node[e]:
	# 						if e_count > 2 and len(up) > 1:
	# 							continue
	# 						e_text = graph.node[e]['child_newick']
	# 					leaf.pop(e_i)
	# 					ew = graph[curNode][e]['weight']
	# 					text = e_text + ":" + str(ew)
	# 					leaves.append(text)
	# 		# add newick text to curNode
	# 		node_text = "(" + ",".join(leaves) + ")"
	# 		last_string = node_text
	# 		graph.node[curNode]['child_newick'] = node_text
	# 		# change curNode to leaf
	# 		cn_i = up.index(curNode)
	# 		up.pop(cn_i)
	# 		leaf.append((curNode, curNode[0]))
	# 	if len(leaf) == 2 and len(up) == 0 and len(graph.nodes()) > 2:
	# 		ew = str(graph[leaf[0][0]][leaf[1][0]]['weight'])
	# 		last_string = "(" + graph.node[leaf[0][0]]['child_newick'] + ":" + ew + "," + graph.node[leaf[1][0]]['child_newick'] + ":" + ew + ")"
	# 	last_string = last_string.replace("(", "(\n")
	# 	last_string = last_string.replace(",", ",\n")
	# 	last_string = last_string.replace(")", ")\n")
	# 	last_string = last_string.rstrip()
	# 	return last_string + ";"

	# def calcMostEdgesToLeaves(self, unprocN, leaf, TG):
	# 	mostLeaves = -1
	# 	retNode = None
	# 	for n in unprocN:
	# 		e_count = 0
	# 		for e in TG[n]:
	# 			for l in leaf:
	# 				if e == l[0]:
	# 					e_count += 1
	# 		if e_count > mostLeaves:
	# 			mostLeaves = e_count
	# 			retNode = n
	# 	return (retNode, mostLeaves)

	def getGenomeToLocus(self):
		return self.genomeToLocus

	def getLocusToGenome(self):
		return self.locusToGenome
