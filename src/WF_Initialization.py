#!/usr/bin/env python

# Gets the filesystem set prepped to run the actual algorithm and writes the commands to do so
# import sys
import string
import random
import os
# import numpy
import pickle
# import re
import logging
# import TreeLib
import networkx as nx
import traceback
import collections


class Tree:
	logger = logging.getLogger("Tree")

	def __init__(self, tree_obj, flow_name, blast_eval, num_cores, alpha, beta, gamma, gain, loss, min_best_hit, syn_dist, homScale, synScale, numHits, minSynFrac, hamming):
		self.tree_obj = tree_obj
		self.genomeToLocusFile = tree_obj.genomeToLocusFile
		self.genomeToLocus = tree_obj.genomeToLocus
		self.locusToGenome = tree_obj.locusToGenome
		self.tree_string = ""
		self.tree = tree_obj.tree
		self.rooted_tree = tree_obj.rooted_tree
		self.blast_eval = blast_eval
		self.num_cores = num_cores
		self.alpha = alpha
		self.beta = beta
		self.gamma = gamma
		self.gain = gain
		self.loss = loss
		self.min_best_hit = min_best_hit
# 		self.cmds_per_job = cmds_per_job
		self.syn_dist = int(syn_dist)
		self.homScale = homScale
		self.synScale = synScale
		self.flow_name = flow_name
		self.num_hits = numHits
		self.min_syn_frac = minSynFrac
		self.hamming = hamming
		self.syn2_path = "#SYNERGY2_PATH"
		Tree.logger.debug("Tree initialized")

	def codeGenomeID(self, genome):
		Tree.logger.debug("".join(traceback.format_stack()))
		tag = ''
		if genome in self.genomeToLocus:
			tag = self.genomeToLocus[genome]
		else:
			# TODO CHECK IF THIS IS EVER CALLED
			Tree.logger.error("If this is called, check when and why to modify the tag generated.")
			import pdb
			pdb.set_trace()
			tag = ''.join(random.choice(string.ascii_uppercase) for x in range(3))
			while tag in self.locusToGenome:
				tag = ''.join(random.choice(string.ascii_uppercase) for x in range(3))
		self.genomeToLocus[genome] = tag
		self.locusToGenome[tag] = genome
		return tag

	def reregisterGenomeID(self, identifier, newChildren):
		oldGenome = self.locusToGenome[identifier]
		newGenome = ";".join(newChildren)
		del self.genomeToLocus[oldGenome]
		self.genomeToLocus[newGenome] = identifier
		self.locusToGenome[identifier] = newGenome

	def writeLocusTagFile(self):
		tag_out = open(self.genomeToLocusFile, 'w')
		for g in self.genomeToLocus:
			line = "\t".join([g, self.genomeToLocus[g]]) + "\n"
			tag_out.write(line)
		tag_out.close()
		Tree.logger.info("Wrote locus tags to locus_tag_file.txt")

	def makePicklesForSingleGenome(self, working_dir, genome, node):
		gdat = open(working_dir + "genomes/" + genome + "/annotation.txt", 'r').readlines()
		# this file might not be necessary in the long-term, in fact these genome directories may not be necessary at all
		# everything should be pulled out and stored in pickles after the repo files are parsed.  Whatever.
		ndat = open(working_dir + "nodes/" + node + "/" + node + ".pep", 'w')
		# print node
		x = gdat[0].split()[1]
		y = "_".join(x.split("_")[:-1])
		if not y == node:
			Tree.logger.warning("%s %s", (node, y))
			Tree.logger.warning("%s" % (gdat[0]))

		# these hashes will be pickles
		genes = {}
		gene_map = {}
		locusToTID = {}
		# this dict is for synteny information, keys are scaffold IDs
		neighbors = {}
		for g in gdat:
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
		root = nodeTier[max(nodeTier)][0]
		stack = []
		queue = collections.deque([[root, 0, []]])  # current_node, current_id, [child1_id, child2_id]
		count = 1
		while len(queue) > 0:
			current = queue.pop()
			stack.append(current)
			for e in self.rooted_tree.edges(current[0]):  # get children
				queue.append([e[1], count, []])  # add child to the queue (~breadth first search)
				current[2].append(count)  # add child to parent dependency
				count += 1
		with open(working_dir + "uger_jobs.sh", "w") as out:
			out.write("#! /bin/bash\n\nTIME=$(date +%s)\n")
			while len(stack) > 0:
				current = stack.pop()
				if current[0][0] != "L":  # not a leaf
					out.write("qsub -N j${TIME}" + str(current[1]))
					if len(current[2]) != 0:
						out.write(" -hold_jid j${TIME}" + str(current[2][0]) + ",j${TIME}" + str(current[2][1]))
					out.write(" " + working_dir + "nodes/" + str(current[0]) + "/" + str(current[0]) + ".sh\n")
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
					# if not (e[1] in nodeTier[0] or e[1] in local_proc_nodes):
						# print "WaitForFile ",e[1]
						# local_proc_nodes.append((next_cmd_id,"wait",e[1]))
						# set_cmd_count+=1
						# next_cmd_id = "1."+str(sets)+"."+str(set_cmd_count)

				# print "NodeFlow",curNode,kids
				self.makeSingleNodeFlow(working_dir, curNode, next_cmd_id, kids)
				# self.makeSingleConsensusFlow(working_dir,curNode,next_cmd_id)
				# self.makeFinalizeNodeFlow(working_dir,curNode,next_cmd_id)
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
		self.makeNodeFlowLauncher(working_dir, serial_sets)

	def makeSingleNodeFlow(self, working_dir, curNode, cmd_id, kids):
		my_dir = working_dir + "nodes/" + curNode + "/"
		child1 = kids[0]
		child2 = kids[1]

		syn2_path = self.syn2_path
# 		config_file = syn2_path + "WF_NodeFlowTemplate.ini"
		# config_file = syn2_path+"WF_NewNodeFlowTemplate.ini"
# 		template_file = syn2_path + "WF_NodeFlowTemplate.xml"
		sh_file = syn2_path + "NewNodeShTemplate.sh"

# 		my_config_file = my_dir + curNode + ".ini"
# 		my_template_file = my_dir + curNode + ".xml"
		my_sh_file = my_dir + curNode + ".sh"

# 		c_file = open(config_file, 'r').read()
# 		c_file = c_file.replace('#SYNERGY2_PATH', syn2_path)
# 		c_file = c_file.replace('#WORKING_DIR', working_dir)
# 		c_file = c_file.replace('#CHILD1', child1)
# 		c_file = c_file.replace('#CHILD2', child2)
# 		c_file = c_file.replace('#NODE', curNode)
# 		c_file = c_file.replace('#ID', cmd_id)
# 		c_file = c_file.replace('#BLAST_EVAL', str(self.blast_eval))
# 		c_file = c_file.replace('#NUM_CORES', str(self.num_cores))
# 		c_file = c_file.replace('#HAMMING', str(self.hamming))
# 		c_file = c_file.replace('#CMDS_PER_JOB', str(self.cmds_per_job))
# 		c_file = c_file.replace('#ALPHA', str(self.alpha))
# 		c_file = c_file.replace('#GAMMA', str(self.gamma))
# 		c_file = c_file.replace('#GAIN', str(self.gain))
# 		c_file = c_file.replace('#LOSS', str(self.loss))
# 		c_file = c_file.replace('#MIN_BEST_HIT', str(self.min_best_hit))
# 		c_file = c_file.replace('#NUM_HITS', str(self.num_hits))
# 		c_file = c_file.replace('#MIN_SYNTENIC_FRACTION', str(self.min_syn_frac))
# 		c_file = c_file.replace('#HOMOLOGY_SCALE', str(self.homScale))
# 		c_file = c_file.replace('#SYNTENY_SCALE', str(self.synScale))
# 		c_file = c_file.replace('#WORKING_DIR', working_dir)

# 		my_conf = open(my_config_file, 'w')
# 		my_conf.write(c_file)
# 		my_conf.close()

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
		s_file = s_file.replace('#HOMOLOGY_SCALE', str(self.homScale))
		s_file = s_file.replace('#SYNTENY_SCALE', str(self.synScale))
		s_file = s_file.replace('#WORKING_DIR', working_dir)

		my_sh = open(my_sh_file, 'w')
		my_sh.write(s_file)
		my_sh.close()

# 		t_file = open(template_file, 'r').read()
# 		t_file = t_file.replace('#NODE', curNode)
# 		t_file = t_file.replace('#ID', cmd_id)
# 		t_file = t_file.replace('#SYNERGY2_PATH', syn2_path)
# 		t_file = t_file.replace('#WORKING_DIR', working_dir)
# 		t_file = t_file.replace('#CHILD1', child1)
# 		t_file = t_file.replace('#CHILD2', child2)
# 		t_file = t_file.replace('#WORKING_DIR', working_dir)
# 		my_temp = open(my_template_file, 'w')
# 		my_temp.write(t_file)
# 		my_temp.close()

	# def makeSingleConsensusFlow(self, working_dir, curNode, prev_cmd_id):
	# 	ids = prev_cmd_id.split(".")
	# 	last = int(ids[-1])+1
	# 	ids[-1] = str(last)
	# 	cmd_id = ".".join(ids)

	# 	s2 = self.syn2_path
	# 	my_dir = working_dir+"nodes/"+curNode+"/"

	# 	config_file = s2+"WF_consensus_seq_flow_template.ini"
	# 	template_file = s2+"WF_consensus_seq_flow_template.xml"

	# 	my_config_file = my_dir+"consensus_seq_flow_config.ini"
	# 	my_template_file = my_dir+"consensus_seq_flow_template.xml"

	# 	c_file = open(config_file, 'r').read()
	# 	c_file = c_file.replace('#WORKING_DIR', working_dir)
	# 	c_file = c_file.replace('#NODE', curNode)
	# 	c_file = c_file.replace('#ID', cmd_id)
	# 	c_file = c_file.replace('#SYNERGY2_PATH', s2)

	# 	t_file = open(template_file, 'r').read()
	# 	t_file = t_file.replace('#ID', cmd_id)

	# 	my_conf = open(my_config_file, 'w')
	# 	my_conf.write(c_file)
	# 	my_conf.close()

	# 	my_temp = open(my_template_file, 'w')
	# 	my_temp.write(t_file)
	# 	my_temp.close()

	# def makeFinalizeNodeFlow(self, working_dir, curNode, prev_cmd_id):
	# 	ids = prev_cmd_id.split(".")
	# 	last = int(ids[-1])+2
	# 	ids[-1] = str(last)
	# 	cmd_id = ".".join(ids)

	# 	s2 = self.syn2_path
	# 	my_dir = working_dir+"nodes/"+curNode+"/"

	# 	config_file = s2+"WF_finalize_node_template.ini"
	# 	template_file = s2+"WF_finalize_node_template.xml"

	# 	my_config_file = my_dir+"finalize_node_config.ini"
	# 	my_template_file = my_dir+"finalize_node_template.xml"

	# 	c_file = open(config_file, 'r').read()
	# 	c_file = c_file.replace('#WORKING_DIR', working_dir)
	# 	c_file = c_file.replace('#NODE', curNode)
	# 	c_file = c_file.replace('#ID', cmd_id)
	# 	c_file = c_file.replace('#SYNERGY2_PATH', s2)

	# 	t_file = open(template_file, 'r').read()
	# 	t_file = t_file.replace('#ID', cmd_id)

	# 	my_conf = open(my_config_file, 'w')
	# 	my_conf.write(c_file)
	# 	my_conf.close()

	# 	my_temp = open(my_template_file, 'w')
	# 	my_temp.write(t_file)
	# 	my_temp.close()

	def makeNodeFlowLauncher(self, working_dir, serial_sets):
		ini_file = open(working_dir + self.flow_name + ".ini", 'w')
		xml_file = open(working_dir + self.flow_name + ".xml", 'w')
		xml_file.write("""<?xml version="1.0" encoding="utf-8"?>""" + "\n")
		xml_file.write("""<commandSetRoot>""" + "\n")
		xml_file.write("""\t<commandSet type="parallel">""" + "\n")
		xml_file.write("""\t\t<name>""" + self.flow_name + """</name>""" + "\n")
		ini_file.write("[1]\nname=" + self.flow_name + "\n")
		instance_commands = []

		for s in serial_sets:
			sub_s = 1
			xml_file.write("""\t\t<commandSet type="serial">""" + "\n")
			xml_file.write("""\t\t\t<configMapId>1.""" + str(s) + """</configMapId>""" + "\n")
			xml_file.write("""\t\t\t<name>subflow1.""" + str(s) + """</name>""" + "\n")
			ini_file.write("[1." + str(s) + "]\n")
			ini_file.write("name=subflow1." + str(s) + "\n\n")

			for n in serial_sets[s]:
				ini_id = n[0]
				if len(n) == 3 and n[1] == "wait":
					# wait command
					xml_file.write("\t\t\t<command>\n")
					xml_file.write("\t\t\t\t<name>WaitForNode</name>\n")
					xml_file.write("\t\t\t\t<configMapId>" + ini_id + "</configMapId>\n")
					xml_file.write("\t\t\t</command>\n")
					ini_file.write("[" + ini_id + "]\n")
					ini_file.write("name=WaitForNode\n")
					ini_file.write("type=RunUnixCommand\n")
					ini_file.write("executable=" + self.syn2_path + "WF_WaitForFile.py\n")
					ini_file.write("arg=" + working_dir + "nodes/" + n[2] + "/\n")
					ini_file.write("arg=NODE_COMPLETE\n\n")
				else:
					# node flow command
					node = n[1]
					path_to_node_dir = working_dir + "nodes/" + node + "/"
					instance_cmd = "RunWorkflow -t " + path_to_node_dir + node + ".xml -c " + path_to_node_dir + node + ".ini -i " + path_to_node_dir + node + "_instance.xml\n"
					# tier = int(ini_id.split(".")[-1])
					dist_to_root = len(nx.shortest_path(self.rooted_tree, self.tree_obj.root, node))
					instance_commands.append((instance_cmd, dist_to_root))
					xml_file.write("\t\t\t<command>\n")
					xml_file.write("\t\t\t\t<name>" + node + "</name>\n")
					xml_file.write("\t\t\t\t<configMapId>" + ini_id + "</configMapId>\n")
					xml_file.write("\t\t\t</command>\n")
					ini_file.write("[" + ini_id + "]\n")
					ini_file.write("name=" + node + "\n")
					ini_file.write("type=RunUnixCommand\n")
					ini_file.write("executable=" + self.syn2_path + "WF_WaitForFile.py\n")
					ini_file.write("arg=" + working_dir + "nodes/" + n[1] + "/\n")
					ini_file.write("arg=NODE_COMPLETE\n\n")

				sub_s += 1
			xml_file.write("""\t\t</commandSet>""" + "\n")

		xml_file.write("""\t</commandSet>""" + "\n")
		xml_file.write("""</commandSetRoot>""")
		xml_file.close()

		instance_command_out = open("instance_commands.txt", 'w')
		instance_commands = sorted(instance_commands, key=lambda tup: tup[1], reverse=True)
		for ic in instance_commands:
			instance_command_out.write(ic[0])
		instance_command_out.close()

	def toNewick(self, graph):
		up = []  # unprocessed
		leaf = []
		for n in graph.nodes():
			if len(graph[n]) > 1:
				up.append(n)
			else:
				leaf.append((n, n[0]))
		curNode = None
		last_string = ""
		if len(graph.nodes()) == 2:
			ew = str(graph[leaf[0][0]][leaf[1][0]]['weight'])
			last_string = "(" + leaf[0][0] + ":" + ew + "," + leaf[1][0] + ":" + ew + ")"
		while len(up) > 0:
			(curNode, e_count) = self.calcMostEdgesToLeaves(up, leaf, graph)
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
						ew = graph[curNode][e]['weight']
						text = e_text + ":" + str(ew)
						leaves.append(text)
			# add newick text to curNode
			node_text = "(" + ",".join(leaves) + ")"
			last_string = node_text
			graph.node[curNode]['child_newick'] = node_text
			# change curNode to leaf
			cn_i = up.index(curNode)
			up.pop(cn_i)
			leaf.append((curNode, curNode[0]))
		if len(leaf) == 2 and len(up) == 0 and len(graph.nodes()) > 2:
			ew = str(graph[leaf[0][0]][leaf[1][0]]['weight'])
			last_string = "(" + graph.node[leaf[0][0]]['child_newick'] + ":" + ew + "," + graph.node[leaf[1][0]]['child_newick'] + ":" + ew + ")"
		last_string = last_string.replace("(", "(\n")
		last_string = last_string.replace(",", ",\n")
		last_string = last_string.replace(")", ")\n")
		last_string = last_string.rstrip()
		return last_string + ";"

	def calcMostEdgesToLeaves(self, unprocN, leaf, TG):
		mostLeaves = -1
		retNode = None
		for n in unprocN:
			e_count = 0
			for e in TG[n]:
				for l in leaf:
					if e == l[0]:
						e_count += 1
			if e_count > mostLeaves:
				mostLeaves = e_count
				retNode = n
		return (retNode, mostLeaves)

# 	def getMinDistancePair(self, nodes, extinct):
# 		min_extinct_pcount = sorted(extinct, key=lambda tup: tup[2])[0][2]
# 		# calc smallest line distance between 2 nodes
# 		# (ln_count,locus,dist,p_count,tp_count, c_paren)
# 		min_dist = nodes[-1][0] - nodes[0][0] + 1.0
# 		i = len(nodes) - 1
# 		while i >= 0:
# 			j = len(nodes) - 1
# 			while j >= 0:
# 				min_row = min(i, j)
# 				max_row = max(i, j)
# 				if i == j:
# 					j -= 1
# 				elif nodes[min_row][5] is True:
# 					j -= 1
# 				elif min(nodes[i][3], nodes[j][3]) < min_extinct_pcount:
# 					j -= 1
# 				elif nodes[min_row][3] <= nodes[max_row][3]:
# 					j -= 1
# 				else:
# 					if abs(nodes[i][0] - nodes[j][0]) < min_dist:
# 						min_dist = abs(nodes[i][0] - nodes[j][0])
# 					j -= 1
# 			i -= 1
# 		# if a node pairing has the minimum distance, add it
# 		md_pairs = []
# 		max_rparen = 0
# 		for n in nodes:
# 			for m in nodes:
# 				if m[0] <= n[0]:
# 					continue
# 				min_row = min(m[0], n[0])
# 				if m[0] == min_row and m[5] is True:
# 					continue
# 				if n[0] == min_row and n[5] is True:
# 					continue
# 				if abs(n[3]-m[3]) > 1:
# 					continue
# 				if abs(n[0]-m[0]) == min_dist:
# 					my_mrp = max(n[3], m[3])
# 					if my_mrp > max_rparen:
# 						max_rparen = my_mrp
# 					md_pairs.append((n, m))
# 		# get pairs that have the maximum rparen counts
# 		max_rparen_pairs = []
# 		for md in md_pairs:
# 			if md[0][3] == max_rparen or md[1][3] == max_rparen:
# 				max_rparen_pairs.append(md)
# 		max_rparen_pairs = sorted(max_rparen_pairs, key=lambda pair: pair[1][4], reverse=True)
# 		return max_rparen_pairs[0]

	def getGenomeToLocus(self):
		return self.genomeToLocus

	def getLocusToGenome(self):
		return self.locusToGenome
