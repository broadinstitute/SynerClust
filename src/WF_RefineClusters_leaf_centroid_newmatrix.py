#!/usr/bin/env python

import sys
import os
import NJ
import pickle
import networkx as nx
import numpy
from collections import Counter
import logging
import subprocess
# import scipy.spatial.distance as distance

DEVNULL = open(os.devnull, 'w')


def usage():
	print """From the rough cluster trees generated by WF_MakeRoughClusters, refine clusters so that all duplication events (paralogs) \
	occur after speciation from the most recent common ancestor or MRCA, which is [node]. [node_dir] is the directory that contains all \
	data about [node]. [flow_id] refers to a step in the workflow; [jobs_per_cmd] is the number of consensus sequence computations \
	distributed to a single node.

	[alpha], [gamma], [gain], and [loss] are parameters that impact the rooting of a tree

	WF_RefineClusters.py [node_dir] [flow_id] [jobs_per_cmd] [node] [alpha] [gamma] [gain] [loss] [children...]
	"""
	sys.exit(1)


# def load_leaves(leavespkls, locus_tag, node):
# 	for leaf in locus_tag[node]:
# 		if leaf[:9] == "L_0000000":
# 			leavespkls[leaf] = None
# 		else:
# 			load_leaves(leavespkls, locus_tag, leaf)


def main(argv):
	node_dir = argv[0]
	repo_path = argv[0][:-6]
	mrca = argv[3]
	node = argv[3]
	alpha = float(argv[4])
	gamma = float(argv[5])
	gain = float(argv[6])
	loss = float(argv[7])
	children = argv[8:]
	my_dir = node_dir + mrca + "/"
# 	NO_BREAK_EW = 0.5
	beta = 0.01

	FORMAT = "%(asctime)-15s %(levelname)s %(module)s.%(name)s.%(funcName)s at %(lineno)d :\n\t%(message)s\n"
	logger = logging.getLogger()
	logging.basicConfig(filename=my_dir + 'RefineClusters_leaf_centroid.log', format=FORMAT, filemode='w', level=logging.DEBUG)
	# add a new Handler to print all INFO and above messages to stdout
	ch = logging.StreamHandler(sys.stdout)
	ch.setLevel(logging.INFO)
	logger.addHandler(ch)
	logger.info('Started')

	if "CLUSTERS_REFINED" in os.listdir(my_dir):
		sys.exit(0)

	# read trees, resolve clusters
	tree_dir = my_dir + "trees/"
	cluster_dir = my_dir + "clusters"
	if "clusters" in os.listdir(my_dir):
		if "old" not in os.listdir(my_dir):
			os.system("mkdir " + my_dir + "old")

		os.system("mv -f " + cluster_dir + "/ " + my_dir + "old/")
	os.system("mkdir " + cluster_dir)
	cluster_dir = cluster_dir + "/"

	cluster_counter = 1  # used to number clusters
	synteny_data = {}
	pickleSeqs = {}
	pickleToCons = {}
	singletons_pep = {}
	pickleMaps = {}
	picklePeps = {}
	childrenpkls = {}
	# print "last_tree", last_tree
	# load locus_mapping files from children
	for c in children:
		mapFile = node_dir + c + "/locus_mappings.pkl"
		pklFile = open(mapFile, 'rb')
		pickleMaps[c] = pickle.load(pklFile)
		pklFile.close()
		synFile = node_dir + c + "/synteny_data.pkl"
		pklFile = open(synFile, 'rb')
		synteny_data[c] = pickle.load(pklFile)
		pklFile.close()
		if c[0] == "L":
			with open(node_dir + c + "/" + c + ".pkl", "r") as f:
				childrenpkls[c] = pickle.load(f)
		else:  # c[0] == "N"
			with open(node_dir + c + "/pep_data.pkl", "r") as f:
				childrenpkls[c] = pickle.load(f)
			with open(node_dir + c + "/singletons_pep_data.pkl", "r") as f:
				childrenpkls[c].update(pickle.load(f))
		
# 	old_orphans = open(tree_dir + "orphan_genes.txt", 'r').readlines()
# 	orphans = open(tree_dir + "orphan_genes.txt", 'r').readlines()
	orphans = []
	ok_trees = []
# 
	locus_tag = {}
	with open(repo_path + "genomes/locus_tag_file.txt", "r") as f:
		for line in f:
			if ";" in line:
				s = line.rstrip().split("\t")
				locus_tag[s[1]] = s[0].split(";")

# 	load_leaves(childrenpkls, locus_tag, node)

	pklfile = None
# 	syntable = {}
# 	for leaf in childrenpkls.keys():
# 		with open(repo_path + "nodes/" + leaf + "/" + leaf + ".pkl", "r") as pklfile:
# 			childrenpkls.update(pickle.load(pklfile))
# 		with open(repo_path + "nodes/" + leaf + "/synteny_table.pkl", "r") as syntablefile:
# 			syntable[leaf] = pickle.load(syntablefile)

# 	with open(repo_path + "nodes/" + node + "/locus_mappings.pkl", 'r') as pklfile:
# 		locus_mapping = pickle.load(pklfile)

# 	cluster_leaf_count = {}
# 	gene_to_cluster = {}
# 	for cluster in locus_mapping:
# 		cluster_leaf_count[cluster] = {}
# 		for gene in locus_mapping[cluster]:
# 			leaf = "_".join(gene.split("_")[:-1])
# 			gene_to_cluster[gene] = cluster
# 			if leaf not in cluster_leaf_count[cluster]:
# 				cluster_leaf_count[cluster][leaf] = 0
# 			cluster_leaf_count[cluster][leaf] += 1
	with open(repo_path + "nodes/" + node + "/trees/gene_to_cluster.pkl", "r") as f:
		gene_to_cluster = pickle.load(f)
	with open(repo_path + "nodes/" + node + "/trees/cluster_to_genes.pkl", "r") as f:
		cluster_to_genes = pickle.load(f)
# 	with open(repo_path + "nodes/" + node + "pep_data.pkl", "r") as f:
# 		pickleSeqs = pickle.load(f)

	# fasttree_cmd =["#MUSCLE_PATH", "-maxiters", "2", "-diags", "-sv", "-distance1", "kbit20_3", "-quiet", "|", "#FASTTREE_PATH", "-quiet"]
	muscle_cmd = ["/home/kamigiri/tools/muscle3.8.31_i86linux64", "-maxiters", "2", "-diags", "-sv", "-distance1", "kbit20_3", "-quiet"]
	fasttree_cmd = ["/home/kamigiri/tools/FastTreeDouble", "-quiet", "-nosupport"]
	ok_trees = []
	
	for cluster in cluster_to_genes:
# 	for clusterID in pickleSeqs:
		if len(cluster_to_genes[cluster]) == 1:
			orphans.append(cluster_to_genes[cluster][0])
			continue
# 		if len(pickleSeqs[clusterID]) == 1:
# 			orphans.append(pickleSeqs[clusterID][0].split(";")[0][1:])
# 			continue
		
			
		stdin_data = ""
# 		add_to_stdin = ""
# 		for gene in cluster_to_genes[cluster]:
# 			stdin_data += ">" + gene + "\n"
# 			stdin_data += childrenpkls[gene] + "\n"
		for gene in cluster_to_genes[cluster]:
# 			stdin_data += ">" + gene + "\n"
			stdin_data += ">" + gene + "\n"
			try:
				if children[0][0] == "L":
					stdin_data += childrenpkls[children[0]][gene] + "\n"
				else:
					stdin_data += childrenpkls[children[0]][gene][0].split("\n")[1] + "\n"					
			except KeyError:
				if children[1][0] == "L":
					stdin_data += childrenpkls[children[1]][gene] + "\n"
				else:
					stdin_data += childrenpkls[children[1]][gene][0].split("\n")[1] + "\n"
			
		process = subprocess.Popen(muscle_cmd, stdin = subprocess.PIPE, stdout = subprocess.PIPE, stderr = DEVNULL)
		output = process.communicate(stdin_data)[0]
		process = subprocess.Popen(fasttree_cmd, stdin = subprocess.PIPE, stdout = subprocess.PIPE, stderr = DEVNULL)
		output = process.communicate(output)[0]
		# if (len(locus_mapping[cluster]) > 5):
		print(output + "\n\n")
		# tree = NJ.NJTree("", "", node, alpha, beta, gamma, gain, loss)

		graph = nx.Graph()
		counter = 1
		leaves = []
		while(True):
			r = output.find(")")
			l = output[:r].rfind("(")

			children_string = output[l + 1:r].split(",")
			if len(children_string) == 1:
				if len(graph.nodes()) == 0:
					ok_trees.append(output.split(":")[0][1:])
				break
			group = "node" + str(counter)
			counter += 1

			if group not in graph.nodes():  # isn't it always a new node?
				graph.add_node(group)
			for child in children_string:
				child = child.split(":")
				if child[0] not in graph.nodes():
					graph.add_node(child[0])
					leaves.append(child[0])
				graph.add_edge(group, child[0], dist=float(child[1]))
			output = output[:l] + group + output[r + 1:]

		assert leaves.sort() == cluster_to_genes[cluster].sort()  # if good, can replace leaves list by graph.nodes()
# 		if 'L_0000000_9VzQfBUR7laObXloLMne2g_000709' not in graph.nodes():
# 			continue
# 		if 'L_0000000_9VzQfBUR7laObXloLMne2g_001369' not in graph.nodes():
# 			continue
		# get distance matrix from this tree, then add synteny distances to it with a certain ratio
# 		matrix = {}
# 		nodes = graph.nodes()  # nodes are headers/gene references
# 		nodes.sort()
		leaves.sort()

		syn = {}
		for n in leaves:  # genes
# 		for n in children:
			syn[n] = []
			leaf = "_".join(n.split("_")[:-1])
			for m in synteny_data[leaf][n]['neighbors']:
				syn[n].append(gene_to_cluster[m])

		hom_matrix = numpy.empty(len(leaves) * (len(leaves) - 1) / 2)
		syn_matrix = numpy.empty(len(leaves) * (len(leaves) - 1) / 2)
		i = 1
		pos = 0
		for m in leaves[1:]:
			syn_m = set(syn[m])
			mSeqs = len(syn[m])
# 			for n in graph.nodes():
			for n in leaves[:i]:
				hom_matrix[pos] = nx.shortest_path_length(graph, n, m, "dist")
				nSeqs = len(syn[n])
				matches = 0
				if mSeqs == 0 or nSeqs == 0:
					syn_matrix[pos] = 1.0  # no neighbors in common if someone has no neighbors  # -= 0 ? does it change anything?
					pos += 1
					continue
				all_neighbors = syn_m & set(syn[n])
				for a in all_neighbors:
					t_m = max(syn[m].count(a), 0)
					t_n = max(syn[n].count(a), 0)
					matches += min(t_m, t_n)
				synFrac = float(matches) / float(min(mSeqs, nSeqs))  # why mSeqs and not len(syn_m) which is a set that removes duplicates?
				syn_matrix[pos] = 1.0 - synFrac
				pos += 1
			i += 1
	# Root, evaluate and split every tree until all trees are OK
# 	hom_file = open(tree_dir + "homology_matrices.dat", "r")
# 	syn_file = open(tree_dir + "synteny_matrices.dat", "r")
# 	hom_line = hom_file.readline()
# 	syn_line = syn_file.readline()
# 	hom_mat = []
# 	syn_mat = []
# 	while hom_line != "":
# 		if hom_line == "//\n":
# 		while(True):
		unchecked_trees = []
		# Could probably replace these lines by directly reading the distance matrix file into the list without needing to create a NJ.NJTree
# 			myTree = NJ.NJTree(hom_mat, syn_mat, mrca, alpha, beta, gamma, gain, loss)
		myTree = NJ.NJTree(mrca, alpha, beta, gamma, gain, loss)
		myTree.buildGraphFromNewDistanceMatrix(hom_matrix, syn_matrix, leaves)
# 		myTree = NJ.NJTree(tree_file, syn_file, mrca, alpha, beta, gamma, gain, loss)
# 			d_lines = myTree.readDistanceMatrix()  # TODO this is actually just hom_mat and the only place it is used
# 			unchecked_trees.append((d_lines, "false"))  # False refers to orphan status
		# If "false" refers to orphan status, why are there also both "true" and "orphan"?
		unchecked_trees.append((myTree, False))
		while len(unchecked_trees) > 0:
			this_tree = unchecked_trees.pop()
			myTree = this_tree[0]
			isOrphan = this_tree[1]
# 				myTree = NJ.NJTree("filler", syn_mat, mrca, alpha, beta, gamma, gain, loss)
			# single node tree genes are added to orphans
			if isOrphan == "orphan":
				logger.critical("Need to handle orphan case")
# 				for n in nodes:  # can there be more than one node if it's an orphan? If not, can use nodes[0] and assert len(nodes) == 0
# 					o = n[1].split(":")[0]
# 					o = o.replace("(", "")
# 					orphans.append(o)
					# last_tree = "orphan"
				continue
			# multiple node trees continue
			else:
# 					bigNode = myTree.buildGraphFromDistanceMatrix(uTree)
				myleaves = myTree.bigNode.split(";")
				mysources = set([])  # sources are the child species contributing to this tree
				for m in myleaves:
					mysources.add("_".join(m.split("_")[:-1]))
				# a valid tree has genes from both children, single source trees are broken into individual genes and added to the orphan list
				if len(mysources) == 1:
					for m in myleaves:
						orphans.append(m)
				# if the tree has >1 source, it is rooted and evaluated
				else:
					root = myTree.rootTree()
					myTree.checkTree(root)
					# tree is valid, added to resolved clusters
					if myTree.OK == "true":
						format_nodes = []
						for n in myTree.graph.nodes():
							if n.count(";") == 0:
								format_nodes.append(n)
						ok_trees.append(format_nodes)
					# tree is invalid, added to unchecked trees unless it is an orphan
					else:
						# additional orphan exit
						if myTree.OK == "orphan":
							unchecked_trees.append((NJ.NJTree.toNewick(myTree.graph).split("\n"), myTree.OK))
# 								unchecked_trees.append((NJ.NJTree.splitNewTree(myTree), myTree.OK)) # need to return both subtrees + myTree.OK in list
						else:
# 								(myNewicks, myMatrices) = myTree.splitTree(root)
							(new_trees, new_root_edges) = myTree.splitNewTree(root)
# 								for m in myMatrices:
							for new_tree in new_trees:
								unchecked_trees.append((new_tree, myTree.OK))
# 			hom_mat = []
# 			syn_mat = []
# 		
# 			hom_mat.append(hom_line)
# 			syn_mat.append(syn_line)
# 		hom_line = hom_file.readline()
# 		syn_line = syn_file.readline()
# # 	for t in os.listdir(tree_dir):
# 	hom_file.close()
# 	syn_file.close()


	newPickleMap = {}  # this will turn into the locus mappings for this node
	newSyntenyMap = {}
	childToCluster = {}  # child gene/og --> og.id for this node

# 	for o in old_orphans:
# 		orphans.append(o)

	for o in orphans:
		ok_trees.insert(0, [o.rstrip()])

	blast_pep = {}
	for c in children:
		my_blast_pep = open(my_dir + c + ".blast.fa", 'r').readlines()
		curBlast = ""
		curPep = ""
		for m in my_blast_pep:
			m = m.rstrip()
			if len(m) < 1:
				continue
			if m.find(">") > -1:
				if len(curBlast) > 0:
					blast_pep[curBlast].append(curPep)
					curPep = ""
				line = m[1:]
				curBlast = line.split(";")[0]
				if curBlast not in blast_pep:
					blast_pep[curBlast] = []
			else:
				curPep += m
				# blast_pep[curBlast] += m
		if len(curPep) > 0:
			blast_pep[curBlast].append(curPep)

	# sys.exit()
	# make control files for consensus sequence formation
	# cons_cmds = []
	singletons = cluster_dir + "singletons.cons.pep"
	singles = open(singletons, 'w')
	sum_stats = my_dir + "summary_stats.txt"
	sstats = open(sum_stats, 'w')
	for ok in ok_trees:
		c = str(cluster_counter)
		clusterID = ""
		while len(c) < 6:
			c = "0" + c
		clusterID = mrca + "_" + c
		newPickleMap[clusterID] = []
		newSyntenyMap[clusterID] = {'count': 0, 'neighbors': [], 'children': []}

		# get list of leaf sequences to pull and organize in treeSeqs
		treeSeqs = {}
		tree_seq_count = 0
		# leafSeqs = {}
		child_leaves = {}
		taxa = set([])
		taxa_map = {}
		for g in ok:
			child = "_".join(g.split("_")[:-1])
			newSyntenyMap[clusterID]['children'].append(g)
			childToCluster[g] = clusterID
			leafKids = pickleMaps[child][g]
			if child not in child_leaves:
				child_leaves[child] = 0
			child_leaves[child] += len(leafKids)
# TODO in this part, children pickle are opened every time, could probably only open each one
			for l in leafKids:
				newPickleMap[clusterID].append(l)
				lKid = "_".join(l.split("_")[:-1])
				taxa.add(lKid)
				if lKid not in taxa_map:
					taxa_map[lKid] = 0
				taxa_map[lKid] += 1
				if lKid not in pickleSeqs:
					seqFile = node_dir + lKid + "/" + lKid + ".pkl"
					pklFile = open(seqFile, 'rb')
					pickleSeqs[lKid] = pickle.load(pklFile)
					pklFile.close()
				seq = pickleSeqs[lKid][l]
				if seq not in treeSeqs:
					treeSeqs[seq] = []
				treeSeqs[seq].append(l)
				tree_seq_count += 1
				newSyntenyMap[clusterID]['count'] += 1

		my_lengths = []
		min_taxa = len(taxa)
		max_taxa = 0
		for tm in taxa_map:
			tm_val = taxa_map[tm]
			if tm_val > max_taxa:
				max_taxa = tm_val
			if tm_val < min_taxa:
				min_taxa = tm_val
		for seq in treeSeqs:
			for me in treeSeqs[seq]:
				my_lengths.append(len(seq))
		avg = numpy.average(my_lengths)
		std = numpy.std(my_lengths)
		std_avg = std / avg
		out_dat = [clusterID, str(len(my_lengths)), str(len(taxa)), str(min_taxa), str(max_taxa), str(min(my_lengths)), str(max(my_lengths)), str(avg), str(std), str(std_avg)]

# 		sys.exit()
		sstats.write("\t".join(out_dat) + "\n")

		if tree_seq_count == 1:
			for seq in treeSeqs:
				seqlen = str(len(seq))
				singles.write(">" + clusterID + ";" + seqlen + "\n" + seq + "\n")
			seq = treeSeqs.keys()[0]
			singletons_pep[clusterID] = [">" + clusterID + ";" + str(len(seq)) + "\n" + seq + "\n"]
		elif len(ok) == 1:
			for bseq in blast_pep[ok[0]]:
				seqlen = str(len(bseq))
				singles.write(">" + clusterID + ";" + seqlen + "\n" + bseq + "\n")
			bseq = blast_pep[ok[0]][0]
			singletons_pep[clusterID] = [">" + clusterID + ";" + str(len(bseq)) + "\n" + bseq + "\n"]
		else:
			#  temp_pep = cluster_dir + clusterID + ".pep"
			#  pepOut = open(temp_pep, 'w')
			pickleToCons[clusterID] = []
			for seq in treeSeqs:
				seqlen = str(len(seq))
				identifier = treeSeqs[seq][0]  # TODO output ALL IDs from seq because there might be more than a single ID, and this is NOT a .cons.pep file, just a .pep file
# 				pepOut.write(">" + identifier + ";" + seqlen + "\n" + seq + "\n")
				pickleToCons[clusterID].append(">" + identifier + ";" + seqlen + "\n" + seq + "\n")
# 			pepOut.close()
		cluster_counter += 1
	singles.close()
	sstats.close()

	pklPep = my_dir + "pep_data.pkl"
	sdat = open(pklPep, 'wb')
	pickle.dump(pickleToCons, sdat)
	sdat.close()

	with open(my_dir + "singletons_pep_data.pkl", "w") as f:
		pickle.dump(singletons_pep, f)

	# update synteny data
	for clust in newSyntenyMap:
		for child in newSyntenyMap[clust]['children']:
			lc = "_".join(child.split("_")[:-1])
			# logger.debug("%s splitted to %s" % (child, lc))
			for neigh in synteny_data[lc][child]['neighbors']:
				# logger.debug("newSyntenyMap[%s]['neighbors'].append(childToCluster[%s]" % (clust, neigh))
				newSyntenyMap[clust]['neighbors'].append(childToCluster[neigh])
	# pickle synteny data
	pklSyn = my_dir + "synteny_data.pkl"
	sdat = open(pklSyn, 'wb')
	pickle.dump(newSyntenyMap, sdat)
	sdat.close()

	# pickle the locus mappings
	pklMap = my_dir + "locus_mappings.pkl"
	sdat = open(pklMap, 'wb')
	pickle.dump(newPickleMap, sdat)
	sdat.close()

	# script complete call
	clusters_done_file = my_dir + "CLUSTERS_REFINED"
	cr = open(clusters_done_file, 'w')
	cr.write("Way to go!\n")
	cr.close()

if __name__ == "__main__":
	if len(sys.argv) == 1:
		usage()
	else:
		main(sys.argv[1:])
