#!/usr/bin/env python

import sys
import os
# import getopt
import argparse
import logging
import COBRA_Repo_Handling
import WF_Initialization
import TreeLib
import networkx as nx

def usage():
	print """
	python WF_runSynergy2.py [OPTIONS]
	
	-h, --help
		Prints this usage and exits
	REQUIRED:
	-r, --repo [file]
		where [file] is the complete path to the data repository containing your genomic data
	-w, --working [path]
		where [path] is the complete path to the working directory for this analysis
	-t, --tree [file]
		where [file] is a species tree relating all of the genomes to be analyzed
	-f, --flow [string]
		where [string] is any string you choose, and will be the template of your Workflow instance
		
	OPTIONAL ALGORITHM VARIABLES:
	-g, --gamma [float]
		where [float] is a factor in the rooting equation, default = 10.0
	-G, --gain [float]
		where [float] is the probability of a gain event occurring, range (0.0, 1.0], default = 0.05
	-L, --loss [float]
		where [float] is the probability of a loss event occurring, range (0.0, 1.0], default = 0.05
	-b, --blast_eval [float]
		where [float] is the e-value used for all blast analysis, default = #BLAST_EVAL_DEFAULT
	-m, --min_best_hit [float]
		where [float] is the minimum score a blast hit must have to define an edge in the initial clusters, range (0.0, 1.0], default = 0.5
	-s, --synteny_window [int]
		where [int] is the distance in base pairs that will contribute to upstream and downstream to syntenic fraction.  The total window
		size is [int]*2. default = 6000
	-S, --synteny_scale [float]
		where [float] represents the contribution of synteny to branch weight for cluster trees, range (0.0, 1.0], default = 1.0
	-H, --homology_scale [float]
		where [float] represents the contribution of homology to branch weight for cluster trees, range (0.0, 1.0], default = 1.0
	-T, --top_hits [int]
		where [int] is the number of blast hits that will be considered for each gene, default=5
	-F, --min_syntenic_fraction [float]
		where [float] is the minimum syntenic fraction required for two genes from the same species to be considered paralogs
		range [0.0,1.0], default=0.5
	-l, --locus [file]
		where [file] is a locus_tag_file.txt that corresponds to the data in this repository
	-D, --hamming [float]
		where [float] is the maximum hamming distance between a representative sequence and a sequence being represented
		range [0.0,1.0], default=0.6
	
	OPTIONAL PERFORMANCE VARIABLES:
	-n, --num_cores [int]
		where [int] is the number of cores used for blast analysis (-a flag), default = #NUM_CORES_DEFAULT

	"""
	
def main(argv):
	FORMAT = "%(asctime)-15s %(levelname)s %(module)s.%(name)s.%(funcName)s at %(lineno)d :\n\t%(message)s\n"
	logger = logging.getLogger()
	logging.basicConfig(filename='runSynergy2.log', format = FORMAT, filemode='w', level=logging.DEBUG)
	# add a new Handler to print all INFO and above messages to stdout
	ch = logging.StreamHandler(sys.stdout)
	ch.setLevel(logging.INFO)
	logger.addHandler(ch)
	logger.info('Started')

# 	cobra_repo = ""
# 	cobra_repo_path = ""
# 	working_dir = ""
# 	species_tree = ""
# 	flow_name = ""
# 	genome_dir = ""	
# 	blast_eval = #BLAST_EVAL_DEFAULT
# 	num_cores = #NUM_CORES_DEFAULT
# 	cmds_per_job = 500
# 	synteny_window = 6000
# 	alpha = 10.0
# 	beta = 10.0
# 	gamma = 10.0
# 	gain = 0.05
# 	loss = 0.05
# 	min_best_hit = 0.5 
# 	homScale = 1.0
# 	synScale = 1.0
# 	numHits = 5
# 	minSynFrac = 0.5
# 	hamming = 0.6
# 	locus_file = ""	
# 	try:
# 		opts, args = getopt.getopt(argv, "r:w:t:f:g:G:L:l:m:b:c:n:s:H:S:l:T:F:D:h",["repo=","working=","tree=","flow=","gamma=","gain=","loss=","min_best_hit=","blast_eval=","cmds_per_job=","num_cores=","synteny_window=","homology_scale=","synteny_scale=","locus=","top_hits=","min_syntenic_fraction=","hamming=","help"])
# 
# 	except getopt.GetoptError:
# 		usage()
# 		sys.exit(2)

	distribute = 1

	usage = "usage: WF_RefineCluster_leaf_centroid_newmatrix.py [options]"
	parser = argparse.ArgumentParser(usage)
	parser.add_argument('-t', '--tree', dest="species_tree", required=True, help="Species tree. (Required)")
	parser.add_argument('-r', '--repo', dest="cobra_repo", required=True, help="Cobra repository. (Required)")
	parser.add_argument('-w', '--working', dest="working_dir", required=True, help="Working directory. (Required)")
	parser.add_argument('-f', '--flow', dest="flow_name", default="default_flow", help="Flow Name.")
	parser.add_argument('-a', '--alpha', type=float, dest="alpha", default=10.0, help="Synteny weight")
	parser.add_argument('-b', '--beta', type=float, dest="beta", default=10.0, help="Homology weight")
	parser.add_argument('-g', '--gamma', type=float, dest="gamma", default=10.0, help="Gain/Loss weight.")
	parser.add_argument('-G', '--gain', type=float, dest="gain", default=0.05, help="Duplication rate for Poisson distribution.")
	parser.add_argument('-L', '--loss', type=float, dest="loss", default=0.05, help="Loss rate for Poisson distribution.")
	parser.add_argument('-m', '--min_best_hit', type=float, dest="min_best_hit", default=0.5, help="Minimal % of match length for Blastp hits compared to best one.")
	parser.add_argument('-B', '--blast_eval', type=float, dest="blast_eval", default=#BLAST_EVAL_DEFAULT, help="Minimal e-value for Blastp hits.")
	parser.add_argument('-l', '--locus', dest="locus_file", default="", help="Locus file.")
	parser.add_argument('-n', '--num_cores', type=int, dest="num_cores", default=#NUM_CORES_DEFAULT, help="Locus file.")
	parser.add_argument('-F', '--min_syntenic_fraction', type=float, dest="minSynFrac", default=0.5, help="Minimal syntenic fraction.")
	parser.add_argument('-D', '--hamming', type=float, dest="hamming", default=0.6, help="Hamming distance for representative selection.")
	parser.add_argument('-s', '--synteny_window', type=int, dest="synteny_window", default=6000, help="Size of the syntenic window.")
	args = parser.parse_args()
	
	if not args.working_dir[-1] == "/":
		args.working_dir = args.working_dir + "/"
	genome_dir = args.working_dir + "genomes/"
	#set up genomes folder in synergy2 directory
	if not "genomes" in os.listdir(args.working_dir):
		os.system( "mkdir " + genome_dir)
	node_dir = args.working_dir + "nodes/"
	if not "nodes" in os.listdir(args.working_dir):
		os.system("mkdir " + node_dir)
	
	cobra_repo_path = "/".join(args.cobra_repo.split("/")[0:-1])+"/"
	
	#read COBRA repository and set up file system
	myRepo = COBRA_Repo_Handling.RepoParse(args.cobra_repo)
	myRepo.parseRepoFile(cobra_repo_path)
	if len(args.locus_file)>0:
		myRepo.readLocusTagFile(args.locus_file)
		myRepo.assignLocusTags()
		myRepo.writeLocusTagFile(genome_dir)
	elif "locus_tag_file.txt" in os.listdir(genome_dir):
		myRepo.readLocusTagFile(genome_dir+"locus_tag_file.txt")
		myRepo.assignLocusTags()
		myRepo.writeLocusTagFile(genome_dir)
	else:
		myRepo.assignLocusTags()
		myRepo.writeLocusTagFile(genome_dir)
	
	retVal = myRepo.makeGenomeDirectories(genome_dir, distribute, args.synteny_window)
	if retVal:
		#currently this will never happen, but if you are doing a whole bunch of genomes it may become necessary
		sys.exit("Launch this command file on the grid: "+ retVal)
		
		
	#read species tree
	logger.debug("init tree lib")
	myTree = TreeLib.Tree(args.species_tree, genome_dir+"locus_tag_file.txt")
	logger.info("reading genome to locus")
	myTree.readGenomeToLocusFile()
	logger.info("reading tree")
	myTree.readTree()
	logger.info("parsing tree")
	myRootEdge = myTree.parseTree()
	logger.info(myRootEdge)
	rn = myRootEdge.split(",")
	root_edge2 = (rn[0].split(":")[0],rn[1].split(":")[0])
	my_re = []
	for re in root_edge2:
		if re in myTree.genomeToLocus:
			my_re.append(myTree.genomeToLocus[re])
		else:
			my_re.append(re)
	root_edge = (my_re[0],my_re[1])
	
	logger.info("rooting")
	logger.info(root_edge)
	myTree.rootTree(root_edge)
	logger.info("initializing")
	myInitTree = WF_Initialization.Tree(myTree, args.flow_name, args.blast_eval, args.num_cores, args.alpha, args.beta, args.gamma, args.gain, args.loss, args.min_best_hit, args.synteny_window, args.homScale, args.synScale, args.minSynFrac, args.hamming)
	logger.info("dependicizing")
	myInitTree.calculateNodeDependencies(args.working_dir)
	myInitTree.writeLocusTagFile()

	logger.info('Finished')
		
if __name__ == "__main__":
	if len(sys.argv) == 1:
		sys.exit(usage())
	else:
		main(sys.argv[1:])