#!/usr/bin/env python

import sys, os, getopt
import COBRA_Repo_Handling, WF_Initialization, TreeLib
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
	cobra_repo = ""
	cobra_repo_path = ""
	working_dir = ""
	species_tree = ""
	flow_name = ""
	genome_dir = ""
	blast_eval = #BLAST_EVAL_DEFAULT
	num_cores = #NUM_CORES_DEFAULT
	cmds_per_job = 500
	synteny_window = 6000
	alpha = 10.0
	gamma = 10.0
	gain = 0.05
	loss = 0.05
	min_best_hit = 0.5 
	homScale = 1.0
	synScale = 1.0
	numHits = 5
	minSynFrac = 0.5
	distribute = 1
	hamming = 0.6
	locus_file = ""	
	try:
		opts, args = getopt.getopt(argv, "r:w:t:f:g:G:L:l:m:b:c:n:s:H:S:l:T:F:D:h",["repo=","working=","tree=","flow=","gamma=","gain=","loss=","min_best_hit=","blast_eval=","cmds_per_job=","num_cores=","synteny_window=","homology_scale=","synteny_scale=","locus=","top_hits=","min_syntenic_fraction=","hamming=","help"])

	except getopt.GetoptError:
		usage()
		sys.exit(2)
	for opt, arg in opts:
		if opt in ("-t","--tree"):
			species_tree = arg
		elif opt in ("-r","--repo"):
			cobra_repo = arg	
			cobra_repo_path = "/".join(cobra_repo.split("/")[0:-1])+"/"
		elif opt in ("-w","--working"):
			working_dir = arg	
			if not working_dir[-1] == "/":
				working_dir = working_dir+"/"
			genome_dir = working_dir + "genomes/"
			#set up genomes folder in synergy2 directory
			if not "genomes" in os.listdir(working_dir):
				os.system( "mkdir "+genome_dir)
			node_dir = working_dir + "nodes/"
			if not "nodes" in os.listdir(working_dir):
				os.system("mkdir "+node_dir)
		elif opt in ("-f","--flow"):
			flow_name = arg	
		elif opt in ("-a","--alpha"):
			alpha = float(arg)
		elif opt in ("-g","--gamma"):
			gamma = float(arg)
		elif opt in ("-G","--gain"):
			gain = float(arg)
		elif opt in ("-L","--loss"):
			loss = float(arg)
		elif opt in ("-m","--min_best_hit"):
			min_best_hit = float(arg)
		elif opt in ("-b","--blast_eval"):
			blast_eval = float(arg)
		elif opt in ("-l","--locus"):
			locus_file = arg
		elif opt in ("-c","--cmds_per_job"):
			cmds_per_job = int(arg)
		elif opt in ("-n","--num_cores"):
			num_cores = int(arg)	
		elif opt in ("-F","--min_syntenic_fraction"):
			minSynFrac = float(arg)
		elif opt in ("-D","--hamming"):
			hamming = float(arg)
		elif opt in ("-T","--top_hits"):
			numHits = int(arg)	
		elif opt in ("-s","--synteny_window"):
			synteny_window = int(arg)
		elif opt in ("-H","--homology_scale"):
			myArg = float(arg)
			if not (myArg>0.0):
				sys.exit("-H or --homology_scale must be greater than 0.0 and less than or equal to 1.0.")
			else:
				homScale = myArg
		elif opt in ("-S","--synteny_scale"):
			myArg = float(arg)
			if not (myArg>0.0):
				sys.exit("-S or --synteny_scale must be greater than 0.0 and less than or equal to 1.0.")
			else:
				synScale = myArg
		elif opt in ("-h","--help"):
			sys.exit(usage())
		
	#read COBRA repository and set up file system
	myRepo = COBRA_Repo_Handling.RepoParse(cobra_repo)
	myRepo.parseRepoFile(cobra_repo_path)
	if len(locus_file)>0:
		myRepo.readLocusTagFile(locus_file)
		myRepo.assignLocusTags()
		myRepo.writeLocusTagFile(genome_dir)
	elif "locus_tag_file.txt" in os.listdir(genome_dir):
		myRepo.readLocusTagFile(genome_dir+"locus_tag_file.txt")
		myRepo.assignLocusTags()
		myRepo.writeLocusTagFile(genome_dir)
	else:
		myRepo.assignLocusTags()
		myRepo.writeLocusTagFile(genome_dir)
	
	retVal = myRepo.makeGenomeDirectories(genome_dir,distribute,synteny_window)
	if retVal:
		#currently this will never happen, but if you are doing a whole bunch of genomes it may become necessary
		sys.exit("Launch this command file on the grid: "+ retVal)
		
		
	#read species tree
	print "init tree lib"
	myTree = TreeLib.Tree(species_tree, genome_dir+"locus_tag_file.txt")
	print "reading genome to locus"
	myTree.readGenomeToLocusFile()
	print "reading tree"
	myTree.readTree()
	print "parsing tree"
	myRootEdge = myTree.parseTree()
	print myRootEdge
	rn = myRootEdge.split(",")
	root_edge2 = (rn[0].split(":")[0],rn[1].split(":")[0])
	my_re = []
	for re in root_edge2:
		if re in myTree.genomeToLocus:
			my_re.append(myTree.genomeToLocus[re])
		else:
			my_re.append(re)
	root_edge = (my_re[0],my_re[1])
	
	print "rooting"
	print root_edge
	myTree.rootTree(root_edge)
	#~ myTree.rootByMidpoint()
	print "initializing"
	myInitTree = WF_Initialization.Tree(myTree, flow_name,blast_eval, num_cores, alpha, gamma, gain, loss, min_best_hit, cmds_per_job,synteny_window,homScale,synScale,numHits, minSynFrac,hamming)
	print "dependicizing"
	myInitTree.calculateNodeDependencies(working_dir)
	myInitTree.writeLocusTagFile()
		
if __name__ == "__main__":
	if len(sys.argv) == 1:
		sys.exit(usage())
	else:
		main(sys.argv[1:])