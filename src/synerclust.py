#!/usr/bin/env python

import sys
import os
import argparse
import logging
import subprocess
import shlex
import COBRA_Repo_Handling
import WF_Initialization
import TreeLib


def main():
	FORMAT = "%(asctime)-15s %(levelname)s %(module)s.%(name)s.%(funcName)s at %(lineno)d :\n\t%(message)s\n"
	logger = logging.getLogger()
	logging.basicConfig(filename='run_SynerClust.log', format=FORMAT, filemode='w', level=logging.DEBUG)
	# add a new Handler to print all INFO and above messages to stdout
	ch = logging.StreamHandler(sys.stdout)
	ch.setLevel(logging.INFO)
	logger.addHandler(ch)
	logger.info('Started')

	# usage = "usage: WF_RefineCluster_leaf_centroid_newmatrix.py [options]"
	parser = argparse.ArgumentParser()
	parser.add_argument('-t', '--tree', dest="species_tree", required=True, help="Species tree relating all of the genomes to be analyzed. (Required)")
	parser.add_argument('-r', '--repo', dest="cobra_repo", required=True, help="Complete path to data_catalog in the repository containing your genomic data. (Required)")
	parser.add_argument('-w', '--working', dest="working_dir", required=True, help="Complete path to the working directory for this analysis. (Required)")
	parser.add_argument('-m', '--min-best-hit', type=float, dest="min_best_hit", default=0.8, help="Minimal fraction of match length for Blastp hits compared to best one. (default = 0.8)")
	parser.add_argument('-B', '--blast_eval', type=float, dest="blast_eval", default=#BLAST_EVAL_DEFAULT, help="Minimal e-value for Blastp hits. (default = #BLAST_EVAL_DEFAULT)")
	parser.add_argument('-i', '--min-percent-identity', type=float, dest="min_percent_identity", default=0.5, help="Minimal fraction of identity filtering in the protein Blast+ results. (default = 0.5)")
	parser.add_argument('-c', '--min-match-coverage', type=float, dest="min_match_coverage", default=0.5, help="Minimal fraction of the query length that the protein Blast+ match needs to cover. (default = 0.5)")
	parser.add_argument('-l', '--locus', dest="locus_file", default="", help="A locus_tag_file.txt that corresponds to the data in this repository")
	parser.add_argument('-N', '--newick_tag', dest="coded_nwk_file", default="coded_tree.nwk", help="Output file for the newick tree using tag names and number of genomes as distances.")
	parser.add_argument('-n', '--num-threads', type=int, dest="num_cores", default=#NUM_CORES_DEFAULT, help="The (maximal) number of threads used for each node. (default = #NUM_CORES_DEFAULT)")
	parser.add_argument('-nmin', '--min-num-threads', type=int, dest="min_cores", default=None, help="If specified, only the root will use as many threads as provided in -n and other nodes will use less based on how far from the root they are. This is the minimal number of threads used for other nodes.")
	parser.add_argument('-F', '--min_syntenic_fraction', type=float, dest="minSynFrac", default=0.7, help="Minimum common syntenic fraction required for two genes from the same species to be considered paralogs, range [0.0,1.0], default=0.7")
	parser.add_argument('-D', '--dist', type=float, dest="dist", default=1.2, help="Maximum FastTree distance between a representative sequence and sequences being represented for representative selection. (default = 1.2)")
	parser.add_argument('-s', '--synteny_window', type=int, dest="synteny_window", default=6000, help="Distance in base pairs that will contribute to upstream and downstream to syntenic fraction. The total window size is [int]*2. (default = 6000")
	parser.add_argument('-p', '--project', type=str, dest="project", default=None, help="Project name for UGER queue system priority. (unused if not run in UGER mode)")
	parser.add_argument('-qs', '--queue-short', type=str, dest="shortq", default="short", help="Name of the queue used for jobs with a level lower to that defined in --queue-threshold. (unused if not run in UGER mode)(default = \"short\")")
	parser.add_argument('-ql', '--queue-long', type=str, dest="longq", default="long", help="Name of the queue used for jobs with a level higher or equal to that defined in --queue-threshold. (unused if not run in UGER mode)(default = \"long\")")
	parser.add_argument('-qt', '--queue-threshold', type=int, dest="thresq", default=25, help="Branching depth threshold (the XX in N_000000XX) to determine whether to use the short or long queue (unused if not run in UGER mode)(default = 25)")
	parser.add_argument('--no-synteny', dest="synteny", default=True, action='store_false', required=False, help="Disable use of synteny (required if information not available).")
	parser.add_argument('--run', default="none", type=str, dest="run", choices=["none", "single", "uger"], help="Specify if you want all computation to be run directly. Use \"single\" to run the local machine or \"uger\" to submit to the UGE grid.")
	parser.add_argument('--alignment', default="none", type=str, dest="alignment", choices=["none", "scc", "all"], help="Specify if you want cluster alignments using MUSCLE to be computed and written for the root node. Use \"all\" if you want all clusters to be aligned or \"scc\" if you only want Single Copy Core clusters to be aligned.")
	args = parser.parse_args()

	if not args.min_cores:
		args.min_cores = args.num_cores
	elif args.min_cores > args.num_cores:
		exit("Error: Minimal number of threads specified lesser than maximal number specified.")

	args.working_dir = os.path.abspath(args.working_dir) + "/"
	# if not args.working_dir[-1] == "/":
	# 	args.working_dir = args.working_dir + "/"
	genome_dir = args.working_dir + "genomes/"
	# set up genomes folder in synergy2 directory
	if "genomes" not in os.listdir(args.working_dir):
		os.system("mkdir " + genome_dir)
	node_dir = args.working_dir + "nodes/"
	if "nodes" not in os.listdir(args.working_dir):
		os.system("mkdir " + node_dir)
	if "logs" not in os.listdir(args.working_dir):
		os.system("mkdir " + args.working_dir + "logs/")

	# cobra_repo_path = "/".join(args.cobra_repo.split("/")[0:-1]) + "/"
	cobra_repo_path = args.cobra_repo[0:args.cobra_repo.rfind("/") + 1]
	args.cobra_repo = os.path.abspath(args.cobra_repo)
	cobra_repo_path = os.path.abspath(cobra_repo_path) + "/"

	# read COBRA repository and set up file system
	myRepo = COBRA_Repo_Handling.RepoParse(args.cobra_repo)
	myRepo.parseRepoFile(cobra_repo_path)
	if len(args.locus_file) > 0:
		myRepo.readLocusTagFile(args.locus_file)
	elif "locus_tag_file.txt" in os.listdir(genome_dir):
		myRepo.readLocusTagFile(genome_dir + "locus_tag_file.txt")

	myRepo.assignLocusTags()
	myRepo.writeLocusTagFile(genome_dir)

	retVal = myRepo.makeGenomeDirectories(genome_dir, args.synteny_window)

	# read species tree
	logger.debug("init tree lib")
	myTree = TreeLib.Tree(args.species_tree, genome_dir + "locus_tag_file.txt")
	logger.info("reading genome to locus")
	myTree.readGenomeToLocusFile()
	logger.info("reading tree")
	myTree.readTree()
	logger.info("parsing tree")
	myRootEdge = myTree.parseTree()
	logger.info(myRootEdge)
	rn = myRootEdge.split(",")
	root_edge2 = (rn[0].split(":")[0], rn[1].split(":")[0])
	my_re = []
	for re in root_edge2:
		if re in myTree.genomeToLocus:
			my_re.append(myTree.genomeToLocus[re])
		else:
			my_re.append(re)
	root_edge = (my_re[0], my_re[1])

	logger.info("rooting")
	logger.info(root_edge)
	myTree.rootTree(root_edge)
	logger.info("initializing")
	myInitTree = WF_Initialization.Tree(myTree, args.blast_eval, args.num_cores, args.min_best_hit, args.min_percent_identity, args.min_match_coverage, args.synteny_window, args.minSynFrac, args.dist, args.synteny)
	logger.info("dependicizing")

	myInitTree.calculateNodeDependencies(args.working_dir, args.shortq, args.longq, args.thresq, args.project, args.num_cores, args.min_cores)
	myInitTree.writeLocusTagFile()
	myInitTree.writeCodedNewick(genome_dir + args.coded_nwk_file)

	os.chdir(args.working_dir)
	cmd = ["#SYNERCLUST_PATHPostProcessingScript.sh", args.alignment]
	process = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	(output, error) = process.communicate()
	if error:
		exit(error)

	logger.info('Finished Initialization')
	if retVal:
		# currently this will never happen, but if you are doing a whole bunch of genomes it may become necessary
		logger.info("Please run annotation extraction before starting the actual computation: " + retVal)

	if args.run == "single":
		# run genome extraction
		if "needed_extractions.cmd.txt" in os.listdir(args.working_dir + "genomes/"):
			logger.info("Starting annotation extraction.\n")
			cmd = [args.working_dir + "genomes/needed_extractions.cmd.txt"]
			process = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			(output, error) = process.communicate()
			if process.returncode != 0:
				exit("Error while extracting genome annotations:\n" + error)
			logger.info("Finished annotation extraction.\n\n")
		# run jobs
		logger.info("Starting computing jobs.\n")
		cmd = [args.working_dir + "jobs.sh"]
		process = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		while process.poll() is None:
			l = process.stdout.readline()
			if l[:5] == "Node ":
				logger.info(l)
		logger.info(process.stdout.read())
		# (output, error) = process.communicate()
		if process.returncode != 0:
			exit("Error while computing jobs. Please check logs.\n")
		logger.info("Finished computing jobs.\n\n")
		# run root postprocessing
		logger.info("Starting postprocessing root.\n")
		cmd = [args.working_dir + "post_process_root.sh"]
		process = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		(output, error) = process.communicate()
		if process.returncode != 0:
			exit("Error while post processing root. Please check logs.\n")
		logger.info("Finished postprocessing root.\n\n")

	elif args.run == "uger":
		if "tmp" not in os.listdir(args.working_dir):
			os.system("mkdir " + args.working_dir + "tmp/")
		# run genome extraction
		if "needed_extractions.cmd.txt" in os.listdir(args.working_dir + "genomes/"):
			logger.info("Starting annotation extraction.\n")
			cmd = shlex.split("#SYNERCLUST_PATHuger_auto_submit_simple.py -f " + args.working_dir + "genomes/needed_extractions.cmd.txt -t 30 -tmp " + args.working_dir + "tmp/" + (" -p " + args.project if args.project else "") + " -err " + args.working_dir + "/logs" + " -log " + args.working_dir + "/logs" + " -hh 2 -mm 00 -ss 00")
			process = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			(output, error) = process.communicate()
			if process.returncode != 0:
				logger.debug(output)
				logger.error(error)
				exit("Error while extracting genome annotations on the grid:\n")
			logger.info("Finished annotation extraction.\n\n")
		# run jobs
		logger.info("Starting computing jobs.\n")
		cmd = shlex.split("#SYNERCLUST_PATHuger_auto_submit.py -f " + args.working_dir + "uger_jobs.sh -t 30 -l 900 -n " + str(args.num_cores))
		process = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		(output, error) = process.communicate()
		if process.returncode != 0:
			logger.debug(output)
			logger.error(error)
			exit("Error while computing jobs on the grid. Please check logs.\n" + error)
		logger.info("Finished computing jobs.\n\n")
		# run root postprocessing
		logger.info("Starting postprocessing root.\n")
		cmd = shlex.split("#SYNERCLUST_PATHuger_auto_submit_simple.py -f " + args.working_dir + "post_process_root.sh -t 30 -q long -tmp " + args.working_dir + "tmp/" + (" -p " + args.project if args.project else "") + " -hh 48 -mm 00 -ss 00")
		process = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		(output, error) = process.communicate()
		if process.returncode != 0:
			logger.debug(output)
			logger.error(error)
			exit("Error while post processing root on the grid. Please check logs.\n")
		logger.info("Finished postprocessing root.\n\n")


if __name__ == "__main__":
	main()
