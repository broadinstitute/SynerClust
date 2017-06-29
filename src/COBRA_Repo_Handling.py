import sys
# import random
# import string
import os
import logging
import hashlib
import base64
# import traceback


class RepoParse:
	logger = logging.getLogger("RepoParse")

	def __init__(self, repo_file):
		self.repo_file = repo_file
		self.genomeToLocus = {}
		self.locusToGenome = {}
		self.nodeChildrenCount = {}
		self.genomes = []
		# self.locusTags = set([])
		RepoParse.logger.debug("RepoParse initialized")

	def parseRepoFile(self, repo_path):
		RepoParse.logger.debug("Parsing repo files")
		with open(self.repo_file) as f:
			# lines = f.readlines()
			# print lines
			curGenome = {}
			curGenome["transl_table"] = '1'
			# for l in lines:
			for l in f:
				if l[0] == "#":
					continue
				if l.find("//") > -1:
					# new genome
					if len(curGenome) > 0:
						if "Genome" in curGenome and "Annotation" in curGenome and "Sequence" in curGenome:  # verify that all 3 important fields have been filed
							if not os.path.isfile(curGenome["Sequence"]):
								RepoParse.logger.error("curGenome[\"Sequence\"] not found at : %s" % (curGenome["Sequence"]))
								sys.exit("No sequence file specified and not found at default location for %s" % (curGenome["Genome"]))
							RepoParse.logger.debug("Creating genome entry with %s %s %s %s" % (curGenome["Genome"], curGenome["Annotation"], curGenome["Sequence"], curGenome["Peptide"]))
							myG = Genome(curGenome["Genome"], curGenome["Annotation"], curGenome["Sequence"], curGenome["Peptide"], curGenome["transl_table"])
							self.genomes.append(myG)
						elif "Genome" in curGenome and "Peptide" in curGenome and curGenome["Peptide"] is not "null":  # might need to swap last 2 conditions position based on direction of evaluation
							RepoParse.logger.debug("Creating genome entry with %s %s" % (curGenome["Genome"], curGenome["Peptide"]))
							myG = Genome(curGenome["Genome"], "", "", curGenome["Peptide"], "")
							self.genomes.append(myG)
						curGenome = {}
						curGenome["transl_table"] = '1'
					continue
				dat = l.split()
				if len(dat) < 2:
					continue
				if dat[0].find("Genome") > -1:
					curGenome["Peptide"] = "null"
					curGenome["Genome"] = dat[1]  # need to use the hardnamed key in case "Genome" found in for example "Genomes"
				elif dat[0].find("Annotation") > -1:
					if not os.path.isabs(dat[1]):  # if not an absolute path
						# if os.path.isabs(dat[1]): # if not an absolute path
						# if not dat[1][0] == "/": # if not an absolute path
						# if not dat[1].find("/") > -1:
							# genome_path = repo_path + curGenome["Genome"] + "/"
							# genome_path = repo_path
							# myfiles = os.listdir(genome_path)
						seq = None
						if os.path.isfile(repo_path + dat[1]):  # path relative to repo, "./" and "../" can be used
							RepoParse.logger.warning("Are you sure \"%s\" is a gff3 annotation file?" % (dat[1]))
							seq = repo_path + curGenome["Genome"] + ".genome.fa"
							dat[1] = repo_path + dat[1]
						elif os.path.isfile(repo_path + dat[1] + ".annotation.gff3"):
							seq = repo_path + curGenome["Genome"] + ".genome.fa"
							dat[1] = repo_path + dat[1] + ".annotation.gff3"
						elif os.path.isfile(repo_path + curGenome["Genome"] + "/" + dat[1]):
							RepoParse.logger.warning("Are you sure \"%s\" is a gff3 annotation file?" % (dat[1]))
							seq = repo_path + curGenome["Genome"] + "/" + curGenome["Genome"] + ".genome.fa"
							dat[1] = repo_path + curGenome["Genome"] + "/" + dat[1]
						elif os.path.isfile(repo_path + curGenome["Genome"] + "/" + dat[1] + ".annotation.gff3"):
							seq = repo_path + curGenome["Genome"] + "/" + curGenome["Genome"] + ".genome.fa"
							dat[1] = repo_path + curGenome["Genome"] + "/" + dat[1] + ".annotation.gff3"
						else:
							RepoParse.logger.error("Specified annotation file not found for %s" % (curGenome["Genome"]))
							sys.exit("Specified annotation file not found for %s" % (curGenome["Genome"]))
						# curGenome["Sequence"] = genome_path + curGenome["Genome"] + ".genome.fa"
						# if not curGenome["Genome"] + ".genome.fa" in myfiles:
						# 	curGenome["Sequence"] = genome_path + curGenome["Genome"] +".genome"
					# else:
						# sys.exit("Specified annotation file not found for : %s" %(dat[1]))
					curGenome["Annotation"] = dat[1]
					if "Sequence" not in curGenome:  # in case the Sequence is provided before the Annotation
						curGenome["Sequence"] = seq
				elif dat[0].find("Sequence") > -1:
					if not os.path.isabs(dat[1]):  # if not an absolute path
						# if not dat[1][0] == "/": # if not an absolute path
						if os.path.isfile(repo_path + dat[1]):
							dat[1] = repo_path + dat[1]
						elif os.path.isfile(repo_path + dat[1] + ".fasta"):
							dat[1] = repo_path + dat[1] + ".fasta"
						elif os.path.isfile(repo_path + dat[1] + ".genome"):
							dat[1] = repo_path + dat[1] + ".genome"
						elif os.path.isfile(repo_path + dat[1] + ".genome.fa"):
							dat[1] = repo_path + dat[1] + ".genome.fa"
						elif os.path.isfile(repo_path + curGenome["Genome"] + "/" + dat[1]):
							dat[1] = repo_path + curGenome["Genome"] + "/" + dat[1]
						elif os.path.isfile(repo_path + curGenome["Genome"] + "/" + dat[1] + ".fasta"):
							dat[1] = repo_path + curGenome["Genome"] + "/" + dat[1] + ".fasta"
						elif os.path.isfile(repo_path + curGenome["Genome"] + "/" + dat[1] + ".genome"):
							dat[1] = repo_path + curGenome["Genome"] + "/" + dat[1] + ".genome"
						elif os.path.isfile(repo_path + curGenome["Genome"] + "/" + dat[1] + ".genome.fa"):
							dat[1] = repo_path + curGenome["Genome"] + "/" + dat[1] + ".genome.fa"
					else:
						if not os.path.isfile(dat[1]):
							sys.exit("Specified sequence file not found : %s" % (dat[1]))
					curGenome["Sequence"] = dat[1]
				elif dat[0].find("Peptide") > -1:
					if not os.path.isabs(dat[1]):
						if os.path.isfile(repo_path + dat[1]):
							dat[1] = repo_path + dat[1]
						elif os.path.isfile(repo_path + curGenome["Genome"] + "/" + dat[1]):
							dat[1] = repo_path + curGenome["Genome"] + "/" + dat[1]
						curGenome["Peptide"] = dat[1]
					else:
						sys.exit("Specified peptide file not found : %s" % (dat[1]))
				elif dat[0].find("transl_table") > -1:
					curGenome["transl_table"] = dat[1]
			if len(curGenome) > 1:
				myG = Genome(curGenome["Genome"], curGenome["Annotation"], curGenome["Sequence"], curGenome["Peptide"], curGenome["transl_table"])
				self.genomes.append(myG)

	def assignLocusTags(self):
		for g in self.genomes:
			if g.genome in self.genomeToLocus:
				g.locus = self.genomeToLocus[g.genome]
			else:
				g.locus = self.assignGenomeLocus(g.genome)
			self.genomeToLocus[g.genome] = g.locus
			self.locusToGenome[g.locus] = g.genome
			self.nodeChildrenCount[g.locus] = 1

	def readLocusTagFile(self, tag_file):
		tags = open(tag_file, 'r').readlines()
		for t in tags:
			t.rstrip()
			line = t.split()
			self.genomeToLocus[line[0]] = line[1]
			self.locusToGenome[line[1]] = line[0]
			self.nodeChildrenCount[line[1]] = int(line[2])
			# self.locusTags.add(line[1])

	def writeLocusTagFile(self, genome_dir):
		tag_out = open(genome_dir + "locus_tag_file.txt", 'w')
		for g in self.genomeToLocus:
			line = "\t".join([g, self.genomeToLocus[g], str(self.nodeChildrenCount[self.genomeToLocus[g]])]) + "\n"
			tag_out.write(line)
		tag_out.close()
		RepoParse.logger.info("Wrote locus tags to locus_tag_file.txt")

	def assignGenomeLocus(self, genome):
		# L for leaf
		# 000000 because it is a leaf/basic level
		# hash of the genome name, the [:-2] is to discard the "==" from the encoding
		code = "L_0000000_" + base64.urlsafe_b64encode(hashlib.md5(genome).digest())[:-2]
		return code

	def makeGenomeDirectories(self, genome_dir, synteny_window):
		command_lists = []
		for g in self.genomes:
			ret = g.setupDirectory(genome_dir, synteny_window)
			if ret:
				command_lists.append(ret)
		if len(command_lists) > 0:
			cat_cmds = open(genome_dir + "needed_extractions.cmd.txt", 'w')
			cat_cmds.write("#!/usr/bin/env bash\n\n")
			for cl in command_lists:
				cat_cmds.write(cl)
			cat_cmds.close()
			os.chmod(genome_dir + "needed_extractions.cmd.txt", 0775)
			return genome_dir + "needed_extractions.cmd.txt"
		elif "needed_extractions.cmd.txt" in os.listdir(genome_dir):
			cmd = "rm " + genome_dir + "needed_extractions.cmd.txt"
			os.system(cmd)
		return False


class Genome:

	logger = logging.getLogger(__name__)

	def __init__(self, genome, annotation, sequence, peptide, transl_table):
		self.genome = genome
		self.locus = ""
		self.annotation = annotation
		self.sequence = sequence
		self.peptide = peptide
		self.directory = ""
		self.transl_table = transl_table
		Genome.logger.debug("Genome Initialized")

	def setupDirectory(self, genomeDir, synteny_window):
		# make a directory for this genome in the genomeDir
		# myCommands = []
		myDir = genomeDir + self.genome + "/"
		self.directory = myDir
		if self.genome not in os.listdir(genomeDir):
			os.mkdir(myDir)
		myDirFiles = os.listdir(myDir)
		annot = " 0 "
		if "annotation.txt" in myDirFiles:
			pass
			# print "Annotation information already written for "+self.genome+", "+self.locus
		else:
			annot = " 1 "
		nodeDir = genomeDir.replace("genomes", "nodes")
		myNodeDir = nodeDir + self.locus + "/"
		if self.locus not in os.listdir(nodeDir):
			os.mkdir(myNodeDir)
		myNodeFiles = os.listdir(myNodeDir)
		pickles = " 0 "

		myDir = os.path.abspath(myDir) + "/"
		if "PICKLES_COMPLETE" in myNodeFiles:
			pass
			# print self.locus+" pickles are complete!"
		else:
			pickles = " 1 "
		if pickles == " 0 " and annot == " 0 ":
			return False
		# if distribute == 1:
		if self.annotation is not "" and self.sequence is not "":
			cmd = "#SYNERCLUST_PATHFormatAnnotation_external.py -gff " + self.annotation + " -seq " + self.sequence + " -name " + self.genome + " --pep " + self.peptide + " -locus " + self.locus + " -out " + myDir + "annotation.txt -synteny " + str(synteny_window) + " --annot " + annot + " --pickle " + pickles + " --transl_table " + self.transl_table + "\n"
		elif self.peptide is not "":
			cmd = "#SYNERCLUST_PATHFormatAnnotation_external.py -name " + self.genome + " --pep " + self.peptide + " -locus " + self.locus + " -out " + myDir + "annotation.txt --annot " + annot + " --pickle " + pickles + "\n"
		return cmd
