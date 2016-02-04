import sys, random, string, os, logging, hashlib, base64
import traceback


class RepoParse:
	logger = logging.getLogger("RepoParse")
	
	def __init__(self, repo_file):
		self.repo_file = repo_file
		self.genomeToLocus = {}
		self.locusToGenome = {}
		self.genomes = []
		self.locusTags = set([])
		RepoParse.logger.debug("RepoParse initialized")
		
	def parseRepoFile(self,repo_path):
		RepoParse.logger.debug("Parsing repo files")
		with open(self.repo_file) as f:
			lines = f.readlines()
			# print lines
		# lines = open(self.repo_file,'r').readlines()
			curGenome = {}
			for l in lines:
				if l.find("#") > -1:
					continue
				if l.find("//") > -1:
					#new genome
					if len(curGenome) > 0:
						if "Genome" in curGenome and "Annotation" in curGenome and "Sequence" in curGenome: # verify that all 3 important fields have been filed
							RepoParse.logger.debug("Creating genome entry with %s %s %s %s" %(curGenome["Genome"], curGenome["Annotation"], curGenome["Sequence"], curGenome["Peptide"]))
							myG = Genome(curGenome["Genome"], curGenome["Annotation"], curGenome["Sequence"], curGenome["Peptide"])
							self.genomes.append(myG)
						curGenome = {}
					continue
				dat = l.split()
				if len(dat) < 2:
					continue
				if dat[0].find("Genome") > -1:
					curGenome["Peptide"] = "null"
					curGenome["Genome"] = dat[1] # need to use the hardnamed key in case "Genome" found in for example "Genomes"
				elif dat[0].find("Annotation") > -1:
					if not dat[1][0] == "/": # if not an absolute path
					# if not dat[1].find("/") > -1:
						# genome_path = repo_path + curGenome["Genome"] + "/"
						# genome_path = repo_path
						# myfiles = os.listdir(genome_path)
						if os.path.isfile(repo_path + dat[1]): # path relative to repo, "./" and "../" can be used
							RepoParse.logger.warning("Are you sure \"%s\" is a gff3 annotation file?" %(dat[1]))
							dat[1] = repo_path + dat[1]
						elif os.path.isfile(repo_path + dat[1] + ".annotation.gff3"):
							dat[1] = repo_path + dat[1] + ".annotation.gff3"
						elif os.path.isfile(repo_path + curGenome["Genome"] + dat[1]):
							RepoParse.logger.warning("Are you sure \"%s\" is a gff3 annotation file?" %(dat[1]))
							dat[1] = repo_path + curGenome["Genome"] + dat[1]
						elif os.path.isfile(repo_path + curGenome["Genome"] + dat[1] + ".annotation.gff3"):
							dat[1] = repo_path + curGenome["Genome"] + dat[1]  + "annotation.gff3"
						else:
							RepoParse.logger.error("Specified annotation file not found for %s" %(curGenome["Genome"]))
							sys.exit("Specified annotation file not found for %s" %(curGenome["Genome"]))
						# curGenome["Sequence"] = genome_path + curGenome["Genome"] + ".genome.fa"
						# if not curGenome["Genome"] + ".genome.fa" in myfiles:
						# 	curGenome["Sequence"] = genome_path + curGenome["Genome"] +".genome"
					curGenome["Annotation"] = dat[1]
				elif dat[0].find("Sequence") > -1:
					if not ".genome" in dat[1] or not ".fasta" in dat[1]: # .genome.fa is included by just .genome
						RepoParse.logger.warning("Are you sure \"%s\" is a fasta file?" %(dat[1]))
					# if not dat[1].find("/") > -1:
					if not dat[1][0] == "/":
						dat[1] = repo_path + dat[1]
					curGenome["Sequence"] = dat[1]
				elif dat[0].find("Peptide") > -1:
					curGenome["Peptide"] = dat[1]
			if len(curGenome) > 0:
				myG = Genome(curGenome["Genome"], curGenome["Annotation"], curGenome["Sequence"], curGenome["Peptide"])
				self.genomes.append(myG)
		
	def assignLocusTags(self):
		for g in self.genomes:
			if g.genome in self.genomeToLocus:
				g.locus = self.genomeToLocus[g.genome]
			else:
				g.locus = self.assignGenomeLocus(g.genome)
			self.genomeToLocus[g.genome] = g.locus
			self.locusToGenome[g.locus] = g.genome
		
	def readLocusTagFile(self,tag_file):
		tags = open(tag_file,'r').readlines()
		for t in tags:
			t.rstrip()
			line = t.split()
			self.genomeToLocus[line[0]] = line[1]
			self.locusToGenome[line[1]] = line[0]
			self.locusTags.add(line[1])

	def writeLocusTagFile(self,genome_dir):
		tag_out = open(genome_dir+"locus_tag_file.txt", 'w')
		for g in self.genomeToLocus:
			line = "\t".join([g, self.genomeToLocus[g]])+"\n"
			tag_out.write(line)
		tag_out.close()
		RepoParse.logger.info("Wrote locus tags to locus_tag_file.txt")
				
	def assignGenomeLocus(self, genome):
		RepoParse.logger.debug("".join(traceback.format_stack()))
		# L for leaf
		# 000000 because it is a leaf/basic level
		# hash of the genome name, the [:-2] is to discard the "==" from the encoding
		code = "L_0000000_" + base64.urlsafe_b64encode(hashlib.md5(genome).digest())[:-2]

		# code = ''.join(random.choice(string.ascii_uppercase) for x in range(3))
		# while code in self.locusTags:
		# 	code = ''.join(random.choice(string.ascii_uppercase) for x in range(3))
		# self.locusTags.add(code)
		return code
			
	def makeGenomeDirectories(self, genome_dir, distribute,synteny_window):
		command_lists = []
		for g in self.genomes:
			ret = g.setupDirectory(genome_dir,distribute,synteny_window)
			if ret:
				command_lists.append(ret)
		if len(command_lists) >0:
			cat_cmds = open(genome_dir+"needed_extractions.cmd.txt",'w')
			for cl in command_lists:
				cat_cmds.write(cl)
			cat_cmds.close()
			return genome_dir+"needed_extractions.cmd.txt"
		elif "needed_extractions.cmd.txt" in os.listdir(genome_dir):
			cmd = "rm "+genome_dir+"needed_extractions.cmd.txt"
			os.system(cmd)
		return False
		
class Genome:
	logger = logging.getLogger(__name__)

	def __init__(self, genome, annotation, sequence,peptide):
		self.genome = genome
		self.locus = ""
		self.annotation = annotation
		self.sequence = sequence
		self.peptide = peptide
		self.directory = ""
		Genome.logger.debug("Genome Initialized")
		
	def setupDirectory(self, genomeDir,distribute,synteny_window):
		#make a directory for this genome in the genomeDir
		myCommands = []
		myDir = genomeDir+self.genome+"/"
		self.directory = myDir
		if not self.genome in os.listdir(genomeDir):
			os.mkdir(myDir)
		myDirFiles = os.listdir(myDir)
		annot = " 0 "
		if "annotation.txt" in myDirFiles:
			pass
			#~ print "Annotation information already written for "+self.genome+", "+self.locus
		else:
			annot = " 1 " 
		nodeDir = genomeDir.replace("genomes","nodes")
		myNodeDir = nodeDir+self.locus+"/"
		if not self.locus in os.listdir(nodeDir):
			os.mkdir(myNodeDir)
		myNodeFiles = os.listdir(myNodeDir)
		pickles = " 0 "
		if "PICKLES_COMPLETE" in myNodeFiles:
			pass
			#~ print self.locus+" pickles are complete!"
		else:
			pickles = " 1 "
		if pickles == " 0 " and annot == " 0 ":
			return False
		if distribute == 1:
			cmd = "#SYNERGY2_PATHFormatAnnotation_external.py "+self.annotation+" "+self.sequence+" "+self.genome+" "+self.peptide+" "+self.locus+" "+myDir+"annotation.txt "+str(synteny_window)+annot+pickles+"\n"
			return cmd
		else:
			cmd = "#SYNERGY2_PATHFormatAnnotation_external.py "+self.annotation+" "+self.sequence+" "+self.genome+" "+self.peptide+" "+self.locus+" "+myDir+"annotation.txt "+str(synteny_window)+annot+pickles+"\n"
			os.system(cmd)
			#~ cmd2 = "cp "+self.sequence+" "+myDir+"."
			#~ cmd3 ="mv "+myDir+self.sequence+" "+myDir+self.genome+".genome"
			
			Genome.logger.info("Wrote annotation information for " + self.genome)
		
		return False