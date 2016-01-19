#!/usr/bin/env python

import sys

class RepoParse:
	def __init__(self, repo_file, repo_dir, outfile):
		self.repo_file = repo_file
		self.repo_dir = repo_dir
		self.outfile = outfile
		self.annotation_files = []
		
	def parseRepoFile(self):
		lines = open(self.repo_file,'r').readlines()
		curG = ""
		for l in lines:
			if l.find("#") > -1:
				continue
			if l.find("//") > -1:
				continue
			l = l.rstrip()
			dat = l.split()
			if len(dat)<2:
				continue
			elif dat[0].find("Genome") > -1:
				curG = dat[1]
			elif dat[0].find("Annotation") > -1:
				annot = self.repo_dir+curG+"/"+dat[1]+".annotation.pep"
				self.annotation_files.append((curG,annot))
				
	def getPeptideSequences(self):
		out = open(self.outfile,'w')
		for afile in self.annotation_files:
			assembly = afile[0]
			annot_file = afile[1]
			af = open(afile[1], 'r').readlines()
			curGene = ""
			seq = ""
			for a in af:
				a= a.rstrip()
				if a.find(">")>-1:
					if len(curGene)>0:
						if seq.count("*") < 5:
							if seq[-1]=="*":
								tseq = seq[:-1]
								seq = tseq
							length = str(len(seq))
							header = ";".join([curGene,assembly,length])
							out.write(header+"\n")
							out.write(seq+"\n")
						seq = ""
						curGene = ""
					line = a.split()
					curGene = line[0]
				else:
					seq+=a
			if len(curGene)>0:
				if seq.count("*") < 5:
					if seq[-1]=="*":
						tseq = seq[:-1]
						seq = tseq
					length = str(len(seq))
					header = ";".join([curGene,assembly,length])
					out.write(header+"\n")
					out.write(seq+"\n")
		out.close()
		
def usage():
	sys.exit("do it right with a repo, dir, and outfile, ok sport?")
		
def main(argv):
	repo = argv[0]
	repo_dir = argv[1]
	outfile = argv[2]
	
	myParser = RepoParse(repo,repo_dir,outfile)
	myParser.parseRepoFile()
	myParser.getPeptideSequences()
	
if __name__ == "__main__":
	if len(sys.argv) == 1:
		sys.exit(usage())
	else:
		main(sys.argv[1:])
	
	