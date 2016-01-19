#!/usr/bin/env python

import sys, os
def main(argv):
	scc_dir = argv[0]
	gblocks_extension = argv[1]
	
	seqs = {}
	
	files = os.listdir(scc_dir)
	gblocks_files = {}
	for f in files:
		genomesInFile = set([])
		if f.find(gblocks_extension)>-1:
			gblocks_files[f] = 0
			block = open(scc_dir+f,'r').readlines()
			curGenome = ""
			for b in block:
				b = b.rstrip()
				if b.find(">")>-1:
					if len(curGenome) > 2 and not curGenome == "THROW_AWAY":
						gblocks_files[f] = len(seqs[curGenome][f])
						#~ print >> sys.stderr, f, len(seqs[curGenome][f])
					curGenome = b.split(";")[1]
					if curGenome in genomesInFile:
						curGenome = "THROW_AWAY"
					genomesInFile.add(curGenome)
					if not curGenome in seqs:
						seqs[curGenome] = {}
				else:
					b = b.replace(" ","")
					if not f in seqs[curGenome]:
						seqs[curGenome][f] = ""
					seqs[curGenome][f] += b
					
	for s in seqs:
		header = ">"+s
		print header
		mySeq = ""
		for g in gblocks_files:
			if not g in seqs[s]:
				#~ print >> sys.stderr, "not in file", s, gblocks_files[g], g
				#~ print g, gblocks_files[g]
				seqs[s][g] = gblocks_files[g]*"-"
			mySeq += seqs[s][g]
		mySeq = mySeq.upper()
		print mySeq
	

if __name__ == "__main__":
	main(sys.argv[1:])