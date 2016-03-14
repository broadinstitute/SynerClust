import operator
from string import maketrans
class Genome:
	def __init__(self, file):
		self.file = file
		self.totalLength = 0
		self.scaffolds = {}
		self.stops = ['TAA','TAG','TGA']
		self.codons = {'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L','ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V','TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P','ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A','TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E','TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGA':'R','CGG':'R','CGC':'R','AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGG':'G','GGA':'G','GGC':'G'}
		
	def readSequence(self):
		genome_in = open(self.file, 'r').readlines()
		seqID = ""
		for l in genome_in:
			l = l.rstrip()
			if l.find('>') > -1:
				line = l.split()
				if len(seqID) > 0:
					self.totalLength += self.scaffolds[seqID].getSequenceLength()
				seqID = line[0][1:]
				scaffold = Scaffold(seqID, self.totalLength)
				self.scaffolds[seqID] = scaffold
			else:
				self.scaffolds[seqID].addSequence(l)
		if len(seqID) > 0:
			self.totalLength+= self.scaffolds[seqID].getSequenceLength()
				
	def readAndFormatSequence(self):
		genome_in = open(self.file, 'r').readlines()
		seqID = ""
		for l in genome_in:
			if l.find('>') > -1:
				l = l.rstrip()
				line = l.split()
				if len(seqID) > 0:
					self.totalLength += self.scaffolds[seqID].getSequenceLength()
				seqID = '_'.join(line)
				#~ print seqID
				scaffold = Scaffold(seqID, self.totalLength)
				self.scaffolds[seqID] = scaffold
			else:
				self.scaffolds[seqID].addRawSequence(l)
				l=l.rstrip()
				self.scaffolds[seqID].addSequence(l)
				
	def calcGenomeStats(self,scaffGeneCount,geneCount):
		lengths = []
		for s in self.scaffolds:
			gc = 0
			if s in scaffGeneCount:
				gc = scaffGeneCount[s]
			tup = (s,(self.scaffolds[s]).getSequenceLength(),gc)
			lengths.append(tup)
		n50_min = self.totalLength*.5
		n90_min = self.totalLength*.9
		g50_min = geneCount*.5
		g90_min = geneCount*.9
		#~ print n50_min, n90_min
		n50 = 0
		n90 = 0
		g50 = 0
		g90 = 0
		tlen = 0
		gcount = 0
		lengths = sorted(lengths, key=lambda tup: tup[1], reverse=True)
		for l in lengths:
			tlen += l[1]
			gcount += l[2]
			if tlen >= n50_min and n50 == 0:
				n50 = l[1]
			if tlen >= n90_min and n90 == 0:
				n90 = l[1]
			if gcount >= g90_min and g90 ==0:
				g90 = l[1]
			if gcount >= g50_min and g50 ==0:
				g50 = l[1]
		print "n50", n50
		print "n90", n90
		retlist = [str(self.totalLength),str(n50),str(n90),str(g50),str(g90)]
		return retlist
			
				
	def writeSequenceToFile(self, outfile):
		outFile = open(outfile, 'w')
		for s in self.scaffolds:
			outFile.write(s)
			outFile.write("\n")
			outFile.write(self.scaffolds[s].raw_sequence)
		outFile.close()
		
	def writeGapsToFile(self, gap_file):
		#gaps are based on a first nucleotide position of 1, not 0
		#format: [scaffold] [type=seq_end, gap_init, gap_term] [position of 1st or last N]
		outFile = open(gap_file, 'w')
		for s in self.scaffolds:
			for f in self.scaffolds[s].seqFeats:
				if f.isGap() or f.isSequenceEnd():
					out = "\t".join([str(s), f.type, str(f.position)]) + "\n"
					outFile.write(out)
		outFile.close()
				
	def featurizeScaffolds(self):
		for s in self.scaffolds:
			self.scaffolds[s].globalStop = self.scaffolds[s].globalStart + self.scaffolds[s].getSequenceLength() -1
			self.scaffolds[s].identifyGaps()
			mySeqStart = SequenceFeature(s,'sequence_end',0,-1,0)
			mySeqStop = SequenceFeature(s,'sequence_end',len(self.scaffolds[s].sequence) +1, -1, len(self.scaffolds[s].sequence) +1)
			self.scaffolds[s].seqFeats.append(mySeqStart)
			self.scaffolds[s].seqFeats.append(mySeqStop)
			self.scaffolds[s].seqFeats.sort(key=operator.attrgetter('position'))
	
	def coordsToGFF3(self, seqID, start, stop, strand):
		gene_structs = [] #list of tuples: (seqID, local_start, local_stop, strand)
		scaffID = seqID
		#check for gap overlaps
		subFeats = []
		gap_count = 0
		gap_encapsulations = []
		first_gap = None
		last_gap_term = None
		for f in self.scaffolds[scaffID].seqFeats:
			if f.position >= start and f.position <= stop:
				subFeats.append(f)
				if f.isGap():
					if f.isGapTerm():
						last_gap_term = f
					if first_gap == None:
						first_gap = f
					elif f.isGapInit():
						encaps = (last_gap_term, f)
						gap_encapsulations.append(encaps)
					gap_count += 1
			elif f.position > stop:
				break

		if gap_count == 0:
			myCDS = self.getCDS(scaffID, start, stop, strand)
			while not ((stop - start + 1) %3 == 0):
				if myCDS[-3:] in self.stops:
					if strand.find('+')>-1:
						start +=1
					else:
						stop -=1
				else:
					if strand.find('+')>-1:
						stop -=1
					else:
						start +=1
				myCDS = self.getCDS(scaffID, start, stop, strand)	
			myProtein = self.translateCDS(myCDS)
			myStruct = (scaffID, start, stop, strand, myCDS, myProtein)
			#print myStruct
			gene_structs.append(myStruct)

		elif gap_count > 0:
			#for f in subFeats:
				#print f.token
			gap_start = 0 #start in a gap?
			gap_stop = 0 #stop in a gap?
			if first_gap.isGapTerm():
				gap_start = 1
				if gap_count > 1 and gap_count % 2 == 0:
					gap_stop = 1
			else:
				#first_gap is gap init
				if gap_count % 2 == 1:
					gap_stop = 1
			if not gap_start:
				#print 'not gap start'
				frameNeeded = start%3 -1
				if frameNeeded == -1:
					frameNeeded = 2
				for f in subFeats:
					if f.frame == frameNeeded:
						myCDS = self.getCDS(scaffID, start, f.position, strand)
						myProtein = self.translateCDS(myCDS)
						myStruct = (scaffID, start, f.position, strand, myCDS, myProtein)
						#print myStruct
						gene_structs.append(myStruct)
						break
			if not gap_stop:
				#print 'not gap stop'
				frameNeeded = stop%3 +1
				if frameNeeded == 3:
					frameNeeded = 0
				i = -1
				while i > -4:
					if subFeats[i].frame == frameNeeded:
						myCDS = self.getCDS(scaffID, subFeats[i].position, stop, strand)
						myProtein = self.translateCDS(myCDS)
						myStruct = (scaffID, subFeats[i].position, stop, strand, myCDS, myProtein)
						#print myStruct
						gene_structs.append(myStruct)
					i-=1
			if len(gap_encapsulations) > 0:
				#print 'gap encapsulations'
				for g in gap_encapsulations:
					term_i = subFeats.index(g[0])
					init_i = subFeats.index(g[1])
					feat_init = subFeats[term_i +1]
					i = term_i + 2
					frameNeeded = feat_init.frame -1
					if frameNeeded == -1:
						frameNeeded = 2
					while i < init_i:
						if subFeats[i].frame == frameNeeded:
							myCDS = self.getCDS(scaffID, feat_init.position, subFeats[i].position, strand)
							myProtein = self.translateCDS(myCDS)
							myStruct = (scaffID, feat_init.position, subFeats[i].position, strand, myCDS, myProtein)
							#print myStruct
							gene_structs.append(myStruct)
						i+=1
		#~ if len(gene_structs) > 1:
			#~ print start, stop, strand
			#~ for g in gene_structs:
				#~ print g
		return gene_structs
		
	def getCDS(self, seqID, start, stop, strand):
		comp_start = start-1
		comp_stop = stop
		cds = self.scaffolds[seqID].sequence[comp_start:comp_stop]
		cds = cds.upper()
		if strand.find('+')>-1:
			return cds
		else:
			intab = "ACTG"
			outtab = "TGAC"
			transtab = maketrans(intab, outtab)
			cds = cds.translate(transtab)
			cds = cds[::-1]
			return cds
			
	def translateCDS(self, CDS):
		protein = ""
		index = 0
		while index < len(CDS):
			#do some stuff
			codon = CDS[index:index+3]
			amino = ''
			if codon.find('N')> -1:
				amino = 'X'
			elif not codon in self.codons:
				amino = 'X'
			else:
				amino = self.codons[codon]
			protein = protein + amino
			index = index+3
		return protein
		
	def translateCDSFile(self):
		for s in self.scaffolds:
			print ">"+s
			cds = self.scaffolds[s].sequence
			cds = cds.upper()
			print self.translateCDS(cds)

class Scaffold:
	def __init__(self, ID, globalStart):
		self.ID = ID
		self.sequence = ""
		self.raw_sequence = ""
		#globalStart and globalStop begin at 0, sequence features begin at 1
		self.globalStart = int(globalStart)
		self.globalStop = 0
		self.seqFeats = []

	def addSequence(self, sequence):
		self.sequence = self.sequence + sequence
		
	def addRawSequence(self, sequence):
		self.raw_sequence = self.raw_sequence + sequence
	
	def getSequenceLength(self):
		return len(self.sequence)
		
	def identifyGaps(self):
		index = self.sequence.find('NNN')
		while index > -1:
			myGapInit = SequenceFeature(self.ID,'gap_init',index+1, -1,index+self.globalStart+2)
			self.seqFeats.append(myGapInit)
			self.makeUpstreamFeatures(myGapInit)
			
			next_nucs = (self.sequence.find('A',index),self.sequence.find('C',index),self.sequence.find('G',index),self.sequence.find('T',index))
			next_nuc = min(next_nucs)
			
			myGapTerm = SequenceFeature(self.ID,'gap_term',next_nuc,-1,next_nuc+self.globalStart+1)
			self.seqFeats.append(myGapTerm)
			self.makeDownstreamFeatures(myGapTerm)
			index = self.sequence.find('NNN',next_nuc)
			
	def makeUpstreamFeatures(self, gapInitFeat):
		pos = gapInitFeat.position
		i = 3
		while i>0:
			newFeatPos = pos - i
			newFeatFrame = newFeatPos % 3
			upFeat = SequenceFeature(self.ID,'upstream',newFeatPos, newFeatFrame,newFeatPos+self.globalStart+1)
			self.seqFeats.append(upFeat)
			i = i-1
			
	def makeDownstreamFeatures(self, gapTermFeat):
		pos = gapTermFeat.position
		i = 3
		while i>0:
			newFeatPos = pos + i
			newFeatFrame = newFeatPos % 3
			downFeat = SequenceFeature(self.ID,'downstream',newFeatPos, newFeatFrame, newFeatPos+self.globalStart+1)
			self.seqFeats.append(downFeat)
			i =i-1
			
			
class SequenceFeature:
	def __init__(self,seqID,type,position,frame,globalPosition):
		self.seqID = seqID
		self.type = type #start,stop,gap_start,gap_stop,sequence_end
		self.position = int(position)
		self.globalPosition = globalPosition
		self.frame = frame
		self.score = float(0.0)
		self.hitCoverage = 0 #number of blast hits that cover this feature, reset for each pancake
		self.token = ';'.join([self.seqID,self.type,str(self.position),str(self.frame)])
		
	def isStart(self):
		if self.type.find('start')> -1:
			return 1
		else:
			return 0
			
	def isGap(self):
		if self.type.find('gap') > -1:
			return 1
		else:
			return 0
			
	def isGapInit(self):
		if self.type.find('gap_init') > -1:
			return 1
		else:
			return 0
	
	def isGapTerm(self):
		if self.type.find('gap_term') > -1:
			return 1
		else:
			return 0
			
	def isStop(self):
		if self.type.find('stop')> -1:
			return 1
		else:
			return 0
			
	def isSequenceEnd(self):
		if self.type.find('sequence_end') > -1:
			return 1
		else:
			return 0
			
	def isUpstream(self):
		if self.type.find('upstream') > -1:
			return 1
		else:
			return 0
			
	def isDownstream(self):
		if self.type.find('downstream') > -1:
			return 1
		else:
			return 0