from string import maketrans


class Genome:
	def __init__(self, file, transl_table=1):
		self.file = file
		self.totalLength = 0
		self.scaffolds = {}
		self.stops = ['TAA', 'TAG', 'TGA']
		self.codons = {'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L', 'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P', 'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E', 'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W', 'CGT': 'R', 'CGA': 'R', 'CGG': 'R', 'CGC': 'R', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R', 'GGT': 'G', 'GGG': 'G', 'GGA': 'G', 'GGC': 'G'}
		if transl_table == 3:
			self.stops = ['TAA', 'TAG']
			self.codons['ATA'] = 'M'
			self.codons['CTT'] = 'T'
			self.codons['CTC'] = 'T'
			self.codons['CTA'] = 'T'
			self.codons['CTG'] = 'T'
			self.codons['TGA'] = 'W'
			self.codons['CGA'] = ' '
			self.codons['CGC'] = ' '
		# elif transl_table == 11:
		# 	pass  # no difference from transl_table == 1 for codons translations and stops, only for initiation which we don't use here

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
			self.totalLength += self.scaffolds[seqID].getSequenceLength()

	def calcGenomeStats(self, scaffGeneCount, geneCount):
		lengths = []
		for s in self.scaffolds:
			gc = 0
			if s in scaffGeneCount:
				gc = scaffGeneCount[s]
			tup = (s, (self.scaffolds[s]).getSequenceLength(), gc)
			lengths.append(tup)
		n50_min = self.totalLength * .5
		n90_min = self.totalLength * .9
		g50_min = geneCount * .5
		g90_min = geneCount * .9
		# print n50_min, n90_min
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
			if gcount >= g90_min and g90 == 0:
				g90 = l[1]
			if gcount >= g50_min and g50 == 0:
				g50 = l[1]
		print "n50", n50
		print "n90", n90
		retlist = [str(self.totalLength), str(n50), str(n90), str(g50), str(g90)]
		return retlist

	def getCDS(self, seqID, start, stop, strand, phase):
		comp_start = start - 1
		comp_stop = stop
		if strand == "+":
			cds = self.scaffolds[seqID].sequence[comp_start + phase:comp_stop]
		elif strand == "-":
			cds = self.scaffolds[seqID].sequence[comp_start:comp_stop - phase]
		else:
			print("Error reading strand from gff3 annotation")
		# cds = self.scaffolds[seqID].sequence[comp_start:comp_stop]
		cds = cds.upper()
		if strand.find('+') > -1:
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
			# do some stuff
			codon = CDS[index:index + 3]
			amino = ''
			if codon.find('N') > -1:
				amino = 'X'
			elif codon not in self.codons:
				amino = 'X'
			else:
				amino = self.codons[codon]
			protein = protein + amino
			index = index + 3
		return protein


class Scaffold:
	def __init__(self, ID, globalStart):
		self.ID = ID
		self.sequence = ""
		self.raw_sequence = ""
		# globalStart and globalStop begin at 0, sequence features begin at 1
		self.globalStart = int(globalStart)
		self.globalStop = 0
		self.seqFeats = []

	def addSequence(self, sequence):
		self.sequence = self.sequence + sequence

	def getSequenceLength(self):
		return len(self.sequence)
