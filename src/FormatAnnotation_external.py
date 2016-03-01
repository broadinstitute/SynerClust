#!/usr/bin/env python
import sys
import SequenceParse
import pickle
import os
# Given a fasta sequence, gff3 file, locus tag and output file, the Synergy2 directory structure is created for this genome


def numberFromIndex(locus, count):
	c = str(count)
	while len(c) < 6:
		c = "0" + c
	tag = locus + "_" + c
	return tag


def makePicklesForSingleGenome(working_dir, genome, node, SYNTENIC_WINDOW):
	gdat = open(working_dir + "genomes/" + genome + "/annotation.txt", 'r').readlines()
	# this file might not be necessary in the long-term, in fact these genome directories may not be necessary at all
	# everything should be pulled out and stored in pickles after the repo files are parsed.  Whatever.
	ndat = open(working_dir + "nodes/" + node + "/" + node + ".pep", 'w')
	# print node
	x = gdat[0].split()[1]
	y = "_".join(x.split("_")[:-1])
	if not y == node:
		print node, y
		print gdat[0]

	# these hashes will be pickles
	genes = {}
	gene_map = {}
	locusToTID = {}
	# this dict is for synteny information, keys are scaffold IDs
	neighbors = {}
	for g in gdat:
		g = g.rstrip()
		l = g.split()
		if len(l) < 8:
			print g
			continue
		if not l[2] in neighbors:
			neighbors[l[2]] = []
		gene_tup = (l[1], int(l[3]), int(l[4]), l[5], int(l[6]))  # scaffold -> locus,start,stop,strand,length
		locusToTID[l[1]] = l[0]
		neighbors[l[2]].append(gene_tup)
		genes[l[1]] = l[7]
		gene_map[l[1]] = [l[1]]
		line = ">" + l[1] + ";" + l[6] + "\n" + l[7] + "\n"
		ndat.write(line)
	ndat.close()

	# make synteny pickle
	makeSyntenyPickle(working_dir, genome, node, neighbors, SYNTENIC_WINDOW)
	# dump pickles!
	pdat = open(working_dir + "nodes/" + node + "/" + node + ".pkl", 'wb')
	pickle.dump(genes, pdat)
	pdat.close()
	ldat = open(working_dir + "nodes/" + node + "/locus_mappings.pkl", 'wb')
	pickle.dump(gene_map, ldat)
	ldat.close()
	tdat = open(working_dir + "genomes/" + genome + "/locusToTranscript.pkl", 'wb')
	pickle.dump(locusToTID, tdat)
	tdat.close()
	pcomp = open(working_dir + "nodes/" + node + "/PICKLES_COMPLETE", 'w')
	pcomp.write("Way to go!")
	pcomp.close()


def makeSyntenyPickle(working_dir, genome, node, neighbors, SYNTENIC_WINDOW):
	# make a pickle!
	gsyn = {}
	nsyn = {}
	MAX_DIST = SYNTENIC_WINDOW
	for n in neighbors:
		# n is a scaffold ID
		genes = neighbors[n]
		genes = sorted(genes, key=lambda tup: tup[1])
		for g in genes:
			# g_i = genes.index(g)
			locus = g[0]
			lend = g[1]
			rend = g[2]
			strand = g[3]
			# length = g[4]
			gsyn[locus] = []
			nsyn[locus] = {}
			nsyn[locus]['neighbors'] = []
			nsyn[locus]['count'] = 1
			gmid = (rend - lend + 1) / 2 + lend
			minmid = lend - MAX_DIST
			maxmid = rend + MAX_DIST
			for h in genes:
				mid = (h[2] - h[1] + 1) / 2 + h[1]
				if h[0] == locus:
					continue
				if mid >= minmid and mid <= maxmid:
					# here, direction AND strand yields stream direction, it's bitwise 'and' or something
					stream = 0  # -1 means upstream, 1 means downstream
					dist = -1
					direction = None
					if h[1] > lend:
						dist = mid - gmid + 1
						direction = "+"
					else:
						dist = gmid - mid + 1
						direction = "-"
					if direction == strand:
						stream = -1
					else:
						stream = 1
					genome_tup = (h[0], h[1], h[2], dist, stream)
					gsyn[locus].append(genome_tup)
					nsyn[locus]['neighbors'].append(h[0])
				if mid > maxmid:
					break
# TODO inspect here
# List of all genes. For each gene, list of all other genes in the form (gene, left end, right end, distance, stream(?))
	gdat = open(working_dir + "genomes/" + genome + "/synteny_data.pkl", 'wb')
	pickle.dump(gsyn, gdat)
	gdat.close()
	ndat = open(working_dir + "nodes/" + node + "/synteny_data.pkl", 'wb')
	pickle.dump(nsyn, ndat)
	ndat.close()


def extractAnnotation(gff3_file, seq_file, genome_name, locus, out_file, stat_file, peptide_file):
	# get fasta sequence into a SequenceParse object
	myGenome = SequenceParse.Genome(seq_file)
	myGenome.readSequence()

	# process annotation file
	gff3 = open(gff3_file, 'r').readlines()
	out = open(out_file, 'w')

	pep_hash = {}
	if not peptide_file == "null":
		peps = open(peptide_file, 'r').readlines()
		for p in peps:
			p = p.rstrip()
			if p.find(">") == 0:
				pid = p[1:]
				pep_hash[pid] = ""
			else:
				pep_hash[pid] += p

	count = 1
	curID = ""
	myCDS = ""
	cds_tups = []
	gstart = ""
	gstop = ""
	gstrand = ""
	scaffold = ""
	scaffGeneCount = {}
	for g in gff3:
		g = g.rstrip()
		line = g.split()
		if len(line) < 7:
			continue
		# for every "gene" in the gff3, get the ID= field, this is required, this becomes curID
		if line[2] == "gene":
			# deal with last gene feature
			if len(cds_tups) > 0:
				# concatenate all of the CDS sequences together, then translate them with SequenceParse
				if cds_tups[0][3] == "+":
					cds_tups = sorted(cds_tups, key=lambda feat: feat[1])
				else:
					cds_tups = sorted(cds_tups, key=lambda feat: feat[1], reverse=True)
				for c in cds_tups:
					myCDS = myCDS + c[4]
				if curID in pep_hash:
					peptide = pep_hash[curID]
				else:
					peptide = myGenome.translateCDS(myCDS)
				length = str(len(peptide))

				if scaffold not in scaffGeneCount:
					scaffGeneCount[scaffold] = 0
				scaffGeneCount[scaffold] += 1

				# assign locus tag
				locus_tag = numberFromIndex(locus, count)
				count += 1

				# write to file
				outline = "\t".join([curID, locus_tag, scaffold, gstart, gstop, gstrand, length, peptide]) + "\n"
				out.write(outline)
			myCDS = ""
			cds_tups = []
			data = line[8].split(";")
			gstart = line[3]
			gstop = line[4]
			gstrand = line[6]
			for d in data:
				if d.find("ID=") == 0:
					curID = d.split("=")[1]
		# for every subsequence "exon" after the "gene", get the CDS using SequenceParse methods
		elif line[2] == "exon":
			scaffold = line[0]
			start = int(line[3])
			stop = int(line[4])
			strand = line[6]
			next_cds = "junk!"
			if curID not in pep_hash:
				next_cds = myGenome.getCDS(scaffold, start, stop, strand)
			my_tup = (scaffold, int(start), stop, strand, next_cds)
			cds_tups.append(my_tup)
	if len(cds_tups) > 0:
		# concatenate all of the CDS sequences together, then translate them with SequenceParse
		if cds_tups[0][3] == "+":
			cds_tups = sorted(cds_tups, key=lambda feat: feat[1])
		else:
			cds_tups = sorted(cds_tups, key=lambda feat: feat[1], reverse=True)
		for c in cds_tups:
			myCDS = myCDS + c[4]
		peptide = myGenome.translateCDS(myCDS)
		length = str(len(peptide))

		# assign locus tag
		locus_tag = numberFromIndex(locus, count)
		if scaffold not in scaffGeneCount:
			scaffGeneCount[scaffold] = 0
		scaffGeneCount[scaffold] += 1
		# write to file
		outline = "\t".join([curID, locus_tag, scaffold, gstart, gstop, gstrand, length, peptide]) + "\n"
		out.write(outline)
	out.close()

	stats = myGenome.calcGenomeStats(scaffGeneCount, count)  # [length, scaff n50, scaff n90, gene n50, gene n90]
	stat_out = open(stat_file, 'w')
	# outline = Genome\tlocus\tlength\tN50\tN90\tgene_count
	stats.insert(0, genome_name)
	stats.insert(1, locus)
	stats.append(str(count))
	outline = "\t".join(stats) + "\n"
	stat_out.write(outline)
	print outline
	stat_out.close()


def main(argv):
	# gff3, seq, and genome name
	gff3_file = argv[0]
	seq_file = argv[1]
	genome_name = argv[2]
	peptide_file = argv[3]
	locus = argv[4]
	out_file = argv[5]
	SYNTENIC_WINDOW = int(argv[6])
	annot = argv[7]
	pickle_juice = argv[8]
	print argv

	if annot == "1":
		stat_file = out_file.replace("annotation.txt", "stats.txt")
		extractAnnotation(gff3_file, seq_file, genome_name, locus, out_file, stat_file, peptide_file)

	if pickle_juice == "1":
		working = out_file.split("genomes")[0]
		nodes_dir = working + "nodes/"
		node_dir = nodes_dir + locus
		if locus not in os.listdir(nodes_dir):
			os.system("mkdir " + node_dir)
		node_dir = node_dir + "/"
		makePicklesForSingleGenome(working, genome_name, locus, SYNTENIC_WINDOW)


if __name__ == "__main__":
	if len(sys.argv) == 0:
		sys.exit("USAGE: python FormatAnnotation.py annotation.gff3 sequence.fa locus_tag /path/to/annotation/file.txt")
	else:
		main(sys.argv[1:])
