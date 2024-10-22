#!/usr/bin/env python

import sys
import SequenceParse
import pickle
import os
import numpy
import argparse
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
	x = gdat[1].split("\t")[1]
	y = "_".join(x.split("_")[:-1])
	if not y == node:
		print node, y
		print gdat[1]
	# these hashes will be pickles
	genes = {}
	gene_map = {}
	locusToTID = {}
	# this dict is for synteny information, keys are scaffold IDs
	neighbors = {}
	for g in gdat[1:]:
		g = g.rstrip()
		l = g.split("\t")
		if len(l) < 8:
			print g
			continue
		if not l[2] in neighbors:
			neighbors[l[2]] = []
		gene_tup = (l[1], int(l[3]), int(l[4]), l[5], int(l[6]))  # scaffold -> generated locus name ,start,stop,strand,length
		locusToTID[l[1]] = l[0]
		neighbors[l[2]].append(gene_tup)
		genes[l[1]] = l[10]
		gene_map[l[1]] = [l[1]]
		line = ">" + l[1] + ";" + l[6] + "\n" + l[10] + "\n"
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
	open(working_dir + "nodes/" + node + "/NODE_COMPLETE", 'a').close()


def makeSyntenyPickle(working_dir, genome, node, neighbors, SYNTENIC_WINDOW):
	# make a pickle!
	gsyn = {}
	nsyn = {}
	# +1 in size of the array because gene id counter starts at 0, so for later it's easier to just use that id without -1, leaving tsyn[0] empty
	# tsyn = numpy.empty(sum([len(neighbors[i]) for i in neighbors]) + 1, dtype=list)
	# tsyn = numpy.empty(int(max([g[len(g) - 1][0] for g in [neighbors[n] for n in neighbors]]).split("_")[-1]) + 1, dtype=list)
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
			# tlocus = int(locus.split("_")[-1])
			# tsyn[tlocus] = []
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
					# tsyn[tlocus].append(h[0])
				if mid > maxmid:
					break
# TODO inspect here
# List of all genes. For each gene, list of all other genes in the form (gene, left end, right end, distance, stream(?))
	gdat = open(working_dir + "genomes/" + genome + "/synteny_data.pkl", 'wb')
	pickle.dump(gsyn, gdat)
	gdat.close()
	gdat = open(working_dir + "nodes/" + node + "/synteny_data.pkl", 'wb')
	pickle.dump(nsyn, gdat)
	gdat.close()
	# gdat = open(working_dir + "nodes/" + node + "/synteny_table.pkl", 'wb')
	# pickle.dump(tsyn, gdat)
	# gdat.close()


def extractAnnotation(gff3_file, seq_file, genome_name, locus, out_file, stat_file, peptide_file):

	out = open(out_file, 'w')
	pep_hash = {}
	count = 1

	if peptide_file != "null":
		out.write(genome_name + "\t" + peptide_file[peptide_file.rfind("/") + 1:] + "\n")
		peps = open(peptide_file, 'r').readlines()
		for p in peps:
			p = p.rstrip()
			if p.find(">") == 0:
				pid = p[1:]
				npid = pid.split("\t")[0].split(":")[1]
				# for tmp in pid.split():
				# 	if tmp[:5] == "Gene:":
				# 		npid = tmp[5:]
				# 	elif tmp[:8] == "ENSEMBL:":
				# 		npid = tmp[8:]
				pep_hash[npid] = ""
			else:
				pep_hash[npid] += p
		for pid in pep_hash:
			outline = "\t".join([pid, numberFromIndex(locus, count), numberFromIndex(genome_name, count), "0", str((len(pep_hash[pid]) * 3) - 1), ".", str(len(pep_hash[pid])), pid, "None", "None", pep_hash[pid]]) + "\n"
			out.write(outline)  # add a buffer that gets printed every modulo x of count, output buffer at the end too
			count += 1

	else:
		# process annotation file
		gff3 = open(gff3_file, 'r').readlines()
		# get fasta sequence into a SequenceParse object
		myGenome = SequenceParse.Genome(seq_file)
		myGenome.readSequence()

		out.write(genome_name + "\t" + gff3_file[gff3_file.rfind("/") + 1:] + "\n")

		curID = ""
		myCDS = ""
		cds_tups = []
		gstart = ""
		gstop = ""
		gstrand = ""
		scaffold = ""
		scaffGeneCount = {}
		name = "None"
		mrna_id = ""
		alias = None
		for g in gff3:
			g = g.rstrip()
			line = g.split("\t")
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

					if mrna_id is "":
						mrna_id = curID
					if alias is None:
						alias = curID

					# write to file
					outline = "\t".join([curID, locus_tag, scaffold, gstart, gstop, gstrand, length, mrna_id, alias, name, peptide]) + "\n"
					out.write(outline)
				myCDS = ""
				cds_tups = []
				data = "-".join(line[8:]).split(";")   # join to read the whole attributes even if spaces
				gstart = line[3]
				gstop = line[4]
				gstrand = line[6]
				name = "None"
				mrna_id = ""
				alias = None
				# cds_tups_primary = False
				# primary_done = False
				for d in data:
					d = d.split(",")[0]
					if d.lower().find("id=") == 0:
						curID = d.split("=")[1]
					if d.lower().find("dbxref=") == 0:  # added for NCBI's genbank in gff3 format
						curID = d.split(":")[1]
					if d.lower().find("alias") == 0:
						alias = d.split("=")[1]
					if d.lower().find("name=") == 0:
						name = d.split("=")[1]
					if d.lower().find("locus_tag=") == 0:
						name = d.split("=")[1]
			# for every subsequence "exon" after the "gene", get the CDS using SequenceParse methods
			# elif line[2] == "exon":
			elif line[2] == "CDS":
				scaffold = line[0]
				start = int(line[3])
				stop = int(line[4])
				strand = line[6]
				if line[7] == ".":
					phase = 0
				else:
					phase = int(line[7])
				next_cds = "junk!"
				if curID not in pep_hash:
					next_cds = myGenome.getCDS(scaffold, start, stop, strand, phase)
				# if not primary_done:
				# 	next_cds = myGenome.getCDS(scaffold, start, stop, strand, phase)
				my_tup = (scaffold, int(start), stop, strand, next_cds)
				cds_tups.append(my_tup)
			elif line[2] == "mRNA":
				# if cds_tups_primary:
				# 	primary_done = True
				# if not primary_done:
				for d in line[8].split(";"):
					if d.lower().find("id=") == 0:
						mrna_id_tmp = d.split("=")[1]
					if d.lower().find("dbxref=") == 0:
						mrna_id_tmp = d.split(":")[1]
					# 	if d.find("VectorBase:") == 0:  # import that in first position for this one as it can also be in the Dbxref field
					# 		if d[-3:] == "-RA":
					# 			cds_tups_primary = True
					# 			mrna_id = mrna_id_tmp
					# 		cds_tups = []
					# if not cds_tups_primary:
					# 	mrna_id += mrna_id_tmp
				mrna_id += mrna_id_tmp
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

			if (length == 0):
				sys.exit("Empty peptide sequence, please make sure coordinates from the GFF3 file match to the FASTA file.")

			# assign locus tag
			locus_tag = numberFromIndex(locus, count)
			if scaffold not in scaffGeneCount:
				scaffGeneCount[scaffold] = 0
			scaffGeneCount[scaffold] += 1

			if mrna_id == "":
				mrna_id = curID
			if alias is None:
				alias = curID

			# write to file
			outline = "\t".join([curID, locus_tag, scaffold, gstart, gstop, gstrand, length, mrna_id, alias, name, peptide]) + "\n"
			out.write(outline)

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
	out.close()


def main(argv):
	usage = "usage: FormatAnnotation_external.py [options] <reference.fasta> <scaffolds.fasta>"
	parser = argparse.ArgumentParser(usage)
	parser.add_argument('-gff', dest="gff3_file", help="GFF3 annotation file.")
	parser.add_argument('-seq', dest="seq_file", help="Sequence file.")
	parser.add_argument('-name', dest="genome_name", required=True, help="Genome name. (Required)")
	parser.add_argument('--pep', dest="peptide_file", help="Peptide file name. (Overrides Sequence + GFF3, no synteny)")
	parser.add_argument('-locus', dest="locus", required=True, help="Locus name. (Required)")
	parser.add_argument('-out', dest="out_file", required=True, help="Output name. (Required)")
	parser.add_argument('-synteny', type=int, dest="syntenic_window", default=6000, help="Syntenic window size. (Required)")
	parser.add_argument('--annot', dest="annot", default='0', help="")
	parser.add_argument('--pickle', dest="pickle_juice", default='0', help="")
	parser.add_argument("--transl_table", type=int, dest="transl_table", default=1, help="Translation table to use. (NCBI IDs)")
	args = parser.parse_args()

	if args.annot == "1":
		stat_file = args.out_file.replace("annotation.txt", "stats.txt")
		extractAnnotation(args.gff3_file, args.seq_file, args.genome_name, args.locus, args.out_file, stat_file, args.peptide_file)

	if args.pickle_juice == "1":
		working = args.out_file.split("genomes")[0]
		nodes_dir = working + "nodes/"
		node_dir = nodes_dir + args.locus
		if args.locus not in os.listdir(nodes_dir):
			os.system("mkdir " + node_dir)
		node_dir = node_dir + "/"
		makePicklesForSingleGenome(working, args.genome_name, args.locus, args.syntenic_window)


if __name__ == "__main__":
	if len(sys.argv) == 0:
		sys.exit("USAGE: python FormatAnnotation.py annotation.gff3 sequence.fa locus_tag /path/to/annotation/file.txt")
	else:
		main(sys.argv[1:])
