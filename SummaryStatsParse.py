#!/usr/bin/env python
import sys, math

def usage():
	print """SummaryStatsParse.py [summary stats file] [number of genomes]"""
	
def main(argv):
	sstats_file = argv[0]
	num_genomes = float(argv[1])
	
	sstat = open(sstats_file,'r').readlines()
	cumSD = 0.0
	cumSDoverMean = 0.0
	cumScore = 0.0
	totalRows = len(sstat)
	nonOrphans = len(sstat)
	for ss in sstat:
		ss = ss.rstrip()
		line = ss.split()
		cluster = line[0]
		genes = float(line[1])
		taxa = float(line[2])
		min_taxa = float(line[3])
		max_taxa = float(line[4])
		std_dev = float(line[8])
		#~ len_mean = float(line[5])
		sd_over_mean = float(line[9])
		
		#~ t_ratio = min_taxa/max_taxa   #simple ratio
		t_ratio = min_taxa/(genes/taxa) #more differentiation ratio
		SD_factor = 1.0 - sd_over_mean
		#~ print line
		#~ print genes, taxa, min_taxa, max_taxa, sd_over_mean, t_ratio, SD_factor
		#~ sys.exit()
		score = 0.0
		#~ score = t_ratio*SD_factor*(genes/num_genomes)
		score = t_ratio*SD_factor*genes
		
		if genes==1.0 or taxa==1.0:
			nonOrphans -= 1
			cumScore+=score
			continue
		#~ elif sd_over_mean>= 0.5:
			#~ print "\t".join([line[0],line[1],line[2],line[5],line[7],str(score)])
			#~ nonOrphans -=1
			#~ continue
		cumScore+=score
		cumSD += std_dev
		cumSDoverMean += sd_over_mean
		
	avg_SD = cumSD/float(nonOrphans)
	avg_SD_Mean = cumSDoverMean/float(nonOrphans)
	print "rows: total, valid", totalRows, nonOrphans
	print "score", cumScore
	print "avg std dev", avg_SD
	print "avg std dev/mean length", avg_SD_Mean
	
	#~ flaggedSD = 3.0*avg_SD_Mean
	#~ flaggedSD = 0.5
	#~ for ss in sstat:
		#~ ss = ss.rstrip()
		#~ line = ss.split()
		#~ if float(line[7])>flaggedSD:
			#~ print ss
		
	
if __name__ =="__main__":
	if len(sys.argv)==1:
		sys.exit(usage())
	else:
		main(sys.argv[1:])