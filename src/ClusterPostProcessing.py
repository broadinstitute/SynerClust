#!/usr/bin/env python
import sys, os, pickle, logging

def usage():
	print "USAGE: ClusterPostProcessing.py [path to genomes/ in synergy2 working dir] [path to locus_mapping.pkl of node to be processed] [number of leaf genomes descendent of node to be processed]"

def main(argv):
	FORMAT = "%(asctime)-15s %(levelname)s %(module)s.%(name)s.%(funcName)s at %(lineno)d :\n\t%(message)s\n"
	logger = logging.getLogger()
	logging.basicConfig(filename='ClusterPostProcessing.log', format = FORMAT, filemode='w', level=logging.DEBUG)
	# add a new Handler to print all INFO and above messages to stdout
	ch = logging.StreamHandler(sys.stdout)
	ch.setLevel(logging.INFO)
	logger.addHandler(ch)
	logger.info('Started')

	#~ node_path = argv[0]
	locus_mapping = argv[1]
	genome_path = argv[0]
	num_genomes = int(argv[2])
	
	pklFile = open(locus_mapping, 'rb')
	locusMap = pickle.load(pklFile)
	pklFile.close()
	
	l_t = {} #locus to transcript
	l_s = {} #locus to sequence
	genomes = os.listdir(genome_path)
	for g in genomes:
		if g.find("locus_tag_file.txt")>-1:
			continue
		dataFile = genome_path+g+"/annotation.txt"
		data = open(dataFile, 'r').readlines()
		for d in data:
			d = d.rstrip()
			line = d.split()
			l_t[line[1]] = line[0]
			l_s[line[1]] = line[7]
	
	totalGenes = 0
	counter = 1
	clusters = {}
	tIDs = set([])
	clusterDist = {}
	cluster_noOrphan = 0
	almost_core = 0
	almost_scc = []
	for l in locusMap:
		clusters[counter] = {'leaves':{},'transcripts':[]}
		leafKids = locusMap[l]
		#~ myFile = open(l+".pep",'w')
		for k in leafKids:
			#~ myFile.write(">"+k+"\n"+l_s[k]+"\n")
			clusters[counter]['transcripts'].append(l_t[k])
			tIDs.add(l_t[k])
			prefix = "_".join(k.split("_")[:-1])
			logger.debug("%s splitted to %s" %(k, prefix))
			if not prefix in clusters[counter]['leaves']:
				clusters[counter]['leaves'][prefix] = 0
			clusters[counter]['leaves'][prefix]+=1
			
		#~ myFile.close()
		if not len(leafKids) in clusterDist:
			clusterDist[len(leafKids)] = 0
		clusterDist[len(leafKids)]+=1
		if len(clusters[counter]['leaves']) ==num_genomes and len(leafKids)>num_genomes:
			dups = []
			keeper = "YES!"
			for p in clusters[counter]['leaves']:
				if clusters[counter]['leaves'][p] > 1:
					dups.append(p)
			for d in dups:
				leafNums = []
				for l in leafKids:
					if l.find(d)>-1:
						num = int(l.split("_")[-1])
						leafNums.append(num)
				leafNums.sort()
				range = leafNums[-1] - leafNums[0] + 1
				avg = float(range)/float(len(leafNums))
				if avg > 1:
					keeper = "no"
					if avg < 3.0:
						keeper = "close"
				#~ print d, leafNums, range, avg, keeper
			if keeper == "YES!":
				almost_core+=1
				almost_scc.append(counter)
		
		totalGenes += len(clusters[counter]['transcripts'])
		if len(clusters[counter]['transcripts']) > 1:
			cluster_noOrphan+= 1
		counter+=1
	print "total genes: ",totalGenes

	
	pairs = set([])
	scc = set([])
	sccPlus = set([])
	multicc = set([])
	scc_count = 0
	mcc_count = 0
	mcc_clusters = []
	count = 0
	dir_split = locus_mapping.split("/")
	dir_split.pop()
	dir_split.append("final_clusters.txt")
	cluster_out = "/".join(dir_split)
	dir_split.pop()
	dir_split.append("clust_to_trans.txt")
	cTt_out = "/".join(dir_split)
	cout = open(cluster_out, 'w')
	ct_out = open(cTt_out, 'w')
	for c in clusters:
		transcripts = clusters[c]['transcripts']
		cid = "Cluster"+str(c)
		t_out = cid+" (taxa: "+str(len(clusters[c]['leaves']))+", genes: "+str(len(transcripts))+")\t"+" ".join(transcripts)+"\n"
		cout.write(t_out)
		for t in transcripts:
			ct_out.write(cid+"\t"+t+"\n")
		if len(transcripts) == 1:
			continue
		else:
			mypairs = []
			count +=1
			for t in transcripts:
				tNum = int(t)
				for s in transcripts:
					if s==t:
						continue
					sNum = int(s)
					tlist = [sNum,tNum]
					tlist.sort()
					tup = (tlist[0],tlist[1])
					mypairs.append(tup)
					pairs.add(tup)
			if not (len(clusters[c]['leaves']) == num_genomes):
				continue
			else:
				mcc_count+=1
				mcc_clusters.append(c)
				for mp in mypairs:
					multicc.add(mp)
				if c in almost_scc:
					for mp in mypairs:
						sccPlus.add(mp)
				if not (len(transcripts) == num_genomes):
					continue
				else:
					scc_count+=1
					for mp in mypairs:
						scc.add(mp)
						sccPlus.add(mp)
	orphan_count = len(clusters) - cluster_noOrphan
	aux_count = cluster_noOrphan - mcc_count
	print "pairs:", len(pairs)
	print "scc:", scc_count
	print "mcc:", mcc_count-scc_count
	print "aux:", aux_count
	print "orphans:", orphan_count
	print "non-orphan clusters:", count
	cout.close()
	ct_out.close()
	
	ds = locus_mapping.split("/")
	ds.pop()
	mydir = "/".join(ds)+"/"
	pair_pkl = mydir+"tuple_pairs.pkl"
	scc_pkl = mydir +"tuple_pairs_scc.pkl"
	sccPlus_pkl = mydir +"tuple_pairs_sccPlus.pkl"
	mcc_pkl = mydir+"tuple_pairs_mcc.pkl"

	pdat = open(pair_pkl, 'wb')
	pickle.dump(pairs, pdat)
	pdat.close()
				
if __name__ == "__main__":
	if len(sys.argv) == 1:
		sys.exit(usage())
	else:
		main(sys.argv[1:])