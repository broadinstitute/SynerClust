#!/usr/bin/env python
import sys, os, pickle, TreeLib
import scipy.spatial.distance as distance

def usage():
	print """Launches the series of commands required to create a consensus sequence from all contributing leaf node genes.
	Runs external programs Muscle and FastTree
	
	WF_MakeConsensusSeq_leaf_centroid.py [peptide seqs]
	"""


def main(argv):
	temp_pep = argv[0]
	cluster_dir = "/".join(temp_pep.split("/")[0:-1])+"/"
	node_name = temp_pep.split("/")[-3]
	clusterID = temp_pep.split("/")[-1].split(".")[0]

	temp_pep_align = cluster_dir+clusterID+".mus"
	temp_align = clusterID+".mus"
	consensus_pep = cluster_dir+clusterID+".cons.pep"
	
	muscle_cmd = "#MUSCLE_PATH -maxiters 1 -diags -sv -distance1 kbit20_3  -in "+temp_pep+" -out "+temp_pep_align

	os.system(muscle_cmd)

	mus_out = open(temp_pep_align,'r').readlines()
	mus_seqs = {}
	total_length = 0
	for mo in mus_out:
		l = mo.rstrip()
		if l.find('>') > -1:
			line = l.split()
			seqID = line[0][1:]
			mus_seqs[seqID] = []
		else:
			for i in l:
				mus_seqs[seqID].append(i)
			total_length+=len(l)
			
	min_dist = total_length
	min_dist_seq = ""
	for ms in mus_seqs:
		dist = 0
		for ns in mus_seqs:
			dist += distance.hamming(mus_seqs[ms], mus_seqs[ns])
		if dist < min_dist:
			min_dist = dist
			min_dist_seq = ms
	print "we choose", min_dist_seq, min_dist
	cons_out = open(consensus_pep,'w')
	cons_seq = "".join(mus_seqs[min_dist_seq])
	cons_seq = cons_seq.replace("-","")
	cons_out.write(">"+clusterID+";"+str(len(cons_seq))+"\n")
	cons_out.write(cons_seq+"*\n")

	#~ os.system("rm "+temp_pep_align)

	
	sys.exit(0)

if __name__ == "__main__":
	if len(sys.argv) == 1:
		sys.exit(usage())
	else:
		main(sys.argv[1:])