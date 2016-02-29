#!/usr/bin/env python

import sys
import os
# import time
import pickle
import logging
import threading
from multiprocessing import Process, Queue
from Queue import Empty
import scipy.spatial.distance as distance
import subprocess
# import io
# import shlex

# manifest_complete = 0

dev_null = open(os.devnull, 'w')

def usage():
	print "After consensus sequences are generated, concatenates them into one big file, NODE.pep, and creates a NODE_COMPLETE file."
	sys.exit(1)
	
def run_muscle(cmd, stdin_data):
# 	muscle_cmd =["#MUSCLE_PATH", "-maxiters", "2", "-diags", "-sv", "-distance1", "kbit20_3", "-quiet"]
# 	muscle_cmd =["/home/kamigiri/tools/muscle3.8.31_i86linux64", "-maxiters", "2", "-diags", "-sv", "-distance1", "kbit20_3", "-quiet"]
# 	cmd = '/home/kamigiri/tools/muscle3.8.31_i86linux64 -maxiters 2 -diags -sv -distance1 kbit20_3 -quiet'
	process = subprocess.Popen(cmd, stdin = subprocess.PIPE, stdout = subprocess.PIPE, stderr = dev_null)
# 	process = subprocess.Popen(shlex.split(cmd), stdin = subprocess.PIPE, stdout = subprocess.PIPE, stderr = dev_null)
	output = process.communicate(stdin_data)[0]
	
# 	process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)

# 	qr, qw = os.fdopen(r,'r',0), os.fdopen(w,'w',0)
	
# 	qr2, qw2 = os.fdopen(r2,'r',0), os.fdopen(w2,'w',0)
# 	w.write(stdin_data)
# 	w.close()
# 	w.flush()
# 	process = subprocess.Popen(cmd, shell = True, stdout = w2, stdin = r)#, stderr = dev_null)
# 	process.wait()
# 	qw2.close()
# 	qr.close()
# 	output = r2.readlines()
# 	qr2.close()

	# r, w = os.fdopen(r,'r',0), os.fdopen(w,'w',0)
# 	r, w = os.pipe()
# 	r2, w2 = os.pipe()
# 	ow = os.fdopen(w, 'w', 0)
# 	orr = os.fdopen(r, 'r', 0)
# 	ow2 = os.fdopen(w2, 'w', 0)
# 	or2 = os.fdopen(r2, 'r', 0)
# 	# r2, w2 = os.fdopen(r2,'r',0), os.fdopen(w2,'w',0)
# 	ow.write(stdin_data)
# 	ow.close()
# 	cmd = '/home/kamigiri/tools/muscle3.8.31_i86linux64 -maxiters 2 -diags -sv -distance1 kbit20_3 -quiet <&{fd1} >&{fd2}'.format(fd1=r, fd2=w2)
# 	print stdin_data[:100] + "\n"
# 	process = subprocess.call(cmd, shell = True, stderr = dev_null)
# # 	process.wait()
# 	orr.close()
# 	ow2.close()
# 	output = or2.readlines()
# 	or2.close()

	return output
# 	process = subprocess.Popen(cmd, stdout = subprocess.PIPE, stdin = subprocess.PIPE, stderr = dev_null)
# 	stdout_stderr = process.communicate(stdin_data)
# 	process.wait()
# 	return stdout_stderr[0]
# 	return process.stdout.readlines()

def makeConsensus(tq, hamm_dist, consensus_pep, output_lock):
	while True:
		try:
			nos = tq.get(block=True, timeout=3)
			# ~ print tq.qsize()
# 			for no in nos:
			pep_data = nos[0]
			clusterID = nos[1]
# 				temp_pep = no[0]
# 				cluster_dir = "/".join(temp_pep.split("/")[0:-1]) + "/"
# 				node_name = temp_pep.split("/")[-3]
# 				clusterID = temp_pep.split("/")[-1].split(".")[0]
			
# 				temp_pep_align = cluster_dir + clusterID + ".mus"
# 				temp_align = clusterID + ".mus"
# 				consensus_pep = cluster_dir + clusterID + ".cons.pep"
			
# 				muscle_cmd = "#MUSCLE_PATH -maxiters 2 -diags -sv -distance1 kbit20_3 -quiet -in "+temp_pep+" -out "+temp_pep_align
# 				muscle_cmd = "#MUSCLE_PATH -maxiters 2 -diags -sv -distance1 kbit20_3 -quiet -in /dev/stdin -out /dev/stdout"  # < " + pep_data
# 			muscle_cmd = "~/tools/muscle3.8.31_i86linux64 -maxiters 2 -diags -sv -distance1 kbit20_3 -quiet -in /dev/stdin -out /dev/stdout"  # < " + pep_data
# 			muscle_cmd = ["#MUSCLE_PATH", "-maxiters", "2", "-diags", "-sv", "-distance1", "kbit20_3", "-quiet", "-in", "/dev/stdin", "-out", "/dev/stdout"]
			muscle_cmd =["#MUSCLE_PATH", "-maxiters", "2", "-diags", "-sv", "-distance1", "kbit20_3", "-quiet"]
# 			muscle_cmd =["/home/kamigiri/tools/muscle3.8.31_i86linux64", "-maxiters", "2", "-diags", "-sv", "-distance1", "kbit20_3", "-quiet", "-in", "/dev/stdin", "-out", "/dev/stdout"]
# 				os.system(muscle_cmd)
# 				mus_out = open(temp_pep_align, 'r').readlines()
			mus_out = run_muscle(muscle_cmd, "".join(pep_data)).split("\n")
			mus_seqs = {}
			mus_str_seqs = {}
			total_length = 0
			unrep = []
			for mo in mus_out:
				l = mo.rstrip()
				if l.find('>') > -1:
					line = l.split()
					seqID = line[0][1:]
					trailer = 1
					while seqID in mus_seqs:
						seqID = line[0][1:] + "." + str(trailer)
						trailer += 1
					mus_seqs[seqID] = []
					unrep.append(seqID)
					mus_str_seqs[seqID] = ""
				else:
					mus_str_seqs[seqID] += l
					for i in l:
						mus_seqs[seqID].append(i)
					total_length += len(l)
					
# 				cons_out = open(consensus_pep,'w')
			output_lock.acquire()
			cons_out = open(consensus_pep, "a")
			while len(unrep) > 0:
				hams = {}
				minDist = float(len(unrep))
				minDistSeq = ""
				for ms in unrep:
					hams[ms] = {}
					m_m = mus_str_seqs[ms]
					m_len = len(m_m)
					m_start = m_len - len(m_m.lstrip("-"))
					m_stop = len(m_m.rstrip("-"))
					dist = 0
# 						mod_dist = 0
					for ns in unrep:
						hams[ms][ns] = distance.hamming(mus_seqs[ms][m_start:m_stop], mus_seqs[ns][m_start:m_stop])
						dist += distance.hamming(mus_seqs[ms], mus_seqs[ns])
					
					if dist < minDist or len(minDistSeq) == 0:
						minDist = dist
						minDistSeq = ms
				cons_seq = "".join(mus_seqs[minDistSeq])
				cons_seq = cons_seq.replace("-", "")
				cons_out.write(">" + clusterID + ";" + str(len(cons_seq)) + "\n")
				cons_out.write(cons_seq + "*\n")
				repped = []
				for ms in unrep:
					if hams[ms][minDistSeq] < hamm_dist:
						repped.append(ms)
				for r in repped:
					unrep.remove(r)
			cons_out.close()
			output_lock.release()
# 				os.system("rm " + temp_pep_align)
			#######
		except Empty:
			break

if __name__ == "__main__":
	argv = []
	if len(sys.argv) == 1:
		usage()
	else:
		argv = sys.argv[1:]
	# after consensus sequences are generated, concatenates them into one big file, NODE.pep, and creates a NODE_COMPLETE file
	node_dir = argv[0]
	node = argv[1]
	fileDir = node_dir + "clusters/"
	nThreads = int(argv[2])
	numThreads = min(nThreads, 4)
	chunk = 25
	h_dist = 0.6
	if len(argv) > 3:
		h_dist = float(argv[3])
	
	FORMAT = "%(asctime)-15s %(levelname)s %(module)s.%(name)s.%(funcName)s at %(lineno)d :\n\t%(message)s\n"
	logger = logging.getLogger()
	logging.basicConfig(filename=node_dir + 'FinalizeNode_threaded.log', format=FORMAT, filemode='w', level=logging.DEBUG)
	# add a new Handler to print all INFO and above messages to stdout
	ch = logging.StreamHandler(sys.stdout)
	ch.setLevel(logging.INFO)
	logger.addHandler(ch)
	logger.info('Started')

	pklPep = node_dir + "pep_data.pkl"
	sdat = open(pklPep, 'rb')
	pickleSeqs = pickle.load(sdat)
	sdat.close()


	# CLEAN UP
# 	manifest = {}
# 	files = os.listdir(fileDir)
# 	
# 	for f in files:
# 		id = f.split(".")[0]
# 		if not id in manifest:
# 			manifest[id] = []
# 		manifest[id].append(f)
	notOKQ = Queue(0)
	output_lock = threading.Lock()
	consensus_pep = node_dir + node + ".pep"
	os.system("rm " + consensus_pep)  # since the file is opened in append mode every time, we need to delete anything from a previous run
	
# 	numThreads = 1
	processes = [Process(target=makeConsensus, args=(notOKQ, h_dist, consensus_pep, output_lock)) for i in range(numThreads)]
# 	notOK = []
# 	processID = 1

	for clusterID in pickleSeqs:
		notOKQ.put((pickleSeqs[clusterID], clusterID))

# 	my_chunk = []
# 	for m in manifest:
# 		ok = 0
# 		for f in manifest[m]:
# 			if f.find("cons.pep")>-1:
# 				ok=1
# 		if ok==0:
# 			if len(my_chunk) == chunk:
# 				notOKQ.put(my_chunk)
# 				my_chunk = []
# 			myFile =  fileDir+m+".pep"
# # 			notOK.append((myFile,m))
# 			my_chunk.append((myFile,m))  # file_name, ID
# 	if len(my_chunk) >0:
# 		notOKQ.put(my_chunk)

	for p in processes:
		p.start()
		logger.info("Starting %s" % (p.pid))
		# print "Starting",p.pid
	for p in processes:
		logger.info("Starting %s" % (p.pid))
		# print "Stopping",p.pid
		p.join()
	
	logger.info("All consensus sequences accounted for.")
	# print "All consensus sequences accounted for."
# 	time.sleep(10)  # allows all files to get synced up
	
# 	filePattern = "cons.pep"
# 	chunk = 500  # number of consensus sequence to concatenate at once
	
# 	dirFiles = os.listdir(fileDir)
# 	matched_files = []
# 	for df in dirFiles:
# 		if df.find(filePattern) > -1:
# 			matched_files.append(fileDir + df)
	
# 	temp_cats = []
# 	index = 0
# 	while len(matched_files) > 0:
# 		temp = matched_files[0:chunk]
# 		del matched_files[0:chunk]
# 		thisCat = fileDir + "tempCat" + str(index) + ".tc.pep"
# 		temp_cats.append(thisCat)
# 		file_string = " ".join(temp)
# 		os.system("cat " + file_string + " > " + thisCat)
# 		index += 1
# 	
# 	os.system("cat " + node_dir + "clusters/*.tc.pep >" + node_dir + node + ".pep")
# 	os.system("rm " + node_dir + "clusters/*.tc.pep")
# 	os.system("rm -r " + node_dir + "clusters")

	os.system("cat " + node_dir + "clusters/singletons.cons.pep >> " + consensus_pep)  # ">>" to append 
	logger.info("Made pep file")
	# print "Made pep file, sleeping 5 seconds"
# 	time.sleep(5)
	
	sys.exit()
	node_done_file = node_dir + "NODE_COMPLETE"
	cr = open(node_done_file, 'w')
	cr.write("Way to go!\n")
	cr.close()
	
