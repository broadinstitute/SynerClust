#!/usr/bin/env python
import sys, os

def usage():
	print """This script sets up Synergy2 to run on your system.
	"""

def main(argv):
	#there shouldn't be arguments
	
	#variables that will be set:
	s2path = ""
	s2_bin = ""
	s2_src = ""
	QT_path = ""
	muscle_path = ""
	blast_path = ""
	grid_local_exec = ""
	wf_type = ""
	wf_queue = ""
	wf_execenv = ""
	#from standard in, get various stuff
	#SYNERGY2_PATH
	s2path = os.getcwd()+"/"
	print s2path
	os.system("mkdir "+s2path+"bin")
	s2_bin = s2path+"bin/"
	s2_src = s2path+"src"
	#QUICKTREE_PATH
	FT_path = raw_input("Enter the path to your QuickTree executable:\n")
	#MUSCLE_PATH
	muscle_path = raw_input("Enter the path to your MUSCLE executable:\n")
	#BLAST_PATH
	# blast_path = raw_input("Enter the path to your blastall executable:\n")
	blast_path = raw_input("Enter the path to your blast installation folder:\n")

	#~ print "The following variables deal with grid submission via WorkFlow."
	#~ grid_local_exec = raw_input("Will jobs be submitted to a grid or executed locally? (grid|local): \n")
	#~ if grid_local_exec.find("grid")>-1:
		#WORKFLOW_TYPE
		#~ wf_type = raw_input("Enter the grid type that workflow will submit jobs to (LSF|SGE|Condor): \n")
		#WORKFLOW_QUEUE
		#~ wf_queue = raw_input("Enter the queue name that workflow should submit jobs to: \n")
	#WORKFLOW_EXECENV
	#~ wf_execenv = raw_input("Enter the path to the exec env workflow should use: \n")

	print "The following variables are defaults, but can be changed for any execution of Synergy2."
	#BLAST_EVAL_DEFAULT
	eval_default = raw_input("Enter the default e-value for blast; this is supposed to be a relaxed cut-off. Suggested=0.01:\n")
	#NUM_CORES_DEFAULT
	num_cores_default = raw_input("Enter the number of cores that blastall should use (-a flag): ")
	
	src_files = os.listdir(s2_src)
	for sfile in src_files:
		sf = open(s2_src+"/"+sfile,'r').read()
		sf = sf.replace('#SYNERGY2_PATH',s2_bin)
		sf = sf.replace('#FASTTREE_PATH',FT_path)
		sf = sf.replace('#MUSCLE_PATH',muscle_path)
		sf = sf.replace('#BLAST_PATH',blast_path)
		#~ sf = sf.replace('#WORKFLOW_TYPE',wf_type)
		#~ if grid_local_exec.find("local")>-1:
			#~ sf = sf.replace('remote-', '')
		#~ sf = sf.replace('#WORKFLOW_EXECENV',wf_execenv)
		
		sf = sf.replace('#BLAST_EVAL_DEFAULT',eval_default)
		sf = sf.replace('#NUM_CORES_DEFAULT',num_cores_default)

		mod_sf = s2_bin+sfile
		msf= open(mod_sf,'w')
		msf.write(sf)
		msf.close()
		
		os.system("chmod ugo+x "+mod_sf)
		
	print "Synergy2 has been installed. It can be run from ~/bin/WF_runSynergy2.py.  See the README for additional information."
				
if __name__ == "__main__":
	#~ if len(sys.argv) == 1:
		#~ sys.exit(usage())
	#~ else:
	main(sys.argv[1:])