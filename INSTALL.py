#!/usr/bin/env python

import os
import argparse
import platform


def usage():
	print """This script sets up SynerClust to run on your system and compiles FastTree if needed.
	"""


def main():
	usage = "usage: INSTALL.py [options]"
	parser = argparse.ArgumentParser(usage)
	parser.add_argument('-b', dest="blast", default="", help="Path to BLAST+ bin directory (only needed if not in the PATH")
	parser.add_argument('-e', dest="evalue", default="0.01", help="Current node name. (Required)")
	parser.add_argument('-n', dest="threads", default="4", help="Current node name. (Required)")
	args = parser.parse_args()

	scpath = os.getcwd() + "/"
	os.system("mkdir " + scpath + "bin")
	sc_bin = scpath + "bin/"
	sc_src = scpath + "src"
	sc_lib = scpath + "lib/"

	# compile FastTree
	if "FastTree" not in os.listdir(sc_lib):
		os.system("gcc -DUSE_DOUBLE -O3 -finline-functions -funroll-loops -Wall -o " + sc_lib + "FastTree " + sc_lib + "FastTree.c -lm ")

	# replace paths in scripts
	src_files = os.listdir(sc_src)
	for sfile in src_files:
		sf = open(sc_src + "/" + sfile, 'r').read()
		sf = sf.replace('#SYNERCLUST_PATH', sc_bin)
		sf = sf.replace('#FASTTREE_PATH', sc_lib + "FastTree")
		if platform.system() == "Linux":
			sf = sf.replace('#MUSCLE_PATH', sc_lib + "muscle3.8.31_i86linux64")
		elif platform.system() == "Darwin":
			sf = sf.replace('#MUSCLE_PATH', sc_lib + "muscle3.8.31_i86darwin64")
		if len(args.blast) > 0 and args.blast[-1] != "/":
			sf = sf.replace('#BLAST_PATH', args.blast + "/")
		else:
			sf = sf.replace('#BLAST_PATH', args.blast)
		sf = sf.replace('#BLAST_EVAL_DEFAULT', args.evalue)
		sf = sf.replace('#NUM_CORES_DEFAULT', args.threads)

		mod_sf = sc_bin + sfile
		with open(mod_sf, 'w') as msf:
			msf.write(sf)

		os.system("chmod ugo+x " + mod_sf)

	print "SynerClust has been installed. It can be run from ~/bin/synerclust.py.  See the README for additional information."


if __name__ == "__main__":
	main()
