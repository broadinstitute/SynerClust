#!/usr/bin/env python

import sys, os, time

def usage():
	print "WF_WaitForFile.py [file's dir] [file] [frequency - default=60s]"
	sys.exit(1)

def main(argv):
	fileDir = argv[0]
	file = argv[1]
	freq = 60
	if len(argv) > 2:
		freq = int(argv[2])
	
	file_found = 0
	while not file_found:
		dirFiles = os.listdir(fileDir)
		if file in dirFiles:
			file_found = 1
		else:
			time.sleep(freq)
	sys.exit(0)


if __name__ == "__main__":
	if len(sys.argv) == 1:
		usage()
	else:
		main(sys.argv[1:])