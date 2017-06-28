#!/usr/bin/env python

import os
import subprocess
import argparse
import time


parser = argparse.ArgumentParser("usage : Use this script to automatically launch new SynerClust jobs on UGER from a file containing the list of jobs.\nYou need to have ran \"use UGER\" before using this script")
parser.add_argument("-l", dest="limit", type=int, default=900, help="Maximum number of slots that can be used (including those taken by other jobs). default = 90")
parser.add_argument("-f", dest="file", default="uger_list.txt", help="File containing the commands to run. default = \"uger_list.txt\"")
parser.add_argument("-t", dest="wait", type=int, default=300, help="Interval of time at which to check the number of available slots (in seconds). default = 300 (5 minutes)")
parser.add_argument("-n", dest="cores", type=int, default=4, help="Number of CPUs requested per job")

args = parser.parse_args()
DEVNULL = open(os.devnull)
timestamp = str(int(time.time()))

with open(args.file) as f:
	commands = f.readlines()

current = 0
njobs = args.limit
cmd1 = ["qstat"]

while current < len(commands):
	checker = subprocess.Popen(cmd1, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=DEVNULL)
	output = checker.communicate()[0]

	# count number of slots currently being used
	njobs = 0
	for line in output.split("\n")[2:-1]:
		tabs = line.split()
		if tabs[4] == "r" or tabs[4] == "dr":
			njobs += int(tabs[8])
		else:
			njobs += int(tabs[7])

	if njobs < args.limit:
		for i in xrange(njobs, args.limit, args.cores):
			# submit new job
			if "qsub" in commands[current]:
				os.system(commands[current].replace("${TIME}", timestamp).rstrip())
				current += 1
				if current == len(commands):
					exit()
			else:
				current += 1

	time.sleep(args.wait)
