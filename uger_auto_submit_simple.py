#!/usr/bin/env python

import os
import subprocess
import argparse
import time


parser = argparse.ArgumentParser("usage : Use this script to automatically launch new SynerClust jobs on UGER from a file containing the list of jobs.\nYou need to have ran \"use UGER\" before using this script")
parser.add_argument("-l", dest="limit", type=int, default=800, help="Maximum number of slots that can be used (including those taken by other jobs). default = 90")
parser.add_argument("-f", dest="file", default="uger_list.txt", help="File containing the commands to run. default = \"uger_list.txt\"")
parser.add_argument("-t", dest="wait", type=int, default=300, help="Interval of time at which to check the number of available slots (in seconds). default = 300 (5 minutes)")
parser.add_argument("-p", dest="project", default="gscid", help="Project on which to submit for priority.")
parser.add_argument("-q", dest="queue", default="short", help="Queue to submit on.")
parser.add_argument("-err", dest="error", default="/dev/null", help="Error output.")
parser.add_argument("-log", dest="log", default="/dev/null", help="Standard output.")
parser.add_argument("-tmp", dest="tmp", required=True, help="Folder where to write temporary submission scripts.")

args = parser.parse_args()
DEVNULL = open(os.devnull)
timestamp = str(int(time.time()))

script_base = """
#! /bin/bash

#$ -cwd
#$ -q """ + args.queue + """
#$ -P """ + args.project + """

#$ -l m_mem_free=2g
#$ -e """ + args.error + """
#$ -o """ + args.log + """

source /broad/software/scripts/useuse
reuse -q Python-2.7


"""

if args.tmp[-1] != "/":
	args.tmp += "/"

with open(args.file) as f:
	commands = f.readlines()

new_commands = []
for i in xrange(len(commands)):
	with open(args.tmp + timestamp + str(i), "w") as f:
		f.write(script_base + commands[i])
	new_commands.append("qsub -N j" + timestamp + str(i) + " " + args.tmp + timestamp + str(i))

current = 0
njobs = args.limit
cmd1 = ["qstat"]

while current < len(new_commands):
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
		for i in xrange(njobs, args.limit):
			# submit new job
			os.system(new_commands[current].replace("#TIMESTAMP", timestamp).rstrip())
			current += 1
			if current == len(new_commands):
				exit()

	time.sleep(args.wait)
