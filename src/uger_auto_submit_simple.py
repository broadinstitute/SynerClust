#!/usr/bin/env python

import os
import subprocess
import argparse
import shlex
import time


def main():
	parser = argparse.ArgumentParser("usage : Use this script to automatically launch new SynerClust jobs on UGER from a file containing the list of jobs.\nYou need to have ran \"use UGER\" before using this script")
	parser.add_argument("-l", dest="limit", type=int, default=800, help="Maximum number of slots that can be used (including those taken by other jobs). default = 90")
	parser.add_argument("-f", dest="file", default="uger_list.txt", help="File containing the commands to run. default = \"uger_list.txt\"")
	parser.add_argument("-t", dest="wait", type=int, default=300, help="Interval of time at which to check the number of available slots (in seconds). default = 300 (5 minutes)")
	parser.add_argument("-p", dest="project", default=None, help="Project on which to submit for priority.")
	parser.add_argument("-q", dest="queue", default="short", help="Queue to submit on.")
	parser.add_argument("-err", dest="error", default="/dev/null", help="Error output.")
	parser.add_argument("-log", dest="log", default="/dev/null", help="Standard output.")
	parser.add_argument("-tmp", dest="tmp", required=True, help="Folder where to write temporary submission scripts.")
	parser.add_argument("-h", dest="hour", type=str, default="48", help="Hours for runtime limit.")
	parser.add_argument("-m", dest="minute", type=str, default="00", help="Minutes for runtime limit.")
	parser.add_argument("-s", dest="second", type=str, default="00", help="Seconds for runtime limit.")

	args = parser.parse_args()
	DEVNULL = open(os.devnull)
	timestamp = str(int(time.time()))

	if args.error != "/dev/null":
		args.error += "/#NODE.err"

	if args.log != "/dev/null":
		args.log += "/#NODE.log"

	script_base = """
#! /bin/bash

#$ -cwd
#$ -q """ + args.queue + ("""
#$ -P """ + args.project if args.project else "") + """

#$ -l h_rt=""" + args.hour + ":" + args.hour + ":" + args.second + """
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
		if commands[i][0] != "#" and commands[i]:
			with open(args.tmp + timestamp + str(i), "w") as f:
				f.write(script_base.replace("#NODE", timestamp + str(i)) + commands[i])
			new_commands.append("qsub -N j" + timestamp + str(i) + " " + args.tmp + timestamp + str(i))

	current = 0
	njobs = args.limit
	cmd1 = ["qstat"]
	jobids = []

	while current < len(new_commands):
		checker = subprocess.Popen(cmd1, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=DEVNULL)
		output = checker.communicate()[0]

		# count number of slots currently being used
		njobs = 0
		lines = output.split("\n")[2:-1]
		for line in lines:
			tabs = line.split()
			njobs += int(tabs[-1])

		if njobs < args.limit:
			for i in xrange(njobs, args.limit):
				# submit new job
				cmd = shlex.split(new_commands[current].replace("${TIME}", timestamp).rstrip())
				submitted = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=DEVNULL)
				output = submitted.communicate()[0]
				jobids.append(output.split(" ")[2])
				# os.system(new_commands[current].replace("#TIMESTAMP", timestamp).rstrip())
				current += 1
				if current == len(new_commands):
					break
					# exit()

		time.sleep(args.wait)

	not_completed = True
	while not_completed:
		not_completed = False
		checker = subprocess.Popen(cmd1, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=DEVNULL)
		output = checker.communicate()[0]
		lines = output.split("\n")[2:-1]
		for line in lines:
			jobid = line.strip().split()[0]
			if jobid in jobids:
				not_completed = True
				time.sleep(args.wait)
				break

	print 'All jobs completed'
	exit(0)


if __name__ == "__main__":
		main()
