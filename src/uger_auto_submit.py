#!/usr/bin/env python

import os
import subprocess
import argparse
import time
import shlex


def main():
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
	jobids = []

	while current < len(commands):
		checker = subprocess.Popen(cmd1, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=DEVNULL)
		output = checker.communicate()[0]

		# count number of slots currently being used
		njobs = 0
		lines = output.split("\n")[2:-1]
		for line in lines:
			tabs = line.split()
			njobs += int(tabs[-1])

		if njobs < args.limit:
			for i in xrange(njobs, args.limit, args.cores):
				# submit new job
				if "qsub" in commands[current]:
					cmd = shlex.split(commands[current].replace("${TIME}", timestamp).rstrip())
					submitted = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=DEVNULL, shell=True)
					output = submitted.communicate()[0]
					if output:
						jobids.append(output.split(" ")[2])
					# os.system(commands[current].replace("${TIME}", timestamp).rstrip())
					current += 1
					if current == len(commands):
						break
						# exit()
				else:
					current += 1

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
