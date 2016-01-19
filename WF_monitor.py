#!/usr/bin/env python
import sys, os, time, string

def usage():
	print """WF_monitor.py [workflow_instance]
	Does this help?
	"""
	
def subflowStatus(instance):
	if len(instance) < 1:
		return ""
	myInstance = open(instance,'r').readlines()
	curSet = ""
	states = {}
	running = ""
	r_count = 0
	p_count = 0
	c_count = 0
	o_count = 0
	for i in myInstance:
		i = i.rstrip()
		if i.find("name") > -1:
			right = i.index(">")+1
			curSet= i[right:i.index("<",right)]
		elif i.find("Command set with name")>-1:
			continue
		elif i.find("state")>-1:
			right = i.index(">")+1
			myState= i[right:i.index("<",right)]
			if not myState in states:
				states[myState] = 0
			states[myState]+=1
			if myState=="running":
				if not curSet == running:
					r_count=0
				running = curSet
				r_count+=1
			elif myState=="complete":
				c_count+=1
			elif myState=="pending":
				p_count+=1
			else:
				o_count+=1
	t_count = o_count+c_count+r_count+p_count
	percent = str(float(c_count)/float(t_count))
	sys.stdout.write(percent[:5]+"\t"+str(r_count)+"\t"+str(p_count)+"\t"+str(o_count)+"\t")
	sstring = ""
	for s in states:
		sstring += s+": "+str(states[s])+"     "
	#~ print sstring
	if len(running) > 1:
		return (running,r_count)
	else:
		return ""
			
	

def main(argv):
	wf_instance =argv[0]
	refresh = 30
	if len(argv) > 1:
		refresh = int(argv[1])
	myLog = {}

	while(1):
		states = {}
		running = {}
		myInstance = open(wf_instance, 'r').readlines()
		curSet = ""
		curID = ""
		curStr = ""
		curFile = ""
		curDispID = ""
		subFlow = ""
		curArg = ""
		curNode = ""
		myState = ""
		waiting = []
		failed = []
		for i in myInstance:
			i = i.rstrip()
			if i.find("name") > -1:
				right = i.index(">")+1
				curSet= i[right:i.index("<",right)]
				curFile = ""
				if curSet.find("subflow") > -1:
					subFlow1 = curSet.split("flow")[-1]
					subFlow = subFlow1.split("&")[0]
			elif i.find("<id>")>-1:
				right = i.index(">")+1
				curID= i[right:i.index("<",right)]
				curDispID = subFlow+";"+curSet
				curStr = subFlow+";"+curSet+";"+curID
				if not curStr in myLog:
					myLog[curStr]=""
			elif i.find("fileName")>-1:
				right = i.index(">")+1
				curFile= i[right:i.index("<",right)]
			elif i.find("arg")>-1:
				right = i.index(">")+1
				temp= i[right:i.index("<",right)]
				if not (temp.find("NODE")>-1 or temp.find("CLUSTERS")>-1):
					curNode = temp.split("/")[-2]
					curArg = temp+curNode+"_instance.xml"
					if not os.path.isfile(curArg):
						curArg = ""
					elif myState == "running":
						myRunning = subflowStatus(curArg)
						curDispID = subFlow+";"+curNode
						print curDispID, myState
						if len(myRunning) >1:
							if not myRunning[0] in running:
								running[myRunning[0]] = 0
							running[myRunning[0]] +=myRunning[1]
			elif i.find("Command set with name")>-1:
				continue
			elif i.find("state")>-1:
				right = i.index(">")+1
				myState= i[right:i.index("<",right)]
				if not curStr.find("subflow")>-1:
					if not myState in states:
						states[myState] = 0
					states[myState]+=1
				if myState=="error" or myState=="failed":
					#~ print curStr, myState
					if curStr.find("subflow")>-1:
						failed.append(subFlow)
				if not myState == myLog[curStr]:
					stata = myLog[curStr]+"---->"+myState
					print curDispID, stata
					myLog[curStr] = myState
				elif myState == "running" or myState == "pending":
					myRunning = ""
					if curStr.find("subflow")>-1:
						continue
					elif curStr.find("WaitForNode")>-1:
						continue
						#~ curDispID = "\t\t\t\t"+curDispID
						#~ waiting.append(subFlow)
						#~ continue
					else:
						myRunning = subflowStatus(curFile)
					#~ print curDispID, myState
					if len(myRunning) >1:
						if not myRunning[0] in running:
							running[myRunning[0]] = 0
						running[myRunning[0]] +=myRunning[1]

		print "\nPending Jobs: "
		os.system("""bjobs| grep "PEND" |wc -l""")		
		print "Running Jobs: "
		os.system("""bjobs | grep "RUN" |wc -l""")		
		for s in states:
			if s=="failed":
				print s+":  "+str(states[s])+"\t"+" ".join(failed)
			else:
				print s+":  "+str(states[s])
		#~ if len(waiting)>0:	
			#~ print "waiting",len(waiting)
			#~ print " ".join(waiting)
		#~ for r in running:
			#~ print r+":   "+str(running[r])

		#~ print "incomplete:",states["incomplete"]
		#~ print "pending:",states["pending"]
		#~ print "running:",states["running"]
		#~ print "complete:",states["complete"]
		#~ print "failed:",states["failed"]
		print "\n"
		time.sleep(refresh)

				
if __name__ == "__main__":
	if len(sys.argv) == 1:
		sys.exit(usage())
	else:
		main(sys.argv[1:])