#!/usr/bin/env python

# By Fabio Pakk

import os
import sys

def usage():
	print "	This script runs a command (arg 1) n times (arg 3)\n\
	storing the outputs of each run in a a single file,\n\
	whose prefix is given by arg 2. The script will append\n\
	the number of the run to each output file plus a .txt\n\
	extension."
	sys.exit(1)

if len(sys.argv) != 4:
	usage()

for i in range(int(sys.argv[3])):
	command = sys.argv[1] + " > " + sys.argv[2] + str(i) + ".txt"
	print "Executing: \"" + command + "\""
	os.system(command)
