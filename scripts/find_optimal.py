#!/usr/bin/python

import os
import sys
import time
import ctypes

#importing ctypes time to have a higher resolution
CLOCK_MONOTONIC_RAW = 4

class timespec(ctypes.Structure):
	_fields_ = [
		('tv_sec', ctypes.c_long),
		('tv_nsec', ctypes.c_long)
	]

librt = ctypes.CDLL('librt.so.1', use_errno=True)
clock_gettime = librt.clock_gettime
clock_gettime.argtypes = [ctypes.c_int, ctypes.POINTER(timespec)]

def monotonic_time():
	t = timespec()
	if clock_gettime(CLOCK_MONOTONIC_RAW, ctypes.pointer(t)) != 0:
		errno_ = ctypes.get_errno()
		raise OSError (errno, os.strerror(errno_))
	return t.tv_sec + t.tv_nsec*1e-9

filename = os.environ.get ('PYTHONSTARTUP')
if filename and os.path.isfile(filename):
	execfile(filename)

if len (sys.argv) != 5:
	print "./find_optimal [path] [range start] [range stop] [output]"
	sys.exit (-1)

#this is our goofit directory to find optimal
path = sys.argv[1]
begin = int(sys.argv[2])
end = int(sys.argv[3])
output = sys.argv[4]

#this a list for our cvs file
strings = []

#we need to make the first string, which will be the header for the csv file
header = "block,"
for y in range (1, 32):
	header += "grain:{}, total,".format (y)

strings.append (header + "\r\n")

#copy our path to /tmp and work out of there
os.system ("cp -r {} /tmp".format(path))

#remove the path, only grab the very end...
wrkdir = "/tmp/{}".format("GooFit-SoA/")

#loop over the possible range
for i in range (begin, end):
        s = "{},".format(str(begin))
	for j in range (1, 32):
		os.chdir (wrkdir)
		#i is our group size
		#j is our grain size

		os.system ("make clean")
		os.system ("make OVERRIDE_GROUPSIZE={} OVERRIDE_GRAINSIZE={}".format(i, j))

		os.chdir ("examples/dalitz")

		os.system ("make clean && make OVERRIDE_GROUPSIZE={} OVERRIDE_GRAINSIZE={}".format(i, j))

		before = monotonic_time ()
		os.system("./dalitz dalitz_toyMC_000.txt > log")
		after = monotonic_time()

		#we are going to do two things: first parse the log file
		f = open ("log", 'r')
		lines = f.readlines ()
		f.close ()

		#we are pulling out 3rd from last line
		line = lines[len(lines) - 4]

		#Now we are going to take the difference between after and begin
		elapsed = after - before

		s += line[18:-10] + "," + str(elapsed) + ","

	strings.append (s + "\r\n")

timing = open(output + ".csv", 'w')

for s in strings:
	timing.write(s)

timing.close()

print 'Done'
