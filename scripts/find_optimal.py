#!/usr/bin/env python

import ctypes
import os
import sys

# importing ctypes time to have a higher resolution
CLOCK_MONOTONIC_RAW = 4


class timespec(ctypes.Structure):
    _fields_ = [("tv_sec", ctypes.c_long), ("tv_nsec", ctypes.c_long)]


librt = ctypes.CDLL("librt.so.1", use_errno=True)
clock_gettime = librt.clock_gettime
clock_gettime.argtypes = [ctypes.c_int, ctypes.POINTER(timespec)]


def monotonic_time():
    t = timespec()
    if clock_gettime(CLOCK_MONOTONIC_RAW, ctypes.pointer(t)) != 0:
        errno_ = ctypes.get_errno()
        raise OSError(errno_, os.strerror(errno_))
    return t.tv_sec + t.tv_nsec * 1e-9


if len(sys.argv) != 6:
    print("./find_optimal [path] [range start] [range stop] [range increment] [output]")
    sys.exit(-1)

# this is our goofit directory to find optimal
path = sys.argv[1]
begin = int(sys.argv[2])
end = int(sys.argv[3])
inc = int(sys.argv[4])
output = sys.argv[5]

# this a list for our cvs file
strings = []

# we need to make the first string, which will be the header for the csv file
header = "block,"
for y in range(4, 11):
    header += f"grain:{y}, total,"

strings.append(header + "\r\n")

# copy our path to /tmp and work out of there
os.system(f"cp -r {path} /tmp")

# remove the path, only grab the very end...
wrkdir = "/tmp/GooFit/"

# go to this directory
os.chdir(wrkdir)

minGrain = 0
minGroup = 0
currentMin = 1000000

# loop over the possible range
for i in range(begin, end):
    s = f"{begin!s},"
    for j in range(4, 11):
        # i is our group size
        # j is our grain size

        os.system("mkdir build")
        os.chdir("build")

        os.system("make clean")
        os.system(
            f"cmake ../ -DGOOFIT_CUDA_OR_GROUPSIZE={i} -DGOOFIT_CUDA_OR_GRAINSIZE={j}"
        )
        os.system("make -j 12")

        os.chdir("examples/dalitz")

        before = monotonic_time()
        os.system("./dalitz > log")
        after = monotonic_time()

        # we are going to do two things: first parse the log file
        with open("log") as f:
            lines = f.readlines()

        # we are pulling out 3rd from last line
        line = lines[len(lines) - 4]

        # Now we are going to take the difference between after and begin
        elapsed = after - before

        if elapsed < currentMin:
            minGroup = i
            minGrain = j

        s += line[18:-10] + "," + str(elapsed) + ","

        os.chdir(wrkdir)
        os.system("rm -rf build")

    strings.append(s + "\r\n")

with open(f"{output}.csv", "w") as timing:
    for s in strings:
        timing.write(s)

print(f"Group: {minGroup} Grain: {minGrain}\n")

print("Done")
