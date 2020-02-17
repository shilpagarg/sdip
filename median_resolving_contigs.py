import sys
import numpy

f = open(sys.argv[1], 'r')
l = []
for line in f:
	fields = line.strip().split(":")
	if fields[0] == "Total":
		continue
	num_contigs = int(fields[0])
	count = int(fields[1].strip())
	l += [num_contigs] * count
print("List:", l)
print("Median:", numpy.median(l))
print("Mean:", numpy.mean(l))