import sys
from collections import defaultdict
import subprocess
from multiprocessing import Pool, Process, Manager
import os
t = 150

def getChunk(gamFilecache, threadN):
	NR = len(gamFilecache)
	if NR % threadN == 0:
		chunkSize = int(NR / threadN)
	else:
		if NR % threadN != 0:
			chunkSize = NR // threadN + 1
		else:
			chunkSize = NR / threadN
	chunks = []
	for i in range(threadN):
		if (i + 1) * chunkSize <= NR: 
			chunks.append((i*chunkSize, (i+1)*chunkSize))
		else:
			chunks.append((i*chunkSize, NR))
	return chunks

def matchLines(chunk, reads, results):
	ovlp = defaultdict(set)
	for line in chunk:
		tokens = line.split()
		if float(tokens[3]) >= 99.4:
			ovlp[tokens[0]].add(tokens[1])
			ovlp[tokens[1]].add(tokens[0])
	out = set()
	for r in reads:
		out = out.union(ovlp[r])
	results.extend(list(out))


dic = sys.argv[1]
ovlpPath = sys.argv[2]
reads = sys.argv[3]

sys.stderr.write('Reading id-readName table...\n')
D = open(dic, 'r').read().split('\n')[:-1]
id2read = {line.split()[0]:line.split()[1] for line in D}
read2id = {line.split()[1]:line.split()[0] for line in D}
sys.stderr.write('id-readName hashes built\n')
readsets = [read2id[i] for i in open(reads, 'r').read().split('\n')[:-1]]
ovlpFiiles = [ovlpPath+i for i in os.listdir(ovlpPath)]

def doOneChunk(ovlp):
	sys.stderr.write('Reading overlaps...\n')
	O = open(ovlp, 'r').read().split('\n')[:-1]
	chunckidx = getChunk(O, t)
	sys.stderr.write('Overlaps read, chunks built\n')

	

	sys.stderr.write('Start to run parallel...\n')
	m = Manager()
	results = m.list()

	p = Pool(24)

	processes = []
	for c in range(t):
		processes.append(p.apply_async(matchLines, (O[chunckidx[c][0]:chunckidx[c][1]], readsets, results)))

	p.close()
	p.join()

	for i in set(results):
		print(id2read[i])


for i in ovlpFiiles:
	doOneChunk(i)
