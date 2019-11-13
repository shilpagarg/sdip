import sys
from collections import defaultdict
import subprocess
from multiprocessing import Pool, Process, Manager
import os
import argparse

def getChunk(gamFilecache, num_chunks):
	NR = len(gamFilecache)
	if NR % num_chunks == 0:
		chunkSize = int(NR / num_chunks)
	else:
		if NR % num_chunks != 0:
			chunkSize = NR // num_chunks + 1
		else:
			chunkSize = NR / num_chunks
	chunks = []
	for i in range(num_chunks):
		if (i + 1) * chunkSize <= NR: 
			chunks.append((i*chunkSize, (i+1)*chunkSize))
		else:
			chunks.append((i*chunkSize, NR))
	return chunks

def doOneChunk(chunk, reads, results):
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

def doOneFile(ovlp, reads, num_chunks, num_threads):
	sys.stderr.write('Reading overlaps...\n')
	O = open(ovlp, 'r').read().split('\n')[:-1]
	chunckidx = getChunk(O, num_chunks)
	sys.stderr.write('Overlaps read, chunks built\n')

	sys.stderr.write('Start to run parallel...\n')
	m = Manager()
	results = m.list()

	p = Pool(num_threads)

	processes = []
	for c in range(num_chunks):
		processes.append(p.apply_async(doOneChunk, (O[chunckidx[c][0]:chunckidx[c][1]], reads, results)))

	p.close()
	p.join()

	return set(results).union(reads)

def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('dict', type = str, help = 'id-to-readname dictionary')
	parser.add_argument('ovlp_path', type = str, help = 'path to overlap files')
	parser.add_argument('reads', type = str, help = 'read list')
	parser.add_argument('--num_chunks', type = int, default=150, help = 'number of chunks (default: 150)')
	parser.add_argument('--num_threads', type = int, default=1, help = 'number of threads (default: 1)')
	parser.add_argument('--min_identity', type = float, default=99.4, help = 'minimum overlap identity (default: 99.4)')
	args = parser.parse_args()

	sys.stderr.write('Reading id-readName table...\n')

	dict_file = open(args.dict, 'r')
	id2read = {}
	read2id = {}
	for line in dict_file:
		fields = line.strip().split()
		id2read[fields[0]] = fields[1]
		read2id[fields[1]] = fields[0]

	sys.stderr.write('id-readName table built\n')

	reads = [read2id[i] for i in open(args.reads, 'r').read().split('\n')[:-1]]
	ovlpFiles = [args.ovlp_path+i for i in os.listdir(args.ovlp_path)]
	for i in ovlpFiles:
		results = doOneFile(i, reads, args.num_chunks, args.num_threads)
		for i in results:
			print(id2read[i])

if __name__ == '__main__':
	main()
