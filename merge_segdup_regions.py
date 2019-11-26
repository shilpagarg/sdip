#!/usr/bin/python

import sys
import argparse
import networkx as nx

def is_overlapping(c1, s1, e1, c2, s2, e2, min_ovlp = 0.5):
	if c1 == None or s1 == None or e1 == None or c2 == None or s2 == None or e2 == None:
		return False
	if c1 == c2:
		if s2 >= s1:
			if s2 < e1:
				overlap = min(e2, e1) - s2
			else:
				return False
		else:
			if e2 > s1:
				overlap = min(e2, e1) - s1
			else:
				return False
	else:
		return False
	if overlap / (e2-s2) > min_ovlp and overlap / (e1-s1) > min_ovlp:
		return True
	else:
		return False
	
def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('segdups', nargs='?', type=argparse.FileType('r'), default=sys.stdin, help = 'segdup track from UCSC (default: stdin)')
	parser.add_argument('--min_frac', type=float, default=0.98, help = 'minimum match fraction (default: 0.98)')
	parser.add_argument('--min_ovlp', type=float, default=0.5, help = 'minimum overlap (default: 0.5)')
	parser.add_argument('--context', type=int, default=15000, help = 'additional context at each side (default: 15000)')
	args = parser.parse_args()

	segdup_file = args.segdups
	G = nx.Graph()
	# Read segdups and add them as nodes in a graph
	for line in segdup_file:
		if line.startswith("#"):
			continue
		fields = line.strip().split()
		chrom = fields[1]
		start = int(fields[2])
		end = int(fields[3])
		chrom2 = fields[7]
		start2 = int(fields[8])
		end2 = int(fields[9])
		uid = int(fields[11])
		fracMatch = float(fields[26])

		# Connect two regions with an edge if similar enough
		if fracMatch > args.min_frac:
			G.add_edge((chrom, start, end), (chrom2, start2, end2))

	# Connect overlapping regions with an edge
	chroms = set([c for c, s, e in G.nodes()])
	for chrom in chroms:
		nodes = [(c, s, e) for c, s, e in G.nodes() if c == chrom]
		for i in range(len(nodes)):
			for j in range(i+1, len(nodes)):
				c1, s1, e1 = nodes[i]
				c2, s2, e2 = nodes[j]
				if is_overlapping(c1, s1, e1, c2, s2, e2, min_ovlp = args.min_ovlp):
					G.add_edge((c1, s1, e1), (c2, s2, e2))
	
	# Extract connected components from graph
	components = sorted(nx.connected_components(G), key=len, reverse=True)
	component_id = 0
	for component in components:
		nodes = sorted(list(component))
		if len(nodes) > 1:
			regions = []
			current_chrom = None
			current_start = None
			current_end = None
			for c, s, e in nodes:
				if c != current_chrom or s >= current_end:
					if current_chrom != None:
						regions.append((current_chrom, current_start, current_end))
					current_chrom = c
					current_start = s
					current_end = e
				current_end = max(current_end, e)
			if current_chrom != None:
				regions.append((current_chrom, current_start, current_end))
			component_id += 1
			for c, s, e in regions:
				print("%d\t%s\t%d\t%d" % (component_id, c, max(0, s - args.context), e + args.context))

if __name__ == '__main__':
	main()
