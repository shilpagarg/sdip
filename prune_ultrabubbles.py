#!/usr/bin/python

# Detecting Superbubbles in Assembly Graphs, Onodera et al 2013

import sys
import argparse
from Gfa import *

def getset(s):
	assert(len(s) == 1)
	for item in s:
		return item

# Fig. 5
def detect_bubble(graph, s):
	if s not in graph.edges: return None
	if len(graph.edges[s]) < 2: return None
	S = [s]
	visited = set()
	seen = set()
	seen.add(s)
	while len(S) > 0:
		v = S.pop()
		assert v in seen
		seen.remove(v)
		assert v not in visited
		visited.add(v)
		if v not in graph.edges: return None
		if len(graph.edges[v]) == 0: return None
		for edge in graph.edges[v]:
			u = edge[0]
			if u[0] == v[0]: return None
			if reverse(u) in visited: return None
			if u == s: return None
			assert u not in visited
			seen.add(u)
			assert reverse(u) in graph.edges
			assert len(graph.edges[reverse(u)]) >= 1
			has_nonvisited_parent = False
			for parent_edge in graph.edges[reverse(u)]:
				parent = reverse(parent_edge[0])
				if parent not in visited: has_nonvisited_parent = True
			if not has_nonvisited_parent: S.append(u)
		if len(S) == 1 and len(seen) == 1 and S[0] == getset(seen):
			t = S.pop()
			if t in graph.edges:
				for edge in graph.edges[t]:
					if edge[0] == s: return None
			return (s, t, len(visited) - 1)
	return None

def get_paths_in_bubble(graph, start, end):
	if not (start[0] in graph.nodes and end[0] in graph.nodes):
		return None
	paths = []
	S = [[start]]
	while len(S) > 0:
		path = S.pop()
		last = path[-1]
		for edge in graph.edges[last]:
			target = edge[0]
			if target == end:
				paths.append(path + [target])
				if len(paths) > 1000:
					return paths
			else:
				S.append(path + [target])
	return paths


def remove_bubble(graph, paths):
	if paths == None:
		return
	sorted_paths = sorted([(len(p), p) for p in paths])
	longest_path = sorted_paths[-1][1]
	nodes_in_longest_path = [node for node, direction in longest_path]
	for path in sorted_paths[:-1]:
		for node, direction in path[1]:
			if not node in nodes_in_longest_path and node in graph.nodes:
				del graph.nodes[node]
	graph.remove_nonexistent_edges()


def get_bubble_length(graph, bubble):
	paths = get_paths_in_bubble(graph, bubble[0], bubble[1])
	if paths == None:
		return None
	max_length = None
	for path in paths:
		length = 0
		for index in range(1, len(path) - 1):
			from_node = path[index-1]
			to_node = path[index]
			overlaps = [overlap for (topos, (overlap, tags)) in graph.edges[from_node] if topos == to_node]
			assert len(overlaps) == 1
			overlap = overlaps[0]
			to_node_length = len(graph.nodes[to_node[0]].nodeseq)
			#print(to_node_length, overlap)
			length += len(graph.nodes[to_node[0]].nodeseq) - overlap
		overlaps = [overlap for (topos, (overlap, tags)) in graph.edges[path[-2]] if topos == path[-1]]
		assert len(overlaps) == 1
		overlap = overlaps[0]
		length -= overlap
		if max_length == None or length > max_length:
			max_length = length
	return max_length


def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('gfa', nargs='?', type=argparse.FileType('r'), default=sys.stdin, help = 'input GFA graph file (default: stdin)')
	parser.add_argument('action', choices=['list', 'remove'])
	parser.add_argument('--min_length', type = int, default = -1, help = 'filter out bubbles shorter than this length (length of longest path between start and end), default: -1 = filter all')
	args = parser.parse_args()

	ingraph = args.gfa

	graph = Graph()
	graph.load(args.gfa)
	graph.remove_nonexistent_edges()

	bubbles = []
	print("Start finding bubbles..", file=sys.stderr)
	for n in graph.nodes:
		bubble = detect_bubble(graph, (n, True))
		if bubble:
			bubbles.append(bubble)
		bubble = detect_bubble(graph, (n, False))
		if bubble:
			bubbles.append(bubble)
	print("Found %d bubbles of length < %d" % (len([bubble for bubble in bubbles if get_bubble_length(graph, bubble) < args.min_length or args.min_length == -1]), args.min_length), file=sys.stderr)

	if args.action == "list":
		for bubble in bubbles:
			max_length = get_bubble_length(graph, bubble)
			print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}".format(bubble[0][0], "+" if bubble[0][1] else "-", bubble[1][0], "+" if bubble[1][1] else "-", bubble[2], max_length))
	elif args.action == "remove":
		for bubble in bubbles:
			length = get_bubble_length(graph, bubble)
			if args.min_length == -1 or (length != None and length < args.min_length):
				paths = get_paths_in_bubble(graph, bubble[0], bubble[1])
				remove_bubble(graph, paths)
				print("Removed", bubble)
		graph.write(sys.stdout)

if __name__ == '__main__':
	main()
