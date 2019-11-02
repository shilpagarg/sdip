#!/usr/bin/python

import sys
import argparse
import json
from Gfa import *
from collections import Counter


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
	paths = []
	S = [[start]]
	while len(S) > 0:
		path = S.pop()
		last = path[-1]
		for edge in graph.edges[last]:
			target = edge[0]
			if target == end:
				paths.append(path + [target])
			else:
				S.append(path + [target])
	return paths


def remove_bubble(graph, paths):
	sorted_paths = sorted([(len(p), p) for p in paths])
	longest_path = sorted_paths[-1][1]
	nodes_in_longest_path = [node for node, direction in longest_path]
	for path in sorted_paths[:-1]:
		for node, direction in path[1]:
			if not node in nodes_in_longest_path and node in graph.nodes:
				del graph.nodes[node]
	graph.remove_nonexistent_edges()


def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('gfa', type = str, help = 'input GFA graph file')
	parser.add_argument('json', type = str, help = 'input json alignment file')
	parser.add_argument('action', choices=['list', 'remove'])
	parser.add_argument('--min_cov', type = int, default = 0, help = 'minimum edge coverage (default: 0 = retain all edges)')
	parser.add_argument('--additional_json', type = str, default = '', help = 'input json alignment file')
	args = parser.parse_args()

	ingraph = args.gfa

	graph = Graph()
	graph.load(ingraph)
	graph.remove_nonexistent_edges()

	json_file = open(args.json, 'r')
	json_content = json_file.read()
	alignments = json.loads(json_content)
	paths = [aln["nodes"] for aln in alignments if len(aln["nodes"]) > 1]
	coverage_counter = Counter()

	for path in paths:
		for i in range(len(path) - 1):
			coverage_counter[(path[i], path[i+1])] += 1

	if args.additional_json != '':
		json_file2 = open(args.additional_json, 'r')
		json_content2 = json_file2.read()
		alignments2 = json.loads(json_content2)
		paths2 = [aln["nodes"] for aln in alignments2 if len(aln["nodes"]) > 1]
		for path in paths2:
			for i in range(len(path) - 1):
				coverage_counter[(path[i], path[i+1])] += 1

	if args.action == "list":
		print(coverage_counter)
	elif args.action == "remove":
		num_removed = 0
		num_kept = 0
		for node_from in graph.edges:
			targets_before = graph.edges[node_from]
			targets_after = set()
			for node_to in targets_before:
				coverage = coverage_counter[(node_from[0], node_to[0][0])] + coverage_counter[(node_to[0][0], node_from[0])]
				if coverage >= args.min_cov:
					print("Keep edge ({0}, {1}) with coverage {2}".format(node_from[0], node_to[0][0], coverage), file=sys.stderr)
					targets_after.add(node_to)
					num_kept += 1
				else:
					print("Remove edge ({0}, {1}) with coverage {2}".format(node_from[0], node_to[0][0], coverage), file=sys.stderr)
					num_removed += 1
			graph.edges[node_from] = targets_after
		print("Done. Removed {0} out of {1} edges.".format(num_removed, num_removed+num_kept), file=sys.stderr)
		graph.write(sys.stdout)

if __name__ == '__main__':
	main()
