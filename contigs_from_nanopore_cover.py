import sys
import argparse
from collections import defaultdict, Counter
import json
import networkx as nx
from Graph import Graph, Node
from copy import deepcopy

def readGFA(gfaFile):
	gfa = open(gfaFile).read().split('\n')[:-1]
	
	nodesSeq = dict()
	edges = list()
	for line in gfa:
		if line[0] == 'S':
			fields = line.split('\t')
			nodesSeq[fields[1]] = fields[2]
		elif line[0] == 'L':
			fields = line.split('\t')
			node1 = fields[1]
			node1dir = fields[2]
			node2 = fields[3]
			node2dir = fields[4]
			ovlp = int(fields[5][:-1])
			cigar = fields[8]
			edges.append((node1, node1dir, node2, node2dir, ovlp, cigar))

	G = Graph()
	nxg = nx.Graph()
	for node1, node1dir, node2, node2dir, ovlp, cigar in edges:
		if node1 not in G.nodemap:
			n1_seq = nodesSeq[node1]
			n1 = Node(node1, len(n1_seq), n1_seq)
		else:
			n1 = G.nodemap[node1]
		if node2 not in G.nodemap:
			n2_seq = nodesSeq[node2]
			n2 = Node(node2, len(n2_seq), n2_seq)
		else:
			n2 = G.nodemap[node2]
		G.addEdge(n1, node1dir, n2, node2dir, ovlp, cigar)
		nxg.add_node(node1)
		nxg.add_node(node2)
		nxg.add_edge(node1, node2)
	return G, nxg

def lookup_nodename(lemon_id, nodemap):
	matching_nodes = [key for key in nodemap.keys() if key.split("/")[1] == lemon_id]
	assert(len(matching_nodes) == 1)
	return matching_nodes[0]

def remove_uncovered_edges(G, cover_dict):
	removed_edges = []
	for edge, cover in cover_dict.items():
		if cover == 0:
			(name1, dir1), (name2, dir2) = edge
			if dir1 == '+':
				G.nodemap[name1].Enodes.discard(name2)
			else:
				G.nodemap[name1].Bnodes.discard(name2)
			if dir2 == '+':
				G.nodemap[name2].Bnodes.discard(name1)
			else:
				G.nodemap[name2].Enodes.discard(name1)
			del G.edgeOvlp[(name1, name2)]
			removed_edges.append(edge)
			#print("Remove edge between", name1, "+" if dir1 else "-", "and", name2, "-" if dir2 else "-", file=sys.stderr)
	for edge in removed_edges:
		del cover_dict[edge]

def getPathsFromPathCover(G, cover_dict):
	remove_uncovered_edges(G, cover_dict)
	paths = []
	while len(G.nodemap) > 0:
		startOrEndNodes = list(G.getStartOrEndNodes())
		#print("Current start/end nodes:", startOrEndNodes, file=sys.stderr)
		first_node = startOrEndNodes[0]
		if G.nodemap[first_node].Enodes == set() and G.nodemap[first_node].Bnodes != set():
			current_path = [(first_node, '-')]
		elif G.nodemap[first_node].Bnodes == set() and G.nodemap[first_node].Enodes != set():
			current_path = [(first_node, '+')]
		else:
			print("Error: Found lonely node", file=sys.stderr)
			break
		while True:
			last_node, last_dir = current_path[-1]
			if last_dir == '+':
				#print("Forward: ", G.nodemap[last_node].Enodes)
				possible_nxt_nodes = []
				for nxt_node in G.nodemap[last_node].Enodes:
					if last_node in G.nodemap[nxt_node].Bnodes:
						nxt_dir = '+'
					elif last_node in G.nodemap[nxt_node].Enodes:
						nxt_dir = '-'
					else:
						continue
					if cover_dict[((last_node, last_dir), (nxt_node, nxt_dir))] > 0:
						possible_nxt_nodes.append((nxt_node, nxt_dir))
				if len(possible_nxt_nodes) > 0:
					nxt_node, nxt_dir = possible_nxt_nodes[0]
					cover_dict[((last_node, last_dir), (nxt_node, nxt_dir))] -= 1
					current_path.append((nxt_node, nxt_dir))
				else:
					#print("Reached end: ", ",".join([node + ("+" if direction else "-") for node, direction in current_path]), file=sys.stderr)
					paths.append(current_path)
					break
			else:
				#print("Backward: ", G.nodemap[last_node].Bnodes)
				possible_nxt_nodes = []
				for nxt_node in G.nodemap[last_node].Bnodes:
					if last_node in G.nodemap[nxt_node].Bnodes:
						nxt_dir = '+'
					elif last_node in G.nodemap[nxt_node].Enodes:
						nxt_dir = '-'
					else:
						continue
					if cover_dict[((last_node, last_dir), (nxt_node, nxt_dir))] > 0:
						possible_nxt_nodes.append((nxt_node, nxt_dir))
				if len(possible_nxt_nodes) > 0:
					nxt_node, nxt_dir = possible_nxt_nodes[0]
					cover_dict[((last_node, last_dir), (nxt_node, nxt_dir))] -= 1
					current_path.append((nxt_node, nxt_dir))
				else:
					#print("Reached end: ", ",".join([node + ("+" if direction else "-") for node, direction in current_path]), file=sys.stderr)
					paths.append(current_path)
					break
		remove_uncovered_edges(G, cover_dict)
		G.removeLonelyNodes()
	return paths

def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('gfa', type = str, help = 'graph in gfa format')
	parser.add_argument('json', type = str, help = 'nanopore alignments in json format')
	parser.add_argument('lemon', type = str, help = 'graph in lemon format')
	parser.add_argument('cover', type = str, help = 'mc-mpc solver output')
	parser.add_argument('--print_paths', action = "store_true", help = 'print paths (i.e. list of node names) on stderr')
	args = parser.parse_args()

	sys.stderr.write('Reading gfa graph...\n')
	G, nxg = readGFA(args.gfa)
	G_copy = deepcopy(G)
	G.removeLonelyNodes()
	sys.stderr.write('Read gfa graph.\n')

	json_file = open(args.json, 'r')
	longreads = json_file.readlines()

	sys.stderr.write('Reading lemon graph...\n')
	lemon_content = open(args.lemon, 'r').read().split('\n')[:-1]
	arcs = dict()
	arcs_section_started = False
	line_index = 0
	while line_index < len(lemon_content):
		if arcs_section_started:
			fields = lemon_content[line_index].strip().split()
			arcs[fields[2]] = (fields[0], fields[1])
		else:
			if lemon_content[line_index].startswith("@arcs"):
				arcs_section_started = True
				line_index += 2
				continue
		line_index += 1
	sys.stderr.write('Lemon graph read.\n')

	sys.stderr.write('Reading cover...\n')
	cover_file = open(args.cover, 'r')
	cover_dict = Counter()
	for line in cover_file:
		fields = line.strip().split()
		arc = fields[0]
		cover = int(fields[1])
		id1_full, id2_full = arcs[arc]
		id1 = id1_full[:-1]
		dir1 = '+' if id1_full[-1] == "1" else '-'
		id2 = id2_full[:-1]
		dir2 = '+' if id2_full[-1] == "1" else '-'
		name1 = lookup_nodename(id1, G.nodemap)
		name2 = lookup_nodename(id2, G.nodemap)
		cover_dict[((name1, dir1), (name2, dir2))] = cover
	sys.stderr.write('Cover read.\n')

	print("Start computing all paths..", file=sys.stderr)
	all_paths = G.getAllPaths()
	print("Computed in total %d paths." % len(all_paths), file=sys.stderr)
	paths_between_tips = defaultdict(list)
	tips = []
	for path in all_paths:
		first_node = path[0][0]
		last_node = path[-1][0]
		if first_node < last_node:
			paths_between_tips[(first_node, last_node)].append(path)
		else:
			paths_between_tips[(last_node, first_node)].append(path)
		if first_node not in tips:
			tips += [first_node]
		if last_node not in tips:
			tips += [last_node]
	
	print("Start extracting paths from nanopore alignments..", file=sys.stderr)
	hap_count = 0
	cover_dict_copy = deepcopy(cover_dict)
	for tip_i in range(len(tips)):
		for tip_j in range(tip_i+1,len(tips)):
			edges1 = list(nx.bfs_edges(nxg, tips[tip_i], depth_limit=2))
			edges2 = list(nx.bfs_edges(nxg, tips[tip_j], depth_limit=2))
			_, threeawaytip1 = edges1[1]
			_, threeawaytip2 = edges2[1]
			for line in longreads:
				if threeawaytip1 in line and threeawaytip2 in line:
					tip1 = min(tips[tip_i], tips[tip_j])
					tip2 = max(tips[tip_i], tips[tip_j])
					print('Found %d paths between %s and %s. Choosing one randomly.' % (len(paths_between_tips[(tip1,tip2)]), tip1, tip2), file=sys.stderr)
					chosen_path = paths_between_tips[(tip1,tip2)][0]
					for path_i in range(len(chosen_path) - 1):
						from_node, from_dir = chosen_path[path_i]
						to_node, to_dir = chosen_path[path_i+1]
						assert cover_dict[((from_node, from_dir), (to_node, to_dir))] >= 1
						cover_dict[((from_node, from_dir), (to_node, to_dir))] -= 1
					hap_count += 1
					path_seq = G_copy.getPathSeq(chosen_path)
					print('>hap%d len=%d' % (hap_count, len(path_seq)))
					print(path_seq)
					if args.print_paths:
						print("Path %d: " % (hap_count) + ",".join([node for node, direction in chosen_path]), file=sys.stderr)
					break
	remove_uncovered_edges(G, cover_dict)
	G.removeLonelyNodes()
	print("Done extracting paths from nanopore alignments.", file=sys.stderr)
	print(cover_dict_copy - cover_dict, file=sys.stderr)

	print("Start extracting paths from path cover..", file=sys.stderr)
	paths = getPathsFromPathCover(G, cover_dict)
	print('Found %d paths from path cover' % (len(paths)), file=sys.stderr)
	
	for path in paths:
		hap_count += 1
		path_seq = G_copy.getPathSeq(path)
		print('>hap%d len=%d' % (hap_count, len(path_seq)))
		print(path_seq)
		if args.print_paths:
			print("Path %d: " % (hap_count) + ",".join([node for node, direction in path]), file=sys.stderr)

if __name__ == '__main__':
	main()
