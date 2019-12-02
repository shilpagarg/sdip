#!/usr/bin/python

import sys
import argparse
from Gfa import *

def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('gfa', nargs='?', type=argparse.FileType('r'), default=sys.stdin, help = 'input GFA graph file (default: stdin)')
	#parser.add_argument('--max_size', type = int, default = -1, help = 'maximum bubble size (number of nodes between start and end), default: -1 = all bubbles')
	args = parser.parse_args()

	ingraph = args.gfa

	graph = Graph()
	graph.load(args.gfa)
	graph.remove_nonexistent_edges()

	edge_removed = True
	print("Start pruning paths between high-degree nodes..", file=sys.stderr)
	while edge_removed:
		edge_removed = False
		high_degree_nodes = []
		for n in graph.nodes:
			targets = graph.edges[(n, True)].union(graph.edges[(n, False)])
			if len(targets) > 2:
				high_degree_nodes.append(n)

		for start in high_degree_nodes:
			if edge_removed:
				break
            #print("Start with node " + start, file=sys.stderr)
			stack = []
			if len(graph.edges[(start, True)]) > 1:
				stack.append([(start, True)])
			if len(graph.edges[(start, False)]) > 1:
				stack.append([(start, False)])
			
			while len(stack) > 0:
				if edge_removed:
					break
				path = stack.pop()
				if len(path) > 1:
					continue
				last_node, last_orientation = path[-1]
				#print("At node " + last_node, file=sys.stderr)
				if last_orientation == True:
					for (nxt, nxt_dir), infos in graph.edges[(last_node, True)]:
						if nxt in high_degree_nodes:
							if len(graph.edges[(nxt, True if nxt_dir == False else False)]) > 1:
								if len(path) == 1:
									graph.edges[(last_node, last_orientation)].remove(((nxt, nxt_dir), infos))
									graph.edges[(nxt, not nxt_dir)].remove(((last_node, not last_orientation), infos))
                                    #print("Remove edge between %s,%s" % (last_node, nxt), file=sys.stderr)
									edge_removed = True
									break
						else:
							stack.append(path + [(nxt, nxt_dir)])
				else:
					for (nxt, nxt_dir), infos in graph.edges[(last_node, False)]:
						if nxt in high_degree_nodes:
							if len(graph.edges[(nxt, True if nxt_dir == False else False)]) > 1:
								if len(path) == 1:
									graph.edges[(last_node, last_orientation)].remove(((nxt, nxt_dir), infos))
									graph.edges[(nxt, not nxt_dir)].remove(((last_node, not last_orientation), infos))
                                    #print("Remove edge between %s,%s" % (last_node, nxt), file=sys.stderr)
									edge_removed = True
									break
						else:
							stack.append(path + [(nxt, nxt_dir)])
					
	print("Done.", file=sys.stderr)
	graph.write(sys.stdout)

if __name__ == '__main__':
	main()
