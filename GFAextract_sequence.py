import sys
import argparse
import networkx as nx
from collections import defaultdict
from Graph import Graph, Node

'''
Usage: python this.py GFA regionLength(bp)
'''

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
            node1len = int(fields[6])
            node2len = int(fields[7])
            edges.append((node1, node1dir, node2, node2dir, ovlp, node1len, node2len))

    G = Graph()
    nxg = nx.Graph()
    for node1, node1dir, node2, node2dir, ovlp, node1len, node2len in edges:
        if node1 not in G.nodemap:
            n1_seq = nodesSeq[node1]
            assert(len(n1_seq) == node1len)
            n1 = Node(node1, node1len, n1_seq)
        else:
            n1 = G.nodemap[node1]
        if node2 not in G.nodemap:
            n2_seq = nodesSeq[node2]
            assert(len(n2_seq) == node2len)
            n2 = Node(node2, node1len, n2_seq)
        else:
            n2 = G.nodemap[node2]
        G.addEdge(n1, node1dir, n2, node2dir, ovlp)
        nxg.add_node(node1)
        nxg.add_node(node2)
        nxg.add_edge(node1, node2)
    return G, nxg

def get_connected_components(graph):
    connected_components = []
    nodes = list(graph.nodes())

    while len(nodes)!=0:
        start_node = nodes.pop()
        queue = [start_node] #FIFO
        visited = [start_node]
        while len(queue)!=0:
            start = queue[0]
            queue.remove(start)
            neighbours = list(graph.neighbors(start))
            #print(neighbours)
            for neighbour in neighbours:
                if neighbour not in visited:
                    queue.append(neighbour)
                    visited.append(neighbour)
                    nodes.remove(neighbour)
        connected_components.append(visited)
        
    return connected_components


def divide(graph):
    components = get_connected_components(graph)
    print('Found', len(components), 'connected components', file=sys.stderr)
    return components


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('gfa', type = str, help = 'Input GFA graph file.')
    args = parser.parse_args()

    G, nxg = readGFA(args.gfa)
    #components = divide(nxg)    
    print("Start computing paths..", file=sys.stderr)
    paths = G.getAllPaths()
    print("Computed %d paths." % len(paths), file=sys.stderr)
    paths_between_tips = defaultdict(list)
    for path in paths:
        first_node = path[0][0]
        last_node = path[-1][0]
        if first_node < last_node:
            paths_between_tips[(first_node, last_node)].append(path)
        else:
            paths_between_tips[(last_node, first_node)].append(path)

    hap_count = 0
    for (tip1, tip2), paths in paths_between_tips.items():
        paths_with_lengths = sorted([(G.getPathSeqLength(path), path) for path in paths], key=lambda entry: entry[0])
        print('Found %d paths between %s and %s, longest sequence length %d' % (len(paths_with_lengths), tip1, tip2, paths_with_lengths[-1][0]), file=sys.stderr)
        longest_path = paths_with_lengths[-1][1]
        hap_count += 1
        input_file_name_prefix = "_".join(args.gfa.split("/")[-1].split(".")[:-1])
        print('>%s_hap%d' % (input_file_name_prefix, hap_count))
        print(G.getPathSeq(longest_path))

if __name__ == '__main__':
    main()
