import sys
import argparse
import networkx as nx
from collections import defaultdict
from Graph import Graph, Node
import json

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
            edges.append((node1, node1dir, node2, node2dir, ovlp))

    G = Graph()
    nxg = nx.Graph()
    for node1, node1dir, node2, node2dir, ovlp in edges:
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
    parser.add_argument('json', type = str, help = 'Input JSON file')
    parser.add_argument('fa', type = str, help = 'Input FA file')
    args = parser.parse_args()
    jsonfile = args.json
    f = open(jsonfile, 'r')
    longreads = f.readlines()
    G, nxg = readGFA(args.gfa)
    #components = divide(nxg)
    print("Start computing all paths..", file=sys.stderr)
    paths = G.getAllPaths()
    print("Computed in total %d paths." % len(paths), file=sys.stderr)
    paths_between_tips = defaultdict(list)
    tips = []
    for path in paths:
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

    hap_count = 0

    out = open(args.fa, 'w')
#    for (tip1, tip2), paths in paths_between_tips.items():
#        paths_with_lengths = sorted([(G.getPathSeqLength(path), path) for path in paths], key=lambda entry: entry[0])
#        print('Found %d paths between %s and %s, longest sequence length %d' % (len(paths_with_lengths), tip1, tip2, paths_with_lengths[-1][0]), file=sys.stderr)
#        longest_path = paths_with_lengths[-1][1]
#        input_file_name_prefix = "_".join(args.gfa.split("/")[-1].split(".")[:-1])
#    hap_count = 0
    input_file_name_prefix = "_".join(args.gfa.split("/")[-1].split(".")[:-1])
    for i in range(len(tips)):
        for j in range(i+1,len(tips)):
            edges1 = list(nx.bfs_edges(nxg, tips[i], depth_limit=2))
            edges2 = list(nx.bfs_edges(nxg, tips[j], depth_limit=2))
            _, threeawaytip1 = edges1[1]
            _, threeawaytip2 = edges2[1]
            for k, line in enumerate(longreads):
                if threeawaytip1 in line and threeawaytip2 in line:
                    tips1 = min(tips[i], tips[j])
                    tips2 = max(tips[i], tips[j])
                    print('Found %d paths between %s and %s' % (len(paths_between_tips[(tips1,tips2)]), tips1, tips2), file=sys.stderr)
                    for path in paths_between_tips[(tips1,tips2)]:
                        hap_count += 1
                        out.write('>%s_hap%d\n' % (input_file_name_prefix, hap_count))
                        out.write(G.getPathSeq(path))
                        out.write('\n')
                    break
    print("Total paths found is: ", hap_count)
    out.close()

if __name__ == '__main__':
    main()
