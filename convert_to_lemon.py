import sys
import networkx as nx
import argparse

def readGFA(gfaFile):
    gfa = open(gfaFile).read().split('\n')[:-1]
    nodes = set()
    edges= list()
    for line in gfa:
        if line[0] == 'S':
            fields = line.split('\t')
            name = fields[1]
            parts = name.split("/")
            run_time = parts[0]
            zmw_hole_nr = parts[1]
            lemon_id = zmw_hole_nr[-4:] + run_time[-4:]
            nodes.add(lemon_id)
        elif line[0] == 'L':
            fields = line.split('\t')
            node1 = fields[1]
            node1_parts = node1.split("/")
            run_time1 = node1_parts[0]
            zmw_hole_nr1 = node1_parts[1]
            node1dir = fields[2]
            id1_full = zmw_hole_nr1[-4:] + run_time1[-4:] + ("1" if node1dir == "+" else "0")

            node2 = fields[3]
            node2_parts = node2.split("/")
            run_time2 = node2_parts[0]
            zmw_hole_nr2 = node2_parts[1]
            node2dir = fields[4]
            id2_full = zmw_hole_nr2[-4:] + run_time2[-4:] + ("1" if node2dir == "+" else "0")

            ovlp = int(fields[5][:-1])
            edges.append((id1_full, id2_full))
    nodes_directed = set()
    for node1, node2 in edges:
        assert(node1[:-1] in nodes)
        assert(node2[:-1] in nodes)
        nodes_directed.add(node1)
        nodes_directed.add(node2)
    return nodes_directed, edges


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('gfa', type = str, help = 'graph in gfa format')
    args = parser.parse_args()
    
    nodes, edges = readGFA(args.gfa)
    G = nx.Graph()
    for node1, node2 in edges:
        G.add_edge(node1, node2)
    if nx.number_connected_components(G) != 2:
        print("Error: Number of connected compontents is %d and not 2." % (nx.number_connected_components(G)), file=sys.stderr)
        return False

    nodes_in_first_component = list(nx.connected_components(G))[0]
    edges_in_first_component = [(n1, n2) for (n1, n2) in edges if (n1 in nodes_in_first_component and n2 in nodes_in_first_component)]

    print('@nodes')
    print('label')
    for i in nodes_in_first_component:
        print(i)

    print('@arcs')
    print('\t\tlabel\tweight')
    count=1
    for node1, node2 in edges_in_first_component:
        print(node1, node2, str(count), str(1))
        #print(edges[i].split(':')[0].split("/")[1], edges[i].split(':')[1].split("/")[1], str(count+1), str(1))
        count+=1

if __name__ == '__main__':
    main()
