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
            lemon_id = zmw_hole_nr[-6:] + run_time[-2:]
            nodes.add(lemon_id)
        elif line[0] == 'L':
            fields = line.split('\t')
            node1 = fields[1]
            node1_parts = node1.split("/")
            run_time1 = node1_parts[0]
            zmw_hole_nr1 = node1_parts[1]
            node1dir = fields[2]
            id1_full = zmw_hole_nr1[-6:] + run_time1[-2:] + ("1" if node1dir == "+" else "0")

            node2 = fields[3]
            node2_parts = node2.split("/")
            run_time2 = node2_parts[0]
            zmw_hole_nr2 = node2_parts[1]
            node2dir = fields[4]
            id2_full = zmw_hole_nr2[-6:] + run_time2[-2:] + ("1" if node2dir == "+" else "0")

            ovlp = int(fields[5][:-1])
            edges.append((id1_full, id2_full))
    nodes_directed = set()
    for node1, node2 in edges:
        assert(node1[:-1] in nodes)
        assert(node2[:-1] in nodes)
        nodes_directed.add(node1)
        nodes_directed.add(node2)
    return nodes_directed, edges


def checkComponents(components):
    """Check whether any component has both directions of the same read."""
    for component in components:
        seen = set()
        for node in component:
            node_without_dir = node[:-1]
            if node_without_dir in seen:
                print("Error: This graph seems to contain cycles. See node %s." % (node_without_dir), file=sys.stderr)
                return False
            seen.add(node_without_dir)
    return True

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('gfa', type = str, help = 'graph in gfa format')
    args = parser.parse_args()
    
    nodes, edges = readGFA(args.gfa)
    G = nx.Graph()
    for node1, node2 in edges:
        G.add_edge(node1, node2)
    components = list(nx.connected_components(G))
    if not checkComponents(components):
        return False

    chosen_nodes = []
    for component in components:
        new_component = True
        for node in component:
            node_without_dir = node[:-1]
            if (node_without_dir + "1") in chosen_nodes or (node_without_dir + "0") in chosen_nodes:
                new_component = False
                break
        if new_component:
            chosen_nodes.extend(component)

    chosen_edges = [(n1, n2) for (n1, n2) in edges if (n1 in chosen_nodes and n2 in chosen_nodes)]

    print('@nodes')
    print('label')
    for i in chosen_nodes:
        print(i)

    print('@arcs')
    print('\t\tlabel\tweight')
    count=1
    for node1, node2 in chosen_edges:
        print(node1, node2, str(count), str(1))
        #print(edges[i].split(':')[0].split("/")[1], edges[i].split(':')[1].split("/")[1], str(count+1), str(1))
        count+=1

if __name__ == '__main__':
    main()
