import sys
import networkx as nx
import argparse

def readGFA(gfaFile):
    gfa = open(gfaFile).read().split('\n')[:-1]
    edges= list()
    key_to_lemon_id = {}
    lemon_id_to_key = {}
    running_index = 1
    for line in gfa:
        if line[0] == 'L':
            fields = line.split('\t')
            node1 = fields[1]
            node1dir = fields[2]
            key1 = (node1, node1dir)
            if not key1 in key_to_lemon_id:
                key_to_lemon_id[key1] = running_index
                lemon_id_to_key[running_index] = key1
                running_index += 1
            id1 = key_to_lemon_id[key1]
            
            node2 = fields[3]
            node2dir = fields[4]
            key2 = (node2, node2dir)
            if not key2 in key_to_lemon_id:
                key_to_lemon_id[key2] = running_index
                lemon_id_to_key[running_index] = key2
                running_index += 1
            id2 = key_to_lemon_id[key2]

            ovlp = int(fields[5][:-1])
            edges.append((id1, id2))
    nodes_directed = set()
    for node1, node2 in edges:
        nodes_directed.add(node1)
        nodes_directed.add(node2)
    return nodes_directed, edges, key_to_lemon_id, lemon_id_to_key


def store_to_file(key_to_lemon_id, id_table):
    with open(id_table, 'w') as f:
        for key, lemon_id in key_to_lemon_id.items():
            nodename, nodedir = key
            print("\t".join([nodename, nodedir, str(lemon_id)]), file=f)

def checkComponents(components, lemon_id_to_key):
    """Check whether any component has both directions of the same read."""
    for component in components:
        seen = set()
        for node in component:
            nodename, nodedir = lemon_id_to_key[node]
            if nodename in seen:
                print("Error: This graph seems to contain cycles. See node %s." % (nodename), file=sys.stderr)
                return False
            seen.add(nodename)
    return True

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('gfa', type = str, help = 'graph in gfa format')
    parser.add_argument('id_table', type = str, help = 'output table for read to lemon id translation')
    args = parser.parse_args()
    
    nodes, edges, key_to_lemon_id, lemon_id_to_key = readGFA(args.gfa)
    store_to_file(key_to_lemon_id, args.id_table)
    
    G = nx.Graph()
    for node1, node2 in edges:
        G.add_edge(node1, node2)
    components = list(nx.connected_components(G))
    if not checkComponents(components, lemon_id_to_key):
        return False

    chosen_nodes = []
    for component in components:
        new_component = True
        for node in component:
            nodename, nodedir = lemon_id_to_key[node]
            if key_to_lemon_id[(nodename, '+')] in chosen_nodes or key_to_lemon_id[(nodename, '-')] in chosen_nodes:
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
