import sys

def readGFA(gfaFile):
    gfa = open(gfaFile).read().split('\n')[:-1]
    nodes = set()
    edges= list()
    for line in gfa:
        if line[0] == 'S':
            fields = line.split('\t')
            name = fields[1]
            id = name.split("/")[1]
            nodes.add(id)
        elif line[0] == 'L':
            fields = line.split('\t')
            node1 = fields[1]
            id1 = node1.split("/")[1]
            node1dir = fields[2]
            id1_full = id1 + "1" if node1dir == "+" else id1 + "0"
            node2 = fields[3]
            id2 = node2.split("/")[1]
            node2dir = fields[4]
            id2_full = id2 + "1" if node2dir == "+" else id2 + "0"
            ovlp = int(fields[5][:-1])
            edges.append((id1_full, id2_full))
    nodes_directed = set()
    for node1, node2 in edges:
        assert(node1[:-1] in nodes)
        assert(node2[:-1] in nodes)
        nodes_directed.add(node1)
        nodes_directed.add(node2)
    return nodes_directed, edges

nodes, edges = readGFA(sys.argv[1])
print('@nodes')
print('label')
for i in nodes:
    print(i)

print('@arcs')
print('\t\tlabel\tweight')
count=1
for node1, node2 in edges:
    print(node1, node2, str(count), str(1))
    #print(edges[i].split(':')[0].split("/")[1], edges[i].split(':')[1].split("/")[1], str(count+1), str(1))
    count+=1
