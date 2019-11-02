import sys

def rev_comp(sequence):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    reverse_complement = "".join(complement.get(base, base) for base in reversed(sequence))
    return reverse_complement

class Graph(object):
    """docstring for Edge"""
    def __init__(self, ):
        self.nodemap = {}
        self.edgeOvlp = {}

    def addEdge(self, node1, node1dir, node2, node2dir, ovlp):
        self.nodemap[node1.name] = node1
        self.nodemap[node2.name] = node2

        if node1dir == '+':
            node1.Enodes.add(node2.name)
        else:
            node1.Bnodes.add(node2.name)
        if node2dir == '+':
            node2.Bnodes.add(node1.name)
        else:
            node2.Enodes.add(node1.name)

        self.edgeOvlp[(node1.name, node2.name)] = ovlp
        self.edgeOvlp[(node2.name, node1.name)] = ovlp

    def getStartOrEndNodes(self):
        startOrEnd = set()
        for node, N in self.nodemap.items():
            if N.Bnodes == set() or N.Enodes == set():
                startOrEnd.add(node)
        return startOrEnd

    def getRandomLongPath(self, read):
        path = []
        startOrEnd = self.getStartOrEndNodes()
        if read != None:
            last = read
        else:
            last = list(startOrEnd)[0]
        if self.nodemap[last].Enodes == set():
            now = list(self.nodemap[last].Bnodes)[0]
        else:
            now = list(self.nodemap[last].Enodes)[0]
        path.append(last)
        path.append(now)
        while now not in startOrEnd:
            if last in self.nodemap[now].Enodes:
                nx = list(self.nodemap[now].Bnodes)[0]
            else:
                nx = list(self.nodemap[now].Enodes)[0]
            if nx in path:
                print('Cycling!!!!!! Break it at', now)
                break
            path.append(nx)
            last = now
            now = nx
        return path

    def getAllPaths(self):
        startOrEnd = self.getStartOrEndNodes()
        #print(len(startOrEnd), file=sys.stderr)
        paths = []
        for tip in startOrEnd:
            #print("Start with tip " + tip, file=sys.stderr)
            last = tip
            stack = []
            if self.nodemap[last].Enodes == set():
                stack.append([(last, '-')])
            else:
                stack.append([(last, '+')])            
            while len(stack) > 0:
                path = stack.pop()
                last_node, last_orientation = path[-1]
                #print("At node " + last_node, file=sys.stderr)
                if last_orientation == '+':
                    if len(self.nodemap[last_node].Enodes) > 0:
                        for nxt in self.nodemap[last_node].Enodes:
                            if (nxt, '+') in path or (nxt, '-') in path:
                                #print("Have seen this node before: " + nxt, file=sys.stderr)
                                continue
                            if last_node in self.nodemap[nxt].Bnodes:
                                stack.append(path + [(nxt, '+')])
                                #print("step forward", file=sys.stderr)
                            else:
                                stack.append(path + [(nxt, '-')])
                                #print("step forward", file=sys.stderr)
                    else:
                        #print("Reached tip " + str(len(paths)), file=sys.stderr)
                        if len(paths) < 10000:
                            paths.append(path)
                        else:
                            return paths
                        #print([(step[0][-7:-4], step[1]) for step in path])
                else:
                    if len(self.nodemap[last_node].Bnodes) > 0:
                        for nxt in self.nodemap[last_node].Bnodes:
                            if (nxt, '+') in path or (nxt, '-') in path:
                                #print("Have seen this node before: " + nxt, file=sys.stderr)
                                continue
                            if last_node in self.nodemap[nxt].Enodes:
                                stack.append(path + [(nxt, '-')])
                                #print("step forward", file=sys.stderr)
                            else:
                                stack.append(path + [(nxt, '+')])
                                #print("step forward", file=sys.stderr)
                    else:
                        #print("Reached tip " + str(len(paths)), file=sys.stderr)
                        if len(paths) < 10000:
                            paths.append(path)
                        else:
                            return paths
                        #print([(step[0][-7:-4], step[1]) for step in path])
        return paths

    def getPathSeqLength(self, path):
        rawLength = 0
        ovlpLength = 0
        for i in range(len(path) - 1):
            rawLength += self.nodemap[path[i][0]].length
            ovlpLength += self.edgeOvlp[(path[i][0], path[i+1][0])]
        rawLength += self.nodemap[path[-1][0]].length
        seqLength = rawLength - ovlpLength
        return seqLength

    def getPathSeq(self, path):
        if len(path) < 1:
            return ""
        if path[0][1] == '+':
            sequence = self.nodemap[path[0][0]].sequence
        else:
            sequence = rev_comp(self.nodemap[path[0][0]].sequence)
        for i in range(1, len(path)):
            node_sequence = self.nodemap[path[i][0]].sequence
            ovlpLength = self.edgeOvlp[(path[i-1][0], path[i][0])]
            if path[i][1] == '+':
                sequence += node_sequence[ovlpLength:]
            else:
                sequence += rev_comp(node_sequence)[ovlpLength:]
        return sequence


class Node(object):
    """docstring for Node"""
    def __init__(self, name, length, sequence):
        self.name = name
        self.length = length
        self.sequence = sequence
        self.Bnodes = set()
        self.Enodes = set()
