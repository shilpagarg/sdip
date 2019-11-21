#!/usr/bin/python


import sys
import argparse
from Gfa import *

def getset(s):
    assert(len(s) == 1)
    for item in s:
        return item

def getStartOrEndNodes(graph):
    startOrEnd = set()
    for node_id, node in graph.nodes.items():
        if graph.edges[(node_id, True)] == set() or graph.edges[(node_id, False)] == set():
            startOrEnd.add(node_id)
    return startOrEnd

def getAllTips(graph, max_length):
    startOrEnd = getStartOrEndNodes(graph)
    #print(len(startOrEnd), file=sys.stderr)
    tips = []
    for tip in startOrEnd:
        #print("Start with tip " + tip, file=sys.stderr)
        last = tip
        stack = []
        if graph.edges[(last, True)] == set():
            stack.append([(last, '-')])
        else:
            stack.append([(last, '+')])
        while len(stack) > 0:
            path = stack.pop()
            last_node, last_orientation = path[-1]
            #print("At node " + last_node, file=sys.stderr)
            if last_orientation == '+':
                if len(graph.edges[(last_node, True)]) > 0:
                    if len(graph.edges[(last_node, True)]) == 1:
                        topos, infos = list(graph.edges[(last_node, True)])[0]
                        nxt, nxt_dir = topos
                        if nxt_dir == True:
                            if len(graph.edges[(nxt, False)]) > 1:
                                if len(path) <= max_length:
                                    tips.append(path)
                                break
                            else:
                                stack.append(path + [(nxt, '+')])
                        elif nxt_dir == False:
                            if len(graph.edges[(nxt, True)]) > 1:
                                if len(path) <= max_length:
                                    tips.append(path)
                                break
                            else:
                                stack.append(path + [(nxt, '-')])
                    else:
                        break
                else:
                    #print("Reached tip " + str(len(paths)), file=sys.stderr)
                    if len(path) <= max_length:
                        tips.append(path)
            else:
                if len(graph.edges[(last_node, False)]) > 0:
                    if len(graph.edges[(last_node, False)]) == 1:
                        topos, info = list(graph.edges[(last_node, False)])[0]
                        nxt, nxt_dir = topos
                        if nxt_dir == True:
                            if len(graph.edges[(nxt, False)]) > 1:
                                if len(path) <= max_length:
                                    tips.append(path)
                                break
                            else:
                                stack.append(path + [(nxt, '+')])
                        elif nxt_dir == False:
                            if len(graph.edges[(nxt, True)]) > 1:
                                if len(path) <= max_length:
                                    tips.append(path)
                                break
                            else:
                                stack.append(path + [(nxt, '-')])
                    else:
                        break
                else:
                    #print("Reached tip " + str(len(paths)), file=sys.stderr)
                    if len(path) <= max_length:
                        tips.append(path)
                    #print([(step[0][-7:-4], step[1]) for step in path])
    return tips

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('gfa', nargs='?', type=argparse.FileType('r'), default=sys.stdin, help = 'input GFA graph file (default: stdin)')
    parser.add_argument('action', choices=['list', 'remove'])
    parser.add_argument('--max_size', type = int, default = 1, help = 'maximum tip size to list or remove, default: 1')
    args = parser.parse_args()

    ingraph = args.gfa

    graph = Graph()
    graph.load(args.gfa)
    graph.remove_nonexistent_edges()

    print("Start finding tips..", file=sys.stderr)
    tips = getAllTips(graph, args.max_size)
    print("Found %d tips of length <= %d" % (len(tips), args.max_size), file=sys.stderr)

    if args.action == "list":
        for path in tips:
            print(",".join([node for node, direction in path]), file=sys.stdout)
    elif args.action == "remove":
        nodes_to_remove = set()
        for path in tips:
            for node, direction in path:
                nodes_to_remove.add(node)
        
        for node in nodes_to_remove:
            del graph.nodes[node]
        graph.remove_nonexistent_edges()
        graph.write(sys.stdout)

if __name__ == '__main__':
    main()
