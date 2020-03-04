#!/usr/bin/python
import re

class Node:
	def __init__(self):
		self.nodeid = 0
		self.nodeseq = ""
		self.length = 0
		self.readcount = 0
		self.frequency = 0
		self.chain = None

def reverse(pos):
	return ((pos[0], not pos[1]))

def reverse_cigar(cigar):
	#Reverse CIGAR string
	ops = re.findall(r'[0-9]+[MDI]', cigar)
	ops_parsed = [(int(op[:-1]), op[-1]) for op in ops]
	reversed_operations = {"M": "M", "D": "I", "I": "D"}
	ops_reversed = [(l, reversed_operations[op]) for l, op in ops_parsed[::-1]]
	cigar_reversed = "".join(["%d%s" % (l, op) for l, op in ops_reversed])
	return cigar_reversed

class Graph:
	def _parse_node(self, parts):
		result = Node()
		result.nodeid = parts[1]
		result.nodeseq = parts[2]
		result.length = 0
		result.readcount = 0
		result.frequency = 0
		for tag in parts[3:]:
			if tag[0:5] == 'LN:i:':
				result.length = int(tag[5:])
			elif tag[0:5] == 'RC:i:':
				result.readcount = int(tag[5:])
			elif tag[0:5] == 'km:f:':
				result.frequency = float(tag[5:])
		return result
	def __init__(self):
		self.nodes = {}
		self.edges = {}
	def load(self, f):
		for l in f:
			parts = l.strip().split('\t')
			if parts[0] == 'S':
				parsed = self._parse_node(parts)
				self.nodes[parsed.nodeid] = parsed
				if (parts[1], True) not in self.edges: self.edges[(parts[1], True)] = set()
				if (parts[1], False) not in self.edges: self.edges[(parts[1], False)] = set()
			if parts[0] == 'L':
				frompos = (parts[1], parts[2] == '+')
				topos = (parts[3], parts[4] == '+')
				overlap = int(parts[5][:-1])
				if frompos not in self.edges: self.edges[frompos] = set()
				if reverse(topos) not in self.edges: self.edges[reverse(topos)] = set()
				tags = parts[6:]
				cigar = tags[2]
				self.edges[reverse(topos)].add((reverse(frompos), (overlap, "\t".join(tags[0:2] + [reverse_cigar(cigar)] + tags[3:]))))
				self.edges[frompos].add((topos, (overlap, "\t".join(tags))))
	def remove_nonexistent_edges(self):
		extra = []
		for edge in self.edges:
			if edge[0] not in self.nodes:
				extra.append(edge)
				continue
			extra_here = []
			for target in self.edges[edge]:
				if target[0][0] not in self.nodes:
					extra_here.append(target)
					continue
			for e in extra_here:
				self.edges[edge].remove(e)
		for e in extra:
			del self.edges[e]
	def get_high_degree_nodes(self, min_degree=3):
		high_degree_nodes = []
		for n in self.nodes:
			targets = self.edges[(n, True)].union(self.edges[(n, False)])
			if len(targets) >= min_degree:
				high_degree_nodes.append(n)
		return high_degree_nodes
	def write(self, f):
		for node in self.nodes:
			n = self.nodes[node]
			line = "S\t" + str(n.nodeid) + "\t" + n.nodeseq + "\tLN:i:" + str(n.length)
			if n.readcount: line += '\tRC:i:' + str(n.readcount)
			if n.length: line += '\tkm:f:' + str(float(n.readcount)/float(n.length))
			if n.chain: line += "\tbc:Z:" + str(n.chain)
			print(line, file = f)
		for edge in self.edges:
			for target in self.edges[edge]:
				line = "L\t" + str(edge[0]) + "\t" + ("+" if edge[1] else "-") + "\t" + str(target[0][0]) + '\t' + ("+" if target[0][1] else "-") + '\t' + str(target[1][0]) + 'M' + "\t" + target[1][1]
				print(line, file = f)
