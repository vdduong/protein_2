# read NMR restraints from mr.txt file
# test the backbone tracing approach by the HN graph

# NOE secondary structure
# alpha helices & 3,10-helix & beta sheet anti parallel & parallel


import os
import re
import random
from routines import *
from simulated_annealing import *
import matplotlib.pyplot as plt
import pylab
from mpl_toolkits.mplot3d import Axes3D
from collections import deque
from auction_algo import *

path_local = os.getcwd()

#nmr_restraints = file(path_local + '/2KA0.mr.txt')
nmr_restraints = file(path_local + '/1CPZ.mr.txt')
#nmr_restraints = file(path_local + '/2KL2.mr.txt')
#nmr_restraints = open(path_local + '/1CPZ.new_restraints.txt')

for nb_line in xrange(25):
	nmr_restraints.readline()

nodes = set()
distance_matrix = {}

for line in nmr_restraints:
	if line!='\n':
		line = line.rstrip('\n\r')
		data_line = re.split(' +', line)
		try:
			int(data_line[0])
			node_1 = data_line[2] + '_' + data_line[0]
			node_2 = data_line[5] + '_' + data_line[3]
			if (node_1 in nodes) and (node_2 not in nodes):
				distance_matrix[node_1][node_2] = float(data_line[6])
				nodes.add(node_2)
				distance_matrix[node_2] = {}
				distance_matrix[node_2][node_1] = distance_matrix[node_1][node_2]
			elif (node_2 in nodes) and (node_1 not in nodes):
				distance_matrix[node_2][node_1] = float(data_line[6])
				nodes.add(node_1)
				distance_matrix[node_1] = {}
				distance_matrix[node_1][node_2] = distance_matrix[node_2][node_1]
			elif (node_1 not in nodes) and (node_2 not in nodes):
				nodes.add(node_1)
				nodes.add(node_2)
				distance_matrix[node_1], distance_matrix[node_2] = {},{}
				distance_matrix[node_1][node_2] = float(data_line[6])
				distance_matrix[node_2][node_1] = float(data_line[6])
			else:
				distance_matrix[node_1][node_2] = float(data_line[6])
				distance_matrix[node_1][node_2] = float(data_line[6])
		except ValueError:
			node_1 = data_line[3] + '_' + data_line[1]
			node_2 = data_line[6] + '_' + data_line[4]
			if (node_1 in nodes) and (node_2 not in nodes):
				distance_matrix[node_1][node_2] = float(data_line[7])
				nodes.add(node_2)
				distance_matrix[node_2] = {}
				distance_matrix[node_2][node_1] = distance_matrix[node_1][node_2]
			elif (node_2 in nodes) and (node_1 not in nodes):
				distance_matrix[node_2][node_1] = float(data_line[7])
				nodes.add(node_1)
				distance_matrix[node_1] = {}
				distance_matrix[node_1][node_2] = distance_matrix[node_2][node_1]
			elif (node_2 not in nodes) and (node_1 not in nodes):
				nodes.add(node_1)
				nodes.add(node_2)
				distance_matrix[node_1], distance_matrix[node_2] = {},{}
				distance_matrix[node_1][node_2] = float(data_line[7])
				distance_matrix[node_2][node_1] = float(data_line[7])
			else:
				distance_matrix[node_1][node_2] = float(data_line[7])
				distance_matrix[node_2][node_1] = float(data_line[7])
	else:
		break

def d_v(s,t,u,v,distance_matrix):
	'''distance function on the nodes of C_uv'''
	try:
		distance_st = distance_matrix[s][t]
		return 0.5
	except KeyError:
		set_s = set(distance_matrix[s].keys())
		set_t = set(distance_matrix[t].keys())
		set_intersect = set_s.intersection(set_t)
		if len(set_intersect) > 2: return 0
		else: return 1

def dispersion(u,v,distance_matrix):
	disp = 0
	set_u = set(distance_matrix[u].keys())
	set_v = set(distance_matrix[v].keys())
	set_commun = set_u.intersection(set_v)
	for s in set_commun:
		for t in set_commun:
			if s!=t:
				disp+=d_v(s,t,u,v,distance_matrix)
	if len(set_commun) > 0:
		return disp/float(len(set_commun))
	else:
		return 0.0

def rec_disp(u,emb,distance_matrix,disp):
	neighbors = []
	set_commun = dict()
	set_u = set(distance_matrix[u].keys())
	for node in distance_matrix[u].keys():
		if node !=u: 
			neighbors.append(node)
			set_node = set(distance_matrix[node].keys())
			set_unode = set_u.intersection(set_node)
			set_commun[node] = set_unode
	
	for node in neighbors: disp[u][node] = 1.0	# reinitialising the dispersion dictionary
	dict_value = dict()
	
	for nb_iteration in xrange(3):
		for v in neighbors:
			sum_nominator = 0.0
			set_uv = set_commun[v]
			for s in set_uv:
				for t in set_uv:
					if s!=t:
						sum_nominator+=2.*d_v(s,t,u,v,distance_matrix)
			for w in neighbors:
				sum_nominator+=disp[u][w]**2.
			try:
				score_v = sum_nominator/float(emb[u][v])
				dict_value[v] = score_v
			except ZeroDivisionError:
				dict_value[v] = 0.0
		for v in neighbors:
			disp[u][v] = dict_value[v]
	max = 0.0
	for v in disp[u].keys():
		if disp[u][v] > max: max = disp[u][v]
	if max == 0.0: 
		pass
	else:
		for v in disp[u].keys():
			disp[u][v] = disp[u][v]/max
	return disp

nodes = connected_components(nodes, distance_matrix)

nodes_del = set(distance_matrix.keys()).difference(nodes)
for node_1 in nodes_del:
		del distance_matrix[node_1]

emb, disp = dict(), dict()
for node in nodes:
	emb[node] = dict()
for node_1 in nodes:
	set_1 = set(distance_matrix[node_1].keys())
	for node_2 in nodes:
		set_2 = set(distance_matrix[node_2].keys())
		set_commun = set_1.intersection(set_2)
		emb[node_1][node_2] = len(set_commun)

for node in nodes:
	disp[node] = dict()
	disp = rec_disp(node,emb,distance_matrix,disp) 

backbone_root = []
for node in nodes:
	if re.split('_',node)[0] == 'HN' or re.split('_', node)[0]=='H': backbone_root.append(node)

backbone_residue = dict((u,[]) for u in backbone_root) # pseudo-residue 

for u in backbone_residue:
	for v in disp[u].keys():
		if disp[u][v]>0.00 and distance_matrix[u][v] <=8.0 and v not in backbone_residue:
			backbone_residue[u].append(v)

for u in backbone_residue:
	backbone_residue[u].append(u)


adj_table = {}
for node in backbone_residue:
	adj_table[node] = {}
	for node_ in backbone_residue:
		adj_table[node][node_] = 0.0

def prob_adj(HN_1, HN_2, distance_matrix, backbone_residue):
	connections = 0.0
	root_1 = backbone_residue[HN_1]
	root_2 = backbone_residue[HN_2]
	for node_1 in root_1:
		for node_2 in root_2:
			try:
				distance = distance_matrix[node_1][node_2]
				connections+=1.
			except KeyError:
				pass
	if len(root_1)*len(root_2) > 0.0:
		adj_table[HN_1][HN_2] = connections/(float(len(root_1)*len(root_2))) 
	else:
		adj_table[HN_1][HN_2] = 0.0
	return adj_table[HN_1][HN_2]

def adjacency(current_node, adj_table, chain_current):
	'''return the adjacency of current node'''
	neighbors = []
	for node in adj_table[current_node]:
		if adj_table[current_node][node] > 0.0 and node!= current_node and node not in chain_current:
			neighbors.append([node, adj_table[current_node][node]])
	neighbors.sort(key=lambda x:x[1], reverse=True)
	if len(neighbors) > 0: return neighbors[0][0]
	else: return '_'

def chain_score(chain, adj_table):
	score = 0.0
	for i in xrange(len(chain)-1):
		score+=adj_table[chain[i]][chain[i+1]]
	return score

for u_1 in backbone_root:
	for u_2 in backbone_root:
		adj_table[u_1][u_2] = prob_adj(u_1, u_2, distance_matrix, backbone_residue)


for u in adj_table:
	adj_table[u][u] = 0.0

for u in adj_table:
	for v in adj_table[u]:
		if adj_table[u][v]>0.0 : #and u!=v:
			nb_u = int(re.split('_',u)[1])
			nb_v = int(re.split('_',v)[1])
			if nb_u > nb_v:
				adj_table[u][v] = 0.0
			#print u, v, adj_table[u][v]

set_a = list(adj_table.keys())
set_b = list(adj_table.keys())
buffer = 'buffer'
adj_table[buffer] = {}

for u in adj_table:
	adj_table[u][buffer] = 0.1
	adj_table[buffer][u] = 0.1
set_a.append(buffer)
set_b.append(buffer)

auction = auction_process(adj_table, set_a, set_b)
print len(auction), 'connections'
set_ = set()
for key in auction:
	print auction[key], key

backbone_order = ['buffer']

while len(backbone_order)<len(set_b):
	for u in auction:
		if u==auction[backbone_order[-1]]:
			backbone_order.append(u)

print backbone_order
