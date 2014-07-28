# read NMR restraints from mr.txt file
# test the backbone tracing approach by the HN graph

import os
import re
import random
from routines import *
import matplotlib.pyplot as plt
import pylab
from mpl_toolkits.mplot3d import Axes3D

path_local = os.getcwd()

#nmr_restraints = file(path_local + '/2KA0.mr.txt')
nmr_restraints = file(path_local + '/1CPZ.mr.txt')
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

nodes = connected_components(nodes, distance_matrix)

nodes_del = set(distance_matrix.keys()).difference(nodes)
for node_1 in nodes_del:
		del distance_matrix[node_1]
#for node in nodes:
#	distance_matrix[node][node] = 0.0
			
adj_table = {} # table of adjacency between HN, expressing the scores of connectivity
nodes_HN = set()

for node in nodes:
	if re.split('_',node)[0]=='H':
		adj_table[node] = {}
		nodes_HN.add(node)

def prob_adj(HN_1, HN_2, distance_matrix):
	neighbors = {}
	neighbors[HN_1] = distance_matrix[HN_1].keys()
	neighbors[HN_2] = distance_matrix[HN_2].keys()
	commun = len(set(neighbors[HN_1]).intersection(set(neighbors[HN_2])))
	flow_out = len(set(neighbors[HN_1]))
	if HN_2 in neighbors[HN_1]:
		commun+=2
	else:pass
	adj_table[HN_1][HN_2] = float(commun)/float(flow_out)
	return adj_table[HN_1][HN_2]


nb_constraints = 0
for HN_1 in nodes_HN:
	for HN_2 in nodes_HN:
		if HN_1 !=HN_2:
			adj_table[HN_1][HN_2] = prob_adj(HN_1, HN_2, distance_matrix)
			#if adj_table[HN_1][HN_2] > 0.0:
			nb_constraints+=1
			#	print HN_1, HN_2, adj_table[HN_1][HN_2]

def chain_score(chain, adj_table):
	score = 0.0
	for i in xrange(len(chain)-1):
		score+=adj_table[chain[i]][chain[i+1]]
	return score

def dfs(start, adj_table, nodes_HN):
	blacknodes = []
	graynodes = [start]
	neighbors = {}
	for node in nodes_HN:
		neighbors[node] = []
	for node_1 in adj_table.keys():
		for node_2 in adj_table[node_1].keys():
			if adj_table[node_1][node_2] > 0.0 and node_1 !=node_2:
				neighbors[node_1].append(node_2)
	while graynodes:
		current = graynodes.pop()
		for neighbor in neighbors[current]:
			if not neighbor in blacknodes+graynodes:
				graynodes.append(neighbor)
		blacknodes.append(current)
	return blacknodes

# recursive adding new node to an established chain ???
def score_adding_node(node, chain, adj_table):
	'''
	return the score sequence according to the addition of node into the sequence
	(depending on position in the sequence)
	'''
	score = [0.0 for i in xrange(len(chain)+1)]
	for i in xrange(0, len(chain)-1):
		score[i+1]+=adj_table[chain[i]][node] + adj_table[node][chain[i+1]] - adj_table[chain[i]][chain[i+1]]
	score[0]=adj_table[node][chain[0]]
	score[len(chain)]=adj_table[chain[len(chain)-1]][node]
	return score

def chains_adding(node, chains, adj_table):
	'''a list of 10 best previous chains
	following the addition of the new node,
	gives a new list of 10 best new chains'''
	old_chains = list(chains)
	list_chains_score = []
	while old_chains:
		chain_test = old_chains.pop()
		score = score_adding_node(node, chain, adj_table)
		score_base = chain_score(chain_test, adj_table)
		for i in xrange(len(score)):
			new_chain = list(chain_test)
			new_chain.insert(i, node)
			list_chains_score.append((new_chain,score_base+score[i]))
	list_chains_score.sort(key=lambda x:x[1], reverse=True)
	#for item in list_chains_score: print item[0], item[1]
	nb_lim_chain = min(len(list_chains_score), 10) # limiting the number of considered chains to 10
	nb_added = 0
	new_chains = []
	while nb_added < nb_lim_chain:
		new_chains.append(list_chains_score[nb_added][0])
		nb_added+=1
	return new_chains

list_alpha = [int(re.split('_',node)[1]) for node in nodes_HN]
list_alpha.sort()
list_beta = ['H_'+str(item) for item in list_alpha]
#print list_beta
#print chain_score(list_beta, adj_table)

#chain = ['H_66', 'H_67'] # beginning chain
#nodes_residual = set(nodes_HN)
nodes_residual = list(list_beta)
nodes_residual = set(nodes_HN)
chain = random.sample(nodes_residual, 2)
#chain = ['H_63','H_66']
for item in chain:
	nodes_residual.remove(item)
	
chains = [chain]
i = 0
while nodes_residual:	
	node = nodes_residual.pop()
	chains = chains_adding(node, chains, adj_table)
	#for item in chains:
	#	print item, chain_score(item, adj_table)

for item in chains:
	print item, chain_score(item, adj_table)
	print '\n'

print chain_score(list_beta, adj_table)



#blacknodes = dfs('H_4', adj_table, nodes_HN)
#list_alpha = [int(re.split('_',node)[1]) for node in nodes_HN]
#list_alpha.sort()
#list_beta = ['H_'+str(item) for item in list_alpha]
#print list_beta
#print chain_score(list_beta, adj_table)





