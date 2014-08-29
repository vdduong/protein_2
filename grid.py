# read NMR restraints from mr.txt file
import os
import re
import random
import matplotlib.pyplot as plt
import pylab
from mpl_toolkits.mplot3d import Axes3D
import itertools
import numpy as np
import math
from routines import *

def breadth_first_search(neighbors, start):
	blacknodes = []
	graynodes = [start]
	while graynodes:
		current = graynodes.pop()
		for neighbor in neighbors[current]:
			if not neighbor in blacknodes+graynodes:
				graynodes.insert(0, neighbor)
		blacknodes.append(current)
	return blacknodes

def connected_components(nodes, distance_matrix):
	#result = []
	nodes = set(nodes)
	max_group = set()
	while nodes:
		start = nodes.pop()
		bfs_group = breadth_first_search(distance_matrix, start)
		nodes = nodes.difference(set(bfs_group))
		if len(bfs_group) > len(max_group): max_group = set(bfs_group)
		else: pass
		#result.append(set(bfs_group))
	if len(max_group) > 0.9*len(nodes): return max_group
	else: 
		print 'Error: network not dense enough'
		return None
		

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

distance_connection = dict()
distance_connection = distance_matrix
for node in nodes:
	distance_matrix[node][node] = 0.0
for node_1 in nodes:
	for node_2 in nodes:
		try:
			distance_ = distance_matrix[node_1][node_2]
		except KeyError:
			distance_matrix[node_1][node_2] = float('inf')
			distance_matrix[node_2][node_1] = float('inf')
print len(nodes), 'nodes'

neighbors = dict()
for node in nodes:
	neighbors[node] = []
for node in nodes:
	for node_prime in nodes:
		if node!=node_prime and 0.0 < distance_matrix[node][node_prime] < float('inf'):
			neighbors[node].append(node_prime)

def listing_next(list_search_growing, neighbors):
	dict_next_point = dict()
	set_commun_neighbor = set()
	for node in list_search_growing:
		for neighbor in neighbors[node]:
			if neighbor not in set(list_search_growing):
				set_commun_neighbor.add(neighbor)
	
	for commun_neighbor in set_commun_neighbor:
		if commun_neighbor not in set(list_search_growing):
			dict_next_point[commun_neighbor] = 0
	
	for node in list_search_growing:
		for neighbor in neighbors[node]:
			if neighbor in set_commun_neighbor:
				dict_next_point[neighbor]+=1
	
	list_ordering = []
	for neighbor in dict_next_point:
		list_ordering.append([neighbor, dict_next_point[neighbor]])
	list_ordering.sort(key=lambda x:x[1], reverse=True)
	try:
		max_to_be_added = list_ordering[0][1]
		list_next_points = []
		for item in list_ordering:
			if item[1] == max_to_be_added:
				list_next_points.append(item[0])
		return list_next_points
	except IndexError:
		return []


		
list_search_growing = ['H_2', 'HA_2', 'HB2_2', 'QE_44']
dict_4_points = dict() # dictionary gathers respective lists of 4 points used for the construction of each next point

def building_list(list_search_growing, dict_4_points, neighbors):
	list_next_nodes = listing_next(list_search_growing, neighbors)
	for next_node in list_next_nodes:
		dict_4_points[next_node] = set()
	for node in list_search_growing:
		for next_node in list_next_nodes:
			if node in neighbors[next_node]:
				dict_4_points[next_node].add(node)
	
	# random adding point to insufficient list
	for next_node in list_next_nodes:
		while len(dict_4_points[next_node]) <4:
			random_node = list_search_growing[random.randint(0,len(list_search_growing)-1)]
			dict_4_points[next_node].add(random_node)
	
	for next_node in list_next_nodes:
		list_search_growing.append(next_node)
	
	return list_search_growing, dict_4_points


def grid_search(atom_5, sub_points, distance_matrix, dict_coord):
	nb_grid = 10
	
	atom_1 = sub_points[0]
	atom_2 = sub_points[1]
	atom_3 = sub_points[2]
	atom_4 = sub_points[3]
	x, y , z = [],[],[]
	x1, y1, z1 = dict_coord[atom_1].x, dict_coord[atom_1].y, dict_coord[atom_1].z
	x.append(x1)
	y.append(y1)
	z.append(z1)
	x2, y2, z2 = dict_coord[atom_2].x, dict_coord[atom_2].y, dict_coord[atom_2].z
	x.append(x2)
	y.append(y2)
	z.append(z2)
	x3, y3, z3 = dict_coord[atom_3].x, dict_coord[atom_3].y, dict_coord[atom_3].z
	x.append(x3)
	y.append(y3)
	z.append(z3)
	x4, y4, z4 = dict_coord[atom_4].x, dict_coord[atom_4].y, dict_coord[atom_4].z
	x.append(x4)
	y.append(y4)
	z.append(z4)
	
	distance_to = []
	d1 = distance_matrix[atom_5][atom_1]
	d2 = distance_matrix[atom_5][atom_2]
	d3 = distance_matrix[atom_5][atom_3]
	d4 = distance_matrix[atom_5][atom_4]	
	distance_to.append(d1)
	distance_to.append(d2)
	distance_to.append(d3)
	distance_to.append(d4)
	
	bord_min_x = max(-distance_to[i]+x[i] for i in xrange(len(x)))
	bord_max_x = min(distance_to[i]+x[i] for i in xrange(len(x)))
	bord_min_y = max(-distance_to[i]+y[i] for i in xrange(len(y)))
	bord_max_y = min(distance_to[i]+y[i] for i in xrange(len(y)))
	bord_min_z = max(-distance_to[i]+z[i] for i in xrange(len(z)))
	bord_max_z = min(distance_to[i]+z[i] for i in xrange(len(z)))
	
	grid_x = [bord_min_x + float(i)*(bord_max_x-bord_min_x)/float(nb_grid) for i in xrange(nb_grid+1)]
	grid_y = [bord_min_y + float(i)*(bord_max_y-bord_min_y)/float(nb_grid) for i in xrange(nb_grid+1)]
	grid_z = [bord_min_z + float(i)*(bord_max_z-bord_min_z)/float(nb_grid) for i in xrange(nb_grid+1)]
	
	x_init, y_init, z_init = grid_x[0], grid_y[0], grid_z[0]
	err_initial = 0.0
	for i in xrange(4):
		err_initial += abs(-(x_init - x[i])**2 - (y_init- y[i])**2 - (z_init - z[i])**2 + (distance_to[i])**2)
	
	x5,y5,z5 = x_init, y_init, z_init
	
	for x_ in grid_x:
		for y_ in grid_y:
			for z_ in grid_z:
				err = 0.0
				for i in xrange(4):
					err+=abs(-(x_ - x[i])**2 - (y_ - y[i])**2 - (z_ - z[i])**2 + (distance_to[i])**2)
				#print x_, y_, z_, err
				if err < err_initial:
					x5, y5, z5 = x_, y_, z_
					err_initial = err
				else: pass
	return x5, y5, z5

def consensus_coord(node_constructed, dict_4_points, distance_matrix, dict_coord):
	sum_x, sum_y, sum_z = 0.0, 0.0, 0.0
	nb_average = 0.0
	for subset in itertools.combinations(dict_4_points[node_constructed],4):
		x_, y_, z_ = grid_search(node_constructed, subset, distance_matrix, dict_coord)
		sum_x+=x_
		sum_y+=y_
		sum_z+=z_
		nb_average+=1.
	if nb_average > 0.0:
		x, y, z = sum_x/nb_average, sum_y/nb_average, sum_z/nb_average
		dict_coord[node_constructed] = Point(node_constructed, x, y, z)
		return dict_coord
	else:
		print 'Error'
	
	
distance_matrix = fw_algo(nodes, distance_matrix)
dict_coord = constructing_4_points(list_search_growing, distance_matrix)

while len(list_search_growing) < len(nodes):
#for i in xrange(1):
	list_search_growing, dict_4_points = building_list(list_search_growing, dict_4_points, neighbors)
for node_constructed in list_search_growing[4:]:
	dict_coord = consensus_coord(node_constructed, dict_4_points, distance_matrix, dict_coord)


for node in list_search_growing:
	print node, dict_coord[node].x, dict_coord[node].y, dict_coord[node].z

model = []
for key in list(nodes):
	x, y, z = dict_coord[key].x, dict_coord[key].y, dict_coord[key].z
	coord = [x,y,z]
	model.append(np.array(coord))
model = np.array(model)
model_c = centroid(model)
model-=model_c

new_matrix = distance_matrix_generation(nodes, model)
harmonix = 0.0
nb_ = 0.0
for node in distance_matrix:
	for node_ in distance_matrix[node]:
		harmonix+= (distance_matrix[node][node_] - new_matrix[node][node_])**2
		nb_+=1
print nb_
print harmonix/nb_





