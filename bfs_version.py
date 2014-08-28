# read NMR restraints from mr.txt file
import os
import re
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
	list_next_points = []
	for key in dict_next_point:
		if dict_next_point[key] >= 3:
			list_next_points.append(key)
	if len(list_next_points) >0:
		return list_next_points
	else:
		for key in dict_next_point:
			if dict_next_point[key] >=2:
				list_next_points.append(key)
		return list_next_points


		
#list_search_growing = ['H_2', 'HA_2', 'HB2_2', 'QE_44']
#for i in xrange(30):
#	list_next_points = listing_next(list_search_growing, neighbors)
#	print list_next_points
#	for item in list_next_points:
#		list_search_growing.append(item)
#print len(list_search_growing)

#distance_matrix = fw_algo(nodes, distance_matrix)	

def iterative_procedure(vertices, distance_matrix, distance_connection):
	list_vertices = list(vertices)
	#list_c = list(vertices)
	dict_models = {}
	nb_model_build = 100
	nb_try = 0
	while nb_try < nb_model_build:
			try: 
				dict_coord = {}
				random_i = random.randint(0, len(vertices)-1)
				random_point = list(vertices)[random_i]
				list_search = breadth_first_search(distance_connection, random_point)
				list_4_points = list_search[:4]
				dict_coord = constructing_4_points(list_4_points, distance_matrix)
				for i in xrange(4, len(list_search)):
					next_point = list_search[i]
					#next_point = list_c[i]
					x_n, y_n, z_n = fifth_point(next_point, list_4_points, distance_matrix, dict_coord)
					dict_coord[next_point] = Point(next_point, x_n, y_n, z_n)
		
				V = []
				for key in list_vertices:
					x, y, z = dict_coord[key].x, dict_coord[key].y, dict_coord[key].z
					coord = [x,y,z]
					V.append(np.array(coord))
				V = np.array(V)
				V_c = centroid(V)
				V-=V_c
				dict_models[nb_try]=V
				nb_try+=1
				#for i in xrange(len(vertices)):
				#	print V[i][0], V[i][1],V[i][2]
			except ValueError: pass
			except ZeroDivisionError: pass
			except np.linalg.linalg.LinAlgError as err:
				if 'Singular matrix' in err.message:
					pass
				else: raise
	return dict_models












