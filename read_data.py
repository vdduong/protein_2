# routines for assignment free protein structure determination

import numpy as np
import re
import os
import math
import itertools
import random
import matplotlib.pyplot as plt
import pylab

class Point(object):
  def __init__(self, name = '', x=0.0, y=0.0, z=0.0):
    self.name=name
    self.x=x
    self.y=y
    self.z=z


def kabsh(P, Q):
  """ Kabsh algorithm to align two point clouds
  return the rmsd between two structures and the aligned structure """
  C = np.dot(np.transpose(P), Q)
  V, S, W = np.linalg.svd(C)
  d = (np.linalg.det(V)*np.linalg.det(W)) < 0.0
  if (d):
    S[-1] = -S[-1]
    V[:,-1] = - V[:,-1]
  U = np.dot(V, W)
  P = np.dot(P, U)
  
  return rmsd(P, Q), P
  
def rmsd(V, W):
  """ return the rmsd from two sets of vectors V and W """
  D = len(V[0])
  N = len(V)
  rmsd = 0.0
  for v, w in zip(V, W):
    rmsd+=sum((v[i]-w[i])**2 for i in range(D))
  return np.sqrt(rmsd/N)

def centroid(X):
  """return the centroid from a vectorset X """
  C = sum(X)/len(X)
  return C

def fw_algo(vertices, dict_distance):
  """run the Ford algorithm to generate the complete triangular distance matrix """
  d = dict(dict_distance)
  for k in vertices:
    for i in vertices:
      for j in vertices:
        d[i][j] = min(d[i][j], d[i][k] + d[k][j])
  return d

def constructing_4_points(list_4_points, dict_distance):
  """contructs the coordinate of the first four points"""
  dict_coord = dict()
  atom_1 = list_4_points[0]
  x1,y1,z1 = 0.0, 0.0, 0.0
  dict_coord[atom_1] = Point(atom_1, x1, y1, z1)
  
  atom_2 = list_4_points[1]
  x2 = float('%.4f'%dict_distance[atom_1][atom_2])
  y2 = 0.0
  z2 = 0.0
  dict_coord[atom_2] = Point(atom_2, x2, y2, z2)
  
  atom_3 = list_4_points[2]
  d31 = dict_distance[atom_1][atom_3]
  d32 = dict_distance[atom_2][atom_3]
  x3 = float('%.4f'%((d31**2.0 - d32**2.0)/(2.0*x2) + x2/2.0))
  y3 = float('%.4f'%((d31**2.0 - x3**2.0)**(0.5)))
  z3 = 0.0
  dict_coord[atom_3] = Point(atom_3, x3, y3, z3)
  
  atom_4 = list_4_points[3]
  d41 = dict_distance[atom_1][atom_4]
  d42 = dict_distance[atom_2][atom_4]
  d43 = dict_distance[atom_3][atom_4]
  x4 = float('%.4f'%((d41**2.0 - d42**2.0)/(2.0*x2) + x2/2.0))
  y4 = float('%.4f'%((d42**2.0 - d43**2.0 - (x4-x2)**2.0 + (x4-x3)**2.0)/(2.0*y3) +y3/2.0))
  z4 = float('%.4f'%((d41**2.0 - x4**2.0 - y4**2.0)**0.5))
  dict_coord[atom_4] = Point(atom_4, x4, y4, z4)
  
  return dict_coord

def fifth_point(atom_5, list_4_points, dict_distance, dict_coord):
  atom_1 = list_4_points[0]
  atom_2 = list_4_points[1]
  atom_3 = list_4_points[2]
  atom_4 = list_4_points[3]
  x1, y1, z1 = dict_coord[atom_1].x, dict_coord[atom_1].y, dict_coord[atom_1].z
  x2, y2, z2 = dict_coord[atom_2].x, dict_coord[atom_2].y, dict_coord[atom_2].z
  x3, y3, z3 = dict_coord[atom_3].x, dict_coord[atom_3].y, dict_coord[atom_3].z
  x4, y4, z4 = dict_coord[atom_4].x, dict_coord[atom_4].y, dict_coord[atom_4].z
  
  d54 = dict_distance[atom_5][atom_4]
  d53 = dict_distance[atom_5][atom_3]
  d52 = dict_distance[atom_5][atom_2]
  d51 = dict_distance[atom_5][atom_1]
  
  a = np.array([[2.0*(x3-x4), 2.0*(y3-y4), 2.0*(z3-z4)],\
        [2.0*(x2-x4), 2.0*(y2-y4), 2.0*(z2-z4)], \
        [2.0*(x1-x4), 2.0*(y1-y4), 2.0*(z1-z4)]])
  b = np.array([d54**2.0 - d53**2.0 + x3**2.0 + y3**2.0 + z3**2.0 - x4**2.0 - y4**2.0 - z4**2.0, \
        d54**2.0 - d52**2.0 + x2**2.0 + y2**2.0 + z2**2.0 - x4**2.0 - y4**2.0 - z4**2.0, \
        d54**2.0 - d51**2.0 + x1**2.0 + y1**2.0 + z1**2.0 - x4**2.0 - y4**2.0 - z4**2.0])
  solve = np.linalg.solve(a,b)
  x5, y5, z5 = solve[0], solve[1], solve[2]
  
  x5 = float('%.4f'%x5)
  y5 = float('%.4f'%y5)
  z5 = float('%.4f'%z5)
  
  return x5, y5, z5

# clustering in 3D ? Given a list of points resulted from the construction
def distance_3d(point_1, point_2):
	sum = 0.0
	for i in xrange(3):
		sum+=(point_1[i] - point_2[i])**2
	return math.sqrt(sum)
	
def clustering(list_points):
	dict_clustering = dict()
	for index_point, point in enumerate(list_points):
		dict_clustering[index_point] = 1
	for index_1, point_1 in enumerate(list_points):
		for index_2, point_2 in enumerate(list_points):
			if index_2!=index_1 and distance_3d(point_1, point_2) <= 1.0: # threshold of 1. angstrom
				dict_clustering[index_1] +=1
	list_ordering = []
	for key_cluster in dict_clustering:
		list_ordering.append([key_cluster, dict_clustering[key_cluster]])
	list_ordering.sort(key=lambda x:x[1], reverse=True)
	try:
		if float(list_ordering[0][1]) >= float(len(list_points))/16.:
			return list_points[key_cluster]
		else:
			return '_'
	except IndexError:
		return '_'
		
def iterative_procedure(vertices, distance_matrix, distance_connection):
	list_vertices = list(vertices)
	dict_models = {}
	nb_model_build = 1
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
				nb_points_build = 0
				list_points = []
				if i<=9:
					for list_4_points_local in itertools.combinations(list_search[:i], 4):
						x_n, y_n, z_n = fifth_point(next_point, list_4_points_local, distance_matrix, dict_coord)
						list_points.append([x_n,y_n,z_n])	
					if clustering(list_points)!='_':
						x_n,y_n,z_n = clustering(list_points)[0],clustering(list_points)[1], clustering(list_points)[2]
						dict_coord[next_point] = Point(next_point,x_n,y_n,z_n)
						print next_point, x_n, y_n, z_n
					else:
						print "WRONG POINT"
						break
					
				else:
					limit_cluster = 200	
					list_random_local = [random.randrange(0,i) for _ in range(0,4)]
					list_4_points_local = [list_search[j] for j in list_random_local]
					if len(list_random_local) == len(set(list_random_local)):
						while nb_points_build < limit_cluster:
							x_n, y_n, z_n = fifth_point(next_point, list_4_points_local, distance_matrix, dict_coord)
							list_points.append([x_n,y_n,z_n])
							nb_points_build+=1
					else: pass
					if clustering(list_points)!='_':
						dict_coord[next_point] = Point(next_point, x_n, y_n, z_n)
						print next_point, x_n, y_n, z_n
					else: 
						print "WRONG POINT"
						break
			if len(dict_coord)==len(list_vertices):
				V =[]
				for key in list_vertices:
					x,y,z = dict_coord[key].x, dict_coord[key].y, dict_coord[key].z
					coord = [x,y,z]
					V.append(np.array(coord))
				V = np.array(V)
				V_c = centroid(V)
				V-=V_c
				dict_models[nb_try] = V
				nb_try+=1
			else:
				print "INCOMPLETE STRUCTURE"
				pass
			
		except ValueError: pass
		except ZeroDivisionError: pass
		except np.linalg.linalg.LinAlgError as err:
 			if 'Singular matrix' in err.message:
    			# your error handling block
  				pass
  			else: raise
	return dict_models
				
def breadth_first_search(distance_matrix, start):
	blacknodes = []
	graynodes = [start]
	neighbors = {}
	for vertex in distance_matrix.keys():
		neighbors[vertex] = []
	for node_1 in distance_matrix.keys():
		for node_2 in distance_matrix[node_1].keys():
			neighbors[node_1].append(node_2)
	while graynodes:
		current = graynodes.pop()
		list_shuffle = []
		for neighbor in neighbors[current]:
			if not neighbor in blacknodes+graynodes:
				list_shuffle.append(neighbor)
		random.shuffle(list_shuffle)
		for node_shuffle in list_shuffle:
			graynodes.insert(0, node_shuffle)
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


def distance_matrix_generation(nodes, model):
	d = {}
	for node in list(nodes):
		d[node] = {}
	for index_1, node_1 in enumerate(list(nodes)):
		for index_2, node_2 in enumerate(list(nodes)):
			distance = sum((model[index_1][i]-model[index_2][i])**2 for i in xrange(3))
			d[node_1][node_2] = math.sqrt(distance)
	return d

#### test module here 
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

distance_connection = dict(distance_matrix)

for node in nodes:
	distance_matrix[node][node] = 0.0
for node_1 in nodes:
	for node_2 in nodes:
		try:
			distance_ = distance_matrix[node_1][node_2]
		except KeyError:
			distance_matrix[node_1][node_2] = float('inf')
			distance_matrix[node_2][node_1] = float('inf')

distance_matrix = fw_algo(nodes, distance_matrix)	

dict_models = iterative_procedure(nodes, distance_matrix, distance_connection)

list_nodes = list(nodes)

def seq_distance(res_u, list_nodes, distance_matrix, model): # print the distances between res_u and res_u+1
	# let's say 35 is the number
	dict_current, dict_next = dict(), dict()
	for index_node, node in enumerate(list_nodes):
		if int(re.split('_', node)[1])==res_u:	dict_current[node] = index_node
		if int(re.split('_', node)[1])==res_u+1: dict_next[node] = index_node
	
	try:
		HA_current =  'HA_' + str(res_u)
		HN_next = 'H_' + str(res_u+1)
		distance_computed = (model[dict_current[HA_current]][0] - model[dict_next[HN_next]][0])**2 + \
			(model[dict_current[HA_current]][1] - model[dict_next[HN_next]][1])**2 + \
			(model[dict_current[HA_current]][2] - model[dict_next[HN_next]][2])**2
		distance_computed = math.sqrt(distance_computed)
		print 'distance HA(i)-HN(i+1) computed: %.4f and initial: %.4f'%(distance_computed, \
					distance_matrix[HA_current][HN_next])
	except KeyError: pass
	
	try:
		HN_current = 'H_' + str(res_u)
		HN_next = 'H_' + str(res_u+1)
		distance_computed = (model[dict_current[HN_current]][0] - model[dict_next[HN_next]][0])**2 + \
			(model[dict_current[HN_current]][1] - model[dict_next[HN_next]][1])**2 + \
			(model[dict_current[HN_current]][2] - model[dict_next[HN_next]][2])**2
		distance_computed = math.sqrt(distance_computed)
		print 'distance HN(i)-HN(i+1) computed: %.4f and initial: %.4f'%(distance_computed, \
					distance_matrix[HN_current][HN_next])
	except KeyError: pass
	
	try:
		HA_current = 'HA_' + str(res_u)
		HA_next = 'HA_' + str(res_u+1)
		distance_computed = (model[dict_current[HA_current]][0] - model[dict_next[HA_next]][0])**2 + \
			(model[dict_current[HA_current]][1] - model[dict_next[HA_next]][1])**2 + \
			(model[dict_current[HA_current]][2] - model[dict_next[HA_next]][2])**2
		distance_computed = math.sqrt(distance_computed)
		print 'distance HA(i)-HA(i+1) computed: %.4f and initial: %.4f'%(distance_computed, \
					distance_matrix[HA_current][HA_next])
	except KeyError: pass
	
	try:
		HN_current = 'H_' + str(res_u)
		HA_next = 'HA_' + str(res_u+1)
		distance_computed = (model[dict_current[HN_current]][0] - model[dict_next[HA_next]][0])**2 + \
			(model[dict_current[HN_current]][1] - model[dict_next[HA_next]][1])**2 + \
			(model[dict_current[HN_current]][2] - model[dict_next[HA_next]][2])**2
		distance_computed = math.sqrt(distance_computed)
		print 'distance HA(i)-HN(i+1) computed: %.4f and initial: %.4f'%(distance_computed, \
					distance_matrix[HN_current][HA_next])
	except KeyError: pass
	
	sum = 0.0
	nb_average = 0.0
	for node_current in dict_current:
		index_current = dict_current[node_current]
		for node_next in dict_next:
			index_next = dict_next[node_next]
			distance_computed = (model[index_current][0] - model[index_next][0])**2 + \
				(model[index_current][1] - model[index_next][1])**2 + \
				(model[index_current][2] - model[index_next][2])**2
			sum+= (math.sqrt(distance_computed)-distance_matrix[node_current][node_next])**2
			nb_average +=1.0
	print 'local residue-residue harmonic score : %.4f'%(sum/nb_average)
	print '_'*10

######### 
## Average model computation begins here

average_model = dict()
for node in list_nodes:
	average_model[node] = Point(node)
nb_average = 0.0
res_u = 35
for key_model in dict_models.keys():
	model = dict_models[key_model]
	d = distance_matrix_generation(nodes, model)
	sum_ = 0.0
	time_ = 0.0
	for node_1 in list(nodes):
		for node_2 in list(nodes):
			sum_ += (distance_matrix[node_1][node_2]-d[node_1][node_2])**2
			time_ +=1.
	harmonic_score = sum_/time_
	
	if harmonic_score <= 25.:
		print 'model score %.4f'%(harmonic_score)
		seq_distance(res_u, list_nodes, distance_matrix, model)
	
		if nb_average == 0.0: # first model to be added ?
			for index_node, node in enumerate(list_nodes):
				average_model[node].x += model[index_node][0]
				average_model[node].y += model[index_node][1]
				average_model[node].z += model[index_node][2]
			nb_average+=1.
			model_1 = model # assign model_1
		else:
			model_inverse = []
			for index_node, node in enumerate(list_nodes):
				coord_inverse = (model[index_node][0], -model[index_node][1], model[index_node][2])
				model_inverse.append(coord_inverse)
			model_inverse = np.array(model_inverse)
			rmsd_k, model = kabsh(model, model_1)
			rmsd_k_inverse, model_inverse = kabsh(model, model_1)
			
			if rmsd_k == min(rmsd_k, rmsd_k_inverse):
				for index_node, node in enumerate(list_nodes):
					average_model[node].x += model[index_node][0]
					average_model[node].y += model[index_node][1]
					average_model[node].z += model[index_node][2]
				nb_average +=1.0
			else:
				for index_node, node in enumerate(list_nodes):
					average_model[node].x += model_inverse[index_node][0]
					average_model[node].y += model_inverse[index_node][1]
					average_model[node].z += model_inverse[index_node][2]
				nb_average +=1.0

for node in list_nodes:
	average_model[node].x = average_model[node].x / nb_average
	average_model[node].y = average_model[node].y / nb_average
	average_model[node].z = average_model[node].z / nb_average


ref_model = dict()
file = open('/Users/james/Documents/protein/H_coord_1cpz')
for line in file:
	line = line.rstrip('\n\r')
	line = re.split(' +', line)
	ref_model[line[1]] = Point(line[1], float(line[2]), float(line[3]), float(line[4]))

dict_x,dict_y,dict_z = dict(),dict(),dict()
list_bb = []
P,P_inverse,Q = [],[],[]
for index_node, node in enumerate(list_nodes):
	if re.split('_', node)[0]=='H':
		list_bb.append(node)
		dict_x[node] = model[index_node][0]
		dict_y[node] = model[index_node][1]
		dict_z[node] = model[index_node][2]
		coord_exp = [dict_x[node], dict_y[node], dict_z[node]]
		coord_exp_inverse = [dict_x[node], -dict_y[node], dict_z[node]]
		P.append(coord_exp)
		P_inverse.append(coord_exp_inverse)
		coord_ref = [ref_model[node].x, ref_model[node].y, ref_model[node].z]
		Q.append(coord_ref)
P = np.array(P)
Q = np.array(Q)
rmsd_k, P = kabsh(P,Q)
rmsd_k_inverse,P_inverse = kabsh(P_inverse,Q)
inversion = False # to know if the structure should be reversed or not ?
if rmsd_k < rmsd_k_inverse:		
	print rmsd_k
	for index_node, node in enumerate(list_bb):
		dict_x[node] = P[index_node][0]
		dict_y[node] = P[index_node][1]
		dict_z[node] = P[index_node][2]
else:
	print rmsd_k_inverse
	inversion = True
	for index_node, node in enumerate(list_bb):
		dict_x[node] = P_inverse[index_node][0]
		dict_y[node] = P_inverse[index_node][1]
		dict_z[node] = P_inverse[index_node][2]
list_ordering = []
for node in list_bb:
	ind = int(re.split('_', node)[1])
	list_ordering.append(ind)
list_ordering.sort()
for index_node, node in enumerate(list_ordering):
	list_ordering[index_node] = 'H_' + str(node)
		
X,Y,Z = [],[],[]
X_ref,Y_ref,Z_ref = [],[],[]
		
for node in list_ordering:
	X.append(float(dict_x[node]))
	Y.append(float(dict_y[node]))
	Z.append(float(dict_z[node]))
	X_ref.append(ref_model[node].x)
	Y_ref.append(ref_model[node].y)
	Z_ref.append(ref_model[node].z)
		
fig = plt.figure()		
ax = Axes3D(fig)
ax.plot(X,Y,Z,c='r',marker='o')
ax.plot(X_ref,Y_ref,Z_ref,c='b', marker='o')
plt.show()		
				
				
				
				
