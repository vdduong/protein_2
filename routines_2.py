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
				#random.shuffle(list_c)
				list_4_points = list_search[:4]
				#list_4_points = list_c[:4]
				dict_coord = constructing_4_points(list_4_points, distance_matrix)
				for i in xrange(4, len(list_search)):
				#for i in xrange(4, len(list_c)):
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

def depth_first_search(distance_matrix, start):
	blacknodes = []
	graynodes = [start]
	neighbors = {}
	for vertex in distance_matrix.keys():
		neighbors[vertex] = []
	for node_1 in distance_matrix.keys():
		for node_2 in distance_matrix.keys():
			neighbors[node_1].append(node_2)
	while graynodes:
		current = graynodes.pop()
		list_shuffle = []
		for neighbor in neighbors[current]:
			if not neighbor in blacknodes+graynodes:
				list_shuffle.append(neighbor)
		random.shuffle(list_shuffle)
		for node_shuffle in list_shuffle:
			graynodes.append(node_shuffle)
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
		


def graphing(nodes, model):
	list_nodes = list(nodes)
	dict_bb = {}
	list_bb = []
	for index, item in enumerate(list_nodes):
		if re.split('_', item)[0]=='H':
			list_bb.append(item)
			dict_bb[item] = model[index]

	list_item = []
	for item in list_bb:
		ind = int(re.split('_', item)[1])
		list_item.append(ind)
	list_item.sort()
	print list_item
	for index, item in enumerate(list_item):
		list_item[index] = 'H_' + str(item)

	print list_item
	X,Y,Z = [],[],[]
	for item in list_item:
		X.append(float(dict_bb[item][0]))
		Y.append(float(dict_bb[item][1]))
		Z.append(float(dict_bb[item][2]))

	fig = plt.figure()
	ax = Axes3D(fig)
	ax.plot(X,Y,Z,c='r', marker='o')
	plt.show()

def distance_matrix_generation(nodes, model):
	d = {}
	for node in list(nodes):
		d[node] = {}
	for index_1, node_1 in enumerate(list(nodes)):
		for index_2, node_2 in enumerate(list(nodes)):
			distance = sum((model[index_1][i]-model[index_2][i])**2 for i in xrange(3))
			d[node_1][node_2] = math.sqrt(distance)
	return d
#def harmonic_difference(nodes, matrix_1, matrix_2):
#	sum_harmonix = 0.0
#	time_harmonix = 0.0
#	for
#	return sum_harmonix/time_harmonix
	
	









