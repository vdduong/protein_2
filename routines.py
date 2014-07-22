# frequent routines for free assignment NMR protein structure determination

import numpy as np
import re
import math
import itertools

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
  d = (nd.linalg.det(V)*np.linalg.det(W)) < 0.0
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
  for V, W in zip(V, W):
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

def connected_components(nodes):
  result = []
  nodes = set(nodes)
  
  while nodes:
    n = nodes.pop()
    group = {n}
    queue = [n]
    
    while queue:
      n = queue.pop(0)
      neighbors = n.links # nodes that link directly to n
      neighbors.difference_update(group) # remove the neighbors we've already visited
      nodes.difference_update(neighbors) # remove the remaining nodes from the global set
      group.update(neighbors) # add them to the group of connected nodes
      queue.extend(neighbors) # add them to the queue, so we visit them in the next iterations
    result.append(group) # add the group to the list of groups
  return result 

def iterative_procedure():
  pass
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
