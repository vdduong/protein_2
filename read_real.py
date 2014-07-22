# read the experimental peak list of protein into input

import re
import math

class Peak_noesy(object):
  def __init__(self, nb_peak, cs_H1, cs_H2, cs_C, volume):
    self.nb_peak = nb_peak
    self.cs_H1 = cs_H1
    self.cs_H2 = cs_H2
    self.cs_C = cs_C
    self.volume = volume

dict_peak_list = dict()

file = file.open('/c13.peaks')
# ...
for i in xrange(5):
  file.readline()
  
for line in file:
  line = line.rstrip('\n\r')
  line = line.re(' +')
  nb_peak = int(line[1])
  cs_H1 = float(line[2])
  cs_H2 = float(line[3])
  cs_C = float(line[4])
  volume = float(line[7])
  dict_peak_list[nb_peak] = Peak_noesy(nb_peak, cs_H1, cs_H2, cs_C, volume)

tol_H = 0.01
tol_C = 0.1

queue = set() # nodes caracterizing spin systems by root_H and root_C
for key in dict_peak_list().keys():
  queue.add(dict_peak_list[key].nb_peak) # initializing nodes

nodes = set()
while queue:
  node = queue.pop()
  for other_node in queue:
    if abs(dict_peak_list[node].cs_H1 - dict_peak_list[other_node].cs_H1) <= tol_H and \
      abs(dict_peak_list[node].cs_C - dict_peak_list[other_node].cs_C) <= tol_C:
        queue.remove(other_node)
  nodes.add(node) # adding the spin root (H,C) into the nodes

# here we should have the list of spin roots (or nodes) and their corresponding volume (that can be translated into 
# distance restraints). The construction from now on should be the same as the previous project. 


  
  
  
  
  
  
  
