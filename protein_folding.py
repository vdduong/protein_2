# protein folding using similated annealing algorithm
# simulated annealing does not assume the function is differentiable or continuous in its variables

# unfolded proteins = amino-acids are on a line.
# here, we assume only two types of amino-acids (H for hydrophobic and P for hydrophilic)
# and that each acid is a cube, all cubes have the same size, each two consecutive aa are connected at a face.

# limit the possible solid angles to 6 possibles values (0-5)

# our goal is arranging the cubes to minimize the number of faces of hydrophobic cubes that are exposed to water

import random
import math
import copy
from nlib import *

class Protein:
  
  moves = {0:(lambda x,y,z: (x+1,y,z)),
          1:(lambda x,y,z: (x-1,y,z)),
          2:(lambda x,y,z: (x,y+1,z)),
          3:(lambda x,y,z: (x,y-1,z)),
          4:(lambda x,y,z: (x,y,z+1)),
          5:(lambda x,y,z: (x,y,z-1))}
  
  def __init__(self, aminoacids):
    self.aminoacids = aminoacids
    self.angles = [0]*(len(aminoacids)-1)
    self.folding = self.compute_folding(self.angles)
    self.energy = self.compute_energy(self.folding)
  
  def compute_folding(self, angles):
    folding = []
    x,y,z = 0,0,0
    k = 0
    folding[x,y,z] = self.aminoacids[k]
    for angle in angles:
      k+=1
      xn,yn,zn = self.moves[angle](x,y,z)
      if (xn,yn,zn) in folding: return None # impossible folding
      folding[xn,yn,zn] = self.aminoacids[k]
      x,y,z = xn,yn,zn
    return folding
  
  def compute_energy(self, folding):
    E = 0
    for x,y,z in folding:
      aminoacid = folding[x,y,z]
      if aminoacid == 'H':
        for face in range(6):
          if not self.moves[face](x,y,z) in folding:
            E = E + 1
    return E
  
  def fold(self, t):
    while True:
      new_angles = copy.copy(self.angles)
      n = random.randint(1, len(self.aminoacids)-2)
      new_angles[n] = random.randint(0,5)
      new_folding = self.compute_folding(new_angles)
      if new_folding: break # found a valid folding
    new_energy = self.compute_energy(new_folding)
    if (self.energy-new_energy) > t*math.log(random.random()):
      self.angles = new_angles
      self.folding = new_folding
      self.energy = new_energy
    return self.energy
  

def main():
  aminoacids = ''.join(random.choice('HP' for k in range(20)))
  protein = Protein(aminoacids)
  t = 10.0
  while t > 1e-5:
    protein.fold(t=t)
    print(protein.energy, protein.angles)
    t = t*0.99 # cool
  
main()

# the moves dictionary is a dictionary of functions. For each solid angle, moves[angle] maps x,y,z, the coord of an aa
# to the coord of the cube at that solid angle
    
    
    
