# backbone building, based on the main chain tracing model ACR for X ray diffraction protein

def distance_between_protons(H_1,H_2):
  '''
  compute distance between two protons
  '''
  return math.sqrt((H_1.x-H_2.x)**2 + (H_1.y-H_2.y)**2 + (H_1.z-H_2.z)**2)

def pairing_HN():
  return None

# consider only the network of HN protons (backbone+sidechain) and the initial distance matrix
# given those information, one creates ordered groups of four connected protons
# then compares those groups to form the longest possible chain of protons ?
