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
adj_table = {}
nodes_HN = set()

for node in nodes:
  if re.split('_',node)[0]=='H':
    adj_table[node]={}
    nodes_HN.add(node)
adj_table = {}
for node in nodes:
  adj_table[node] = {}

def prob_adj(HN_1, HN_2, distance_matrix):
  neighbors = {}
  neighbors[HN_1] = distance_matrix[HN_1].keys()
  neighbors[HN_2] = distance_matrix[HN_2].keys()
  commun = len(set(neighbors[HN_1]).intersection(set(neighbors[HN_2])))
  flow_out = len(set(neighbors[HN_1]))
  if HN_2 in neighbors[HN_1]:
    commun+=2
  adj_table[HN_1][HN_2] = float(commun)/float(flow_out)
  return adj_table[HN_1][HN_2]

def longest_chain(adj_table, nodes):
  '''define the longest chains of HN found from table of adjacency'''
  ### longest transverse ?
  
  return None

for HN_1 in nodes_HN:
  for HN_2 in nodes_HN:
    adj_table[HN_1][HN_2] = prob_adj(HN_1, HN_2, distance_matrix)
    #if adj_table[HN_1][HN_2] > 0.0:
    #  print HN_1, HN_2, adj_table[HN_1][HN_2]
def chain_score(chain,adj_table):
  global len_prot # length of the protein in term of residus
  score = 0.0
  if len(chain) > len_prot: raise ValueError('length of the chain is more than what is necessary')
  else:
    for i in xrange(len(chain)-1):
      score+=adj_table[chain[i]][chain[j]]
  return score

# recursive processus ?
# think about the method of putting another proton into a chain of i protons already ordered 
# nodes_HN, adj_table, chain

def add_node(node, chain, adj_table):
  score = [[] for i in xrange(len(chain)+1)]
  for i in xrange(0, len(chain)-1):
    score[i+1]+=adj_table[chain[i]][node] + adj_table[node][chain[i+1]]
  score[0]=adj_table[node][chain[0]]
  score[len(chain)]=adj_table[chain[len(chain)-1]][node]
  for i in xrange(len(score)):
    print i, score[i]



    
    
