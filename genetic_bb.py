# a genetic look-like algorithm to retrace the backbone ordering

def bacbone_tracing(adj_table, backbone_root):
  size_superposed = 1
  fragment_size = 3
  backbone_family = []
  for node_1 in adj_table:
    for node_2 in adj_table[node_1]:
        if adj_table[node_1][node_2] > 0.0 and node_1!=node_2:
          family_backbone.append([node_1,node_2])
          
  while fragment_size<=len(backbone_root):
    new_population =  []
    for member_1 in backbone_family:
      for member_2 in backbone_family:
        if member_1[size_superposed-1:] == member_2[:size_superposed-1]:
          new_fragment = member_1 + member_2[size_superposed:]
          if check_fragment(new_fragment): new_population.append(new_fragment)
    size_superposed+=1
    fragment_size+=1
    family_backbone = new_population
  return backbone

def check_fragment(fragment):
  '''check if the fragment has any duplicate inside it'''
  set_node_local = set()
  for node in fragment:
    set_node.add(node)
  if len(set_node)==len(fragment):
    return True
  else: return False
