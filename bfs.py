# breadth first search

def breadth_first_search(graph, start):
  vertices, links = graph
  blacknodes = []
  graynodes = [start]
  neighbors = [[] for vertex in vertices]
  for link in links:
    neighbor[link[0]].append(link[1])
  while graynodes:
    current = graynodes.pop()
    for neighbor in neighbors[current]:
      if not neighbor in blacknodes+graynodes:
        graynodes.insert(0, neighbor)
    blacknodes.append(current)
  return blacknodes
