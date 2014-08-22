def backbone_tracing(adj_table):
	size_superposed = 1
	member_size = 3
	population = []
	for node_1 in adj_table:
		for node_2 in adj_table:
			if adj_table[node_1][node_2] > 0.0 and node_1!=node_2:
				population.append([node_1,node_2])
	print len(population)
	while member_size <= len(adj_table):
		new_population = []
		for member_1 in population:
			for member_2 in population:
				if member_1[len(member_1)-size_superposed:]==member_2[:size_superposed]:
					new_member = member_1 + member_2[size_superposed:]
					if check_member(new_member):
						new_population.append(new_member)
		size_superposed+=1
		member_size+=1
		population = list(new_population)
		#for member in population:
		#	print member, chain_score(member,adj_table)
		print len(population)
		print '_________________'
	return population

def check_member(member):
	'''check if the fragment has any duplicate inside'''
	set_node_local = set()
	for node in member:
		set_node_local.add(node)
	if len(set_node_local)==len(member):
		return True
	else: return False
