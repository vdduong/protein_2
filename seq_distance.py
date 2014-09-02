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
