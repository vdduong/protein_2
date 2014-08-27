# read NMR restraints from mr.txt file
import os
import re
from routines import *
import matplotlib.pyplot as plt
import pylab
from mpl_toolkits.mplot3d import Axes3D

path_local = os.getcwd()

#nmr_restraints = file(path_local + '/2KA0.mr.txt')
nmr_restraints = file(path_local + '/1CPZ.mr.txt')
for nb_line in xrange(25):
	nmr_restraints.readline()

nodes = set()
distance_matrix = {}

for line in nmr_restraints:
	if line!='\n':
		line = line.rstrip('\n\r')
		data_line = re.split(' +', line)
		try:
			int(data_line[0])
			node_1 = data_line[2] + '_' + data_line[0]
			node_2 = data_line[5] + '_' + data_line[3]
			if (node_1 in nodes) and (node_2 not in nodes):
				distance_matrix[node_1][node_2] = float(data_line[6])
				nodes.add(node_2)
				distance_matrix[node_2] = {}
				distance_matrix[node_2][node_1] = distance_matrix[node_1][node_2]
			elif (node_2 in nodes) and (node_1 not in nodes):
				distance_matrix[node_2][node_1] = float(data_line[6])
				nodes.add(node_1)
				distance_matrix[node_1] = {}
				distance_matrix[node_1][node_2] = distance_matrix[node_2][node_1]
			elif (node_1 not in nodes) and (node_2 not in nodes):
				nodes.add(node_1)
				nodes.add(node_2)
				distance_matrix[node_1], distance_matrix[node_2] = {},{}
				distance_matrix[node_1][node_2] = float(data_line[6])
				distance_matrix[node_2][node_1] = float(data_line[6])
			else:
				distance_matrix[node_1][node_2] = float(data_line[6])
				distance_matrix[node_1][node_2] = float(data_line[6])
		except ValueError:
			node_1 = data_line[3] + '_' + data_line[1]
			node_2 = data_line[6] + '_' + data_line[4]
			if (node_1 in nodes) and (node_2 not in nodes):
				distance_matrix[node_1][node_2] = float(data_line[7])
				nodes.add(node_2)
				distance_matrix[node_2] = {}
				distance_matrix[node_2][node_1] = distance_matrix[node_1][node_2]
			elif (node_2 in nodes) and (node_1 not in nodes):
				distance_matrix[node_2][node_1] = float(data_line[7])
				nodes.add(node_1)
				distance_matrix[node_1] = {}
				distance_matrix[node_1][node_2] = distance_matrix[node_2][node_1]
			elif (node_2 not in nodes) and (node_1 not in nodes):
				nodes.add(node_1)
				nodes.add(node_2)
				distance_matrix[node_1], distance_matrix[node_2] = {},{}
				distance_matrix[node_1][node_2] = float(data_line[7])
				distance_matrix[node_2][node_1] = float(data_line[7])
			else:
				distance_matrix[node_1][node_2] = float(data_line[7])
				distance_matrix[node_2][node_1] = float(data_line[7])
	else:
		break

nodes = connected_components(nodes, distance_matrix)

nodes_del = set(distance_matrix.keys()).difference(nodes)
for node_1 in nodes_del:
		del distance_matrix[node_1]

distance_connection = dict(distance_matrix)

for node in nodes:
	distance_matrix[node][node] = 0.0
for node_1 in nodes:
	for node_2 in nodes:
		try:
			distance_ = distance_matrix[node_1][node_2]
		except KeyError:
			distance_matrix[node_1][node_2] = float('inf')
			distance_matrix[node_2][node_1] = float('inf')
print len(nodes)
distance_matrix = fw_algo(nodes, distance_matrix)	

dict_models = iterative_procedure(nodes, distance_matrix, distance_connection)

list_nodes = list(nodes)

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

######### 
## Average model computation begins here

average_model = dict()
for node in list_nodes:
	average_model[node] = Point(node)
nb_average = 0.0
res_u = 35
for key_model in dict_models.keys():
	model = dict_models[key_model]
	d = distance_matrix_generation(nodes, model)
	sum_ = 0.0
	time_ = 0.0
	for node_1 in list(nodes):
		for node_2 in list(nodes):
			sum_ += (distance_matrix[node_1][node_2]-d[node_1][node_2])**2
			time_ +=1.
	harmonic_score = sum_/time_
	
	if harmonic_score <= 20.:
		print 'model score %.4f'%(harmonic_score)
		seq_distance(res_u, list_nodes, distance_matrix, model)
	
		if nb_average == 0.0: # first model to be added ?
			for index_node, node in enumerate(list_nodes):
				average_model[node].x += model[index_node][0]
				average_model[node].y += model[index_node][1]
				average_model[node].z += model[index_node][2]
			nb_average+=1.
			model_1 = model # assign model_1
		else:
			model_inverse = []
			for index_node, node in enumerate(list_nodes):
				coord_inverse = (model[index_node][0], -model[index_node][1], model[index_node][2])
				model_inverse.append(coord_inverse)
			model_inverse = np.array(model_inverse)
			rmsd_k, model = kabsh(model, model_1)
			rmsd_k_inverse, model_inverse = kabsh(model, model_1)
			
			if rmsd_k == min(rmsd_k, rmsd_k_inverse):
				for index_node, node in enumerate(list_nodes):
					average_model[node].x += model[index_node][0]
					average_model[node].y += model[index_node][1]
					average_model[node].z += model[index_node][2]
				nb_average +=1.0
			else:
				for index_node, node in enumerate(list_nodes):
					average_model[node].x += model_inverse[index_node][0]
					average_model[node].y += model_inverse[index_node][1]
					average_model[node].z += model_inverse[index_node][2]
				nb_average +=1.0

for node in list_nodes:
	average_model[node].x = average_model[node].x / nb_average
	average_model[node].y = average_model[node].y / nb_average
	average_model[node].z = average_model[node].z / nb_average


ref_model = dict()
file = open('/Users/james/Documents/protein/H_coord_1cpz')
for line in file:
	line = line.rstrip('\n\r')
	line = re.split(' +', line)
	ref_model[line[1]] = Point(line[1], float(line[2]), float(line[3]), float(line[4]))

dict_x,dict_y,dict_z = dict(),dict(),dict()
list_bb = []
P,P_inverse,Q = [],[],[]
for index_node, node in enumerate(list_nodes):
	if re.split('_', node)[0]=='H':
		list_bb.append(node)
		dict_x[node] = model[index_node][0]
		dict_y[node] = model[index_node][1]
		dict_z[node] = model[index_node][2]
		coord_exp = [dict_x[node], dict_y[node], dict_z[node]]
		coord_exp_inverse = [dict_x[node], -dict_y[node], dict_z[node]]
		P.append(coord_exp)
		P_inverse.append(coord_exp_inverse)
		coord_ref = [ref_model[node].x, ref_model[node].y, ref_model[node].z]
		Q.append(coord_ref)
P = np.array(P)
Q = np.array(Q)
rmsd_k, P = kabsh(P,Q)
rmsd_k_inverse,P_inverse = kabsh(P_inverse,Q)
inversion = False # to know if the structure should be reversed or not ?
if rmsd_k < rmsd_k_inverse:		
	print rmsd_k
	for index_node, node in enumerate(list_bb):
		dict_x[node] = P[index_node][0]
		dict_y[node] = P[index_node][1]
		dict_z[node] = P[index_node][2]
else:
	print rmsd_k_inverse
	inversion = True
	for index_node, node in enumerate(list_bb):
		dict_x[node] = P_inverse[index_node][0]
		dict_y[node] = P_inverse[index_node][1]
		dict_z[node] = P_inverse[index_node][2]
list_ordering = []
for node in list_bb:
	ind = int(re.split('_', node)[1])
	list_ordering.append(ind)
list_ordering.sort()
for index_node, node in enumerate(list_ordering):
	list_ordering[index_node] = 'H_' + str(node)
		
X,Y,Z = [],[],[]
X_ref,Y_ref,Z_ref = [],[],[]
		
for node in list_ordering:
	X.append(float(dict_x[node]))
	Y.append(float(dict_y[node]))
	Z.append(float(dict_z[node]))
	X_ref.append(ref_model[node].x)
	Y_ref.append(ref_model[node].y)
	Z_ref.append(ref_model[node].z)
		
fig = plt.figure()		
ax = Axes3D(fig)
ax.plot(X,Y,Z,c='r',marker='o')
ax.plot(X_ref,Y_ref,Z_ref,c='b', marker='o')
plt.show()

output = open('/Users/james/Documents/protein/pre_structure_random','w')
for node in list_nodes:
	x,y,z = average_model[node].x, average_model[node].y, average_model[node].z
	output.write('  %s  %.4f  %.4f  %.4f\n'%(node,x,y,z))
output.close()
