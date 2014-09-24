# auction algorithm

import math
from collections import deque

# given a price table and a value table
# set_a, set_b to be matched




def node_maximize(value_table, price_table, node_i):
	node_max = node_i
	max_table = []
	for node_j in price_table:
		max_table.append([node_j, value_table[node_i][node_j] - price_table[node_j]])
	max_table.sort(key=lambda x:x[1], reverse=True)
	node_max = max_table[0][0]
	return node_max

def auction_process(value_table, set_a, set_b):
	
	price_table = {}
	auction_table = {}
	for node in set_b:
		price_table[node] = 0.0
		auction_table[node] = None

	gamma = 1.0/float(len(set_b)+1)
		
	# iteration
	while set_a:
		node_i = set_a.pop()
		
		# find node_j that maximize value_table[i][j] - price_table[j]
		node_j = node_maximize(value_table, price_table, node_i)
		
		if value_table[node_i][node_j] - price_table[node_j] >= 0.0:
			# enqueue current owner j into set_a:
			if auction_table[node_j]!=None:
				set_a.append(auction_table[node_j])
			else:
				pass
			auction_table[node_j] = node_i
			price_table[node_j]+=gamma
		else: pass
		#print set_a
		#print '_'*20

	return auction_table

if __name__ == '__main__':
	set_a = [1,2,3]
	set_b = [3,2,1]

	value_table = {}
	for node in set_a:
		value_table[node] = {}

	for node in set_a:
		for node_ in set_b:
			value_table[node][node_] = 0.0
	
	value_table[1][2] = 3.
	value_table[2][2] = 3.
	value_table[3][3] = 3.
	value_table[2][1] = 3.
	value_table[3][1] = 4.
	value_table[1][3] = 4.
	auction = auction_process(value_table, set_a, set_b)
	print auction
