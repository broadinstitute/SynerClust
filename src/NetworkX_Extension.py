#!/usr/bin/env python

import networkx as nx

def all_pairs_path_length(G, weight='weight'):
	"""Compute path lengths in the weighted graph (acyclic). Allows for negative weights.
	Returns a dictionary of dictionaries with path[source][target]=L
	"""

	# TODO add a verification that the graph is really acyclic
	
	nodes = G.nodes()
	paths = nx.shortest_path(G)
	distances = {}
	s = 0
	while s < len(nodes):
		t = s
		source = nodes[s]
		while t < len(nodes):
			target = nodes[t]
			path = paths[source][target]
			i = 0
			if not source in distances:
				distances[source] = {}
			if not target in distances:
				distances[target] = {}
			distances[source][target] = 0.0
			attributes = nx.get_edge_attributes(G, weight)
			while i < len(path)-1:
				if (path[i], path[i+1]) in attributes:
					distances[source][target] += attributes[path[i], path[i+1]]
				else:
					distances[source][target] += attributes[path[i+1], path[i]]
				i+=1
			distances[target][source] = distances[source][target]
			t+=1
		s+=1

	return distances
