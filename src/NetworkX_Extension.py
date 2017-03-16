#!/usr/bin/env python

import networkx as nx


def all_pairs_path_length(G, weight=['weight']):
	"""Compute path lengths in the weighted graph (acyclic). Allows for negative weights.
	Returns a dictionary of dictionaries with path[source][target]=L
	"""

	# TODO add a verification that the graph is really acyclic

	nodes = G.nodes()
	paths = nx.shortest_path(G)
	attributes = []
	for w in weight:
		attributes.append(nx.get_edge_attributes(G, w))
	# attributes = nx.get_edge_attributes(G, weight)
	distances = []
	for w in xrange(len(weight)):
		distances.append({})
		s = 0
		while s < len(nodes):
			t = s
			source = nodes[s]
			while t < len(nodes):
				target = nodes[t]
				path = paths[source][target]
				i = 0
				if source not in distances:
					distances[w][source] = {}
				if target not in distances:
					distances[w][target] = {}
				distances[w][source][target] = 0.0
				while i < len(path) - 1:
					if (path[i], path[i + 1]) in attributes[0]:
						distances[w][source][target] += attributes[w][path[i], path[i + 1]]
					else:
						distances[w][source][target] += attributes[w][path[i + 1], path[i]]
					i += 1
				distances[w][target][source] = distances[w][source][target]
				t += 1
			s += 1
	return (distances, paths)
