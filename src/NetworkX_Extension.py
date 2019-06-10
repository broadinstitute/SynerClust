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
				if source not in distances[w]:
					distances[w][source] = {}
				if target not in distances[w]:
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


def merge(graph, old_graph, left, right, middle):
	graph.add_node(middle)
	# commons = old_graph[left].keys() and old_graph[right].keys()
	commons = [l for l in old_graph[left].keys() if l in old_graph[right].keys()]
	for common in commons:
		graph.add_edge(common, middle, rank=max(old_graph[common][left]['rank'], old_graph[common][right]['rank']))
	for target in old_graph[left]:
		if target not in commons and target != right:
			graph.add_edge(target, middle, rank=old_graph[target][left]['rank'], identity=old_graph[target][left]['identity'])
	for target in old_graph[right]:
		if target not in commons and target != left:
			graph.add_edge(target, middle, rank=old_graph[target][right]['rank'], identity=old_graph[target][right]['identity'])
	graph.remove_node(left)
	graph.remove_node(right)


# def merge_lists_identicals(graph, old_graph, left, right, middle):
def merge_lists_identicals(graph, old_graph, nodes, middle):
	graph.add_node(middle)
	commons = set()
	for node in nodes:
		tmp_commons = [l for l in old_graph[node].keys()]
		commons.update(tmp_commons)

	for common in commons:
		max_rank = max([old_graph[common][node]['rank'] for node in nodes if node in old_graph[common]])
		graph.add_edge(common, middle, rank=max_rank)

	for node in nodes:
		for target in old_graph[node]:
			if target not in commons and target not in nodes and middle not in graph[middle]:
				graph.add_edge(target, middle, rank=old_graph[target][node]['rank'], identity=old_graph[target][node]['identity'])

	for node in nodes:
		graph.remove_node(node)
