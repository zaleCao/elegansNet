import networkx as nx
import matplotlib.pyplot as plt
import csv
import pandas
import random
import pickle
import analysis as al
try:
    import pygraphviz
    from networkx.drawing.nx_agraph import graphviz_layout
except ImportError:
    try:
        import pydotplus
        from networkx.drawing.nx_pydot import graphviz_layout
    except ImportError:
        raise ImportError("This example needs Graphviz and either "
                          "PyGraphviz or PyDotPlus")

G = nx.read_graphml("data/c.elegans.herm_pharynx_1.graphml")

timesteps = 50
simulation_no = 3
timedelay_range = 5
probabilityData = {}
median = al.medianDegree(G)
hubFraction = al.nodeDegreeClassification(G, median)

with open('data/results/dTEdata.txt', 'rb') as f:
	pickle.load(f)	

HL = al.HL(G, dTE, timedelay_range)

with open('data/results/activityData.txt', 'rb') as f:
	pickle.load(f)

with open('data/randomResults/HLdata.txt', 'rb') as f:
	pickle.load(f)	

with open('data/results/dieDownTime.txt', 'rb') as f:
	pickle.load(f)

with open('data/results/died.txt', 'rb') as f:
	pickle.load(f)

with open('data/randomResults/frequencies.txt', 'rb') as f:
	pickle.load(f)	

def nodeDegreeClassification(G):
	averageEdges = G.number_of_edges() / G.number_of_nodes()
	highCounter = 0
	for n,nbrs in G.adjacency_iter():
		if G.degree(n) > 11:
			G.node[n]['degreeClass'] = 'High'
			highCounter += 1
		else:
			G.node[n]['degreeClass'] = 'Low'
	hubFraction = float(highCounter) / G.number_of_nodes()
	print hubFraction