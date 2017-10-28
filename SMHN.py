import networkx as nx
import matplotlib.pyplot as plt
import csv
import pandas
import random
import pickle
import analysis as al
import numpy as np

G = nx.read_graphml("data/c.elegans.herm_pharynx_1.graphml")

probabilityData = {}
median = al.medianDegree(G)
hubFraction = al.nodeDegreeClassification(G, median)

def SMHN(G):
	SH = 0
	SN = 0
	MH = 0
	MN = 0
	for i,nbrs1 in G.adjacency_iter():
		if G.node[i]['role'] == 'S' and G.node[i]['degreeClass'] == 'High':
			SH += 1
		elif G.node[i]['role'] == 'S' and G.node[i]['degreeClass'] == 'Low':
			SN += 1
		elif G.node[i]['role'] == 'M' and G.node[i]['degreeClass'] == 'High':
			MH += 1
		elif G.node[i]['role'] == 'M' and G.node[i]['degreeClass'] == 'Low':
			MN += 1
	print SH
	print SN
	print MH
	print MN 		
	return 0

SMHN(G)