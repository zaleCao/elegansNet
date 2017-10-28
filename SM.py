import networkx as nx
import matplotlib.pyplot as plt
import csv
import pandas
import random
import pickle
import analysis as al
import numpy as np

G = nx.read_graphml("data/c.elegans.herm_pharynx_1.graphml")

timesteps = 500
simulation_no = 100
timedelay_range = 20
probabilityData = {}
median = al.medianDegree(G)
hubFraction = al.nodeDegreeClassification(G, median)

def SM_cal(G, dTE, timedelay):
	SM = {}
	for h in range(1, timedelay):
		dTE_S = 0
		dTE_M = 0
		SCount = 0
		MCount = 0 
		SM[h] = {}
		for i,nbrs1 in G.adjacency_iter():
			for j,nbrs2 in G.adjacency_iter():
				if h in dTE and i in dTE[h] and j in dTE[h][i]:
					if G.node[i]['role'] == 'S':
						dTE_S += dTE[h][i][j] 
						SCount += 1
					elif G.node[i]['role'] == 'M':
						dTE_M += dTE[h][i][j]  
						MCount += 1
		average_dTE_S = float(dTE_S)/float(SCount)
		average_dTE_M = float(dTE_M)/float(MCount)
		SM[h] = average_dTE_S - average_dTE_M
		print SM[h]
	return SM


with open('data/randomResults/dTEdata.txt', 'rb') as f:
	dTE = pickle.load(f)

SM = SM_cal(G, dTE, timedelay_range)


with open('data/randomResults/SM_E.txt', 'wb') as f:
	pickle.dump(SM, f)	
