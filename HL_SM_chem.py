import networkx as nx
import matplotlib.pyplot as plt
import csv
import pandas
import random
import pickle
import numpy as np
from statistics import median

G = nx.read_graphml("data/c.elegans.herm_pharynx_1.graphml")

timesteps = 500
simulation_no = 100
timedelay_range = 20
probabilityData = {}


def nodeDegreeClassification(G, median):
	averageEdges = G.number_of_edges() / G.number_of_nodes()
	highCounter = 0
	for n,nbrs in G.adjacency_iter():
		if G.degree(n) > median:
			G.node[n]['degreeClass'] = 'High'
			highCounter += 1
		else:
			G.node[n]['degreeClass'] = 'Low'
	hubFraction = float(highCounter) / G.number_of_nodes()
	return hubFraction, highCounter

def medianDegree(G):
	degrees = [0] * G.number_of_nodes()
	i = 0
	for n,nbrs in G.adjacency_iter():
		degrees[i] = G.degree(n)
		i += 1
	return median(degrees)

median = medianDegree(G)
print median
hubFraction = nodeDegreeClassification(G, median)
with open('data/randomResults/dTEdata_chem.txt', 'rb') as f:
	dTE = pickle.load(f)

def SM(G, dTE, timedelay):
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
		if SCount == 0:
			SM[h] = 0
		else:				
			average_dTE_S = float(dTE_S)/float(SCount)
			average_dTE_M = float(dTE_M)/float(MCount)
			SM[h] = average_dTE_S - average_dTE_M
			print SM[h]
	return SM

SM = SM(G, dTE, timedelay_range)

x = SM.keys()
y = SM.values()

fig, ax = plt.subplots()
line1, = ax.plot(x, y, linewidth=1,)
plt.axis([0, 20, -0.10, 0.50])

fig.suptitle('Sensor - Motor Value ', fontsize=14)
plt.ylabel('SM Value')
plt.xlabel('Time Delay h')

plt.show()

'''
def HL(G, dTE, timedelay):
	HL = {}
	for h in range(1, timedelay):
		dTE_high = 0
		dTE_low = 0
		highCount = 0
		lowCount = 0
		HL[h] = {} 
		for i,nbrs1 in G.adjacency_iter():
			for j,nbrs2 in G.adjacency_iter():
				if h in dTE and i in dTE[h] and j in dTE[h][i]:
					#print G.node[i]['degreeClass']
					if G.node[i]['degreeClass'] == 'High':
						dTE_high += dTE[h][i][j] 
						highCount += 1
					else:
						dTE_low += dTE[h][i][j]  
						lowCount += 1
		if highCount == 0:
			HL[h] = 0
		else:
			average_dTE_high = float(dTE_high)/float(highCount)
			average_dTE_low = float(dTE_low)/float(lowCount)
			HL[h] = average_dTE_high - average_dTE_low
	return HL

HL = HL(G, dTE, timedelay_range)
with open('data/randomResults/HLdata_chem.txt', 'wb') as f:
	pickle.dump(HL, f)

x = HL.keys()
y = HL.values()

fig, ax = plt.subplots()
line1, = ax.plot(x, y, linewidth=1,)
plt.axis([0, 20, -0.20, 0.20])

fig.suptitle('Hub - Non-Hub Value ', fontsize=14)
plt.ylabel('HN Value')
plt.xlabel('Time Delay h')

plt.show()

'''
