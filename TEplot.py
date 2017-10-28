import pickle
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt


with open('data/randomResults/TEdata.txt', 'rb') as f:
	TE = pickle.load(f)

with open('data/randomResults/hopcountData.txt', 'rb') as f:
	HC = pickle.load(f)

G = nx.read_graphml("data/c.elegans.herm_pharynx_1.graphml")

tePlot = [[0 for x in range(19)] for y in range(7)] 
hopcounter = [0] * 7
for h in range(1,20):
	TEadd = 0
	n = 0 
	for i,nbrs1 in G.adjacency_iter():
		for j,nbrs2 in G.adjacency_iter():
			if i in HC and j in HC[i]:
				hop = HC[i][j]
				if hop in range(7):
					tePlot[hop][h-1] +=TE[h][i][j] 
					hopcounter[hop] += 1 
print tePlot
for hop in range(7):
	for h in range(1,20):
		tePlot[hop][h-1] = tePlot[hop][h-1] / hopcounter[hop]

print hopcounter
print tePlot

fig, ax = plt.subplots()
plt.plot(range(1,20), tePlot[1], 'bx')
plt.plot(range(1,20), tePlot[2], 'g.')
plt.plot(range(1,20), tePlot[3], 'r.')
plt.plot(range(1,20), tePlot[4], 'c.')
plt.plot(range(1,20), tePlot[5], 'm.')
plt.plot(range(1,20), tePlot[6], 'k.')
#{'b', 'g', 'r', 'c', 'm', 'y', 'k', 'w'}

plt.axis([0, 20, 0, 0.0045])

fig.suptitle('Transfer entropy over different time delays ', fontsize=14)
plt.ylabel('Averaged transfer entropy (TE) of all node pairs per hopcount')
plt.xlabel('Time Delay h')

labels = []
for i in range(1, 7):
    labels.append(r'%i - hop' %i)

# I'm basically just demonstrating several different legend options here...
plt.legend(labels, ncol=4, loc='upper center', 
           bbox_to_anchor=[0.5, 1.0], 
           columnspacing=1.0, labelspacing=0.0,
           handletextpad=0.0, handlelength=1.5,
           fancybox=True, shadow=True)

plt.show()	



