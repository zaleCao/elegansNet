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

timesteps = 500
simulation_no = 100
timedelay_range = 20
probabilityData = {}
median = al.medianDegree(G)
hubFraction = al.nodeDegreeClassification(G, median)

with open('data/randomResults/activityData_chem.txt', 'rb') as f:
	activitydata = pickle.load(f)

with open('data/randomResults/hopcountData_chem.txt', 'rb') as f:
	hopcountdata = pickle.load(f)

with open('data/randomResults/died_chem.txt', 'rb') as f:
	died = pickle.load(f)
'''
for i in range(simulation_no):
	probabilityData[i] = al.probability_calulations(activitydata[i], timesteps, timedelay_range)
'''

probabilityData = al.all_itr_probability_calulation(G, activitydata, simulation_no, timesteps, timedelay_range, died)
with open('data/randomResults/probabilityData_chem.txt', 'wb') as f:
	pickle.dump(probabilityData, f)

#with open('data/randomResults/probabilityData.txt', 'rb') as f:
#	probabilityData =	pickle.load(f)		

TE = al.TE(G, probabilityData, timedelay_range)
with open('data/randomResults/TEdata_chem.txt', 'wb') as f:
	pickle.dump(TE, f)

#with open('data/randomResults/TEdata.txt', 'rb') as f:
#	TE =	pickle.load(f)		

dTE = al.dTE(G, TE, timedelay_range)
with open('data/randomResults/dTEdata_chem.txt', 'wb') as f:
	pickle.dump(dTE, f)

#with open('data/randomResults/dTEdata.txt', 'rb') as f:
#	dTE =	pickle.load(f)	
HL = al.HL(G, dTE, timedelay_range)
with open('data/randomResults/HLdata_chem.txt', 'wb') as f:
	pickle.dump(HL, f)	

#with open('data/randomResults/HLdata.txt', 'rb') as f:
#	HL=	pickle.load(f)	

#store data	