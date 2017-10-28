import networkx as nx
from math import log
from statistics import median

def all_itr_probability_calulation(G,activitydata, simulation_no, timesteps, timedelay, died):

	#create an nested list for each pair of nodes in each iteration in each simulation and each number of time delay	
	#TE = [[[[[[]for j in range(G.number_of_nodes())]for i in range(G.number_of_nodes())]for z in range(timedelay)]for y in range(timesteps)]for x in range(simulation_no)]
	TE = {} 

	permutations = ([0,0,0],[1,0,0],[0,1,0],[0,0,1],[1,1,0],[1,0,1],[0,1,1],[1,1,1])
	probabilities = {}
	for h in range(1, timedelay):
		probabilities[h] = {}
		for i,nbrs1 in G.adjacency_iter():
			probabilities[h][i] = {}
			for j,nbrs2 in G.adjacency_iter():
				probabilities[h][i][j] = {}

				#initialise case counters
				KLM000 = 0
				KLM100 = 0
				KLM010 = 0
				KLM001 = 0
				KLM110 = 0
				KLM101 = 0
				KLM011 = 0
				KLM111 = 0

				LM00 = 0
				LM10 = 0
				LM01 = 0
				LM11 = 0

				KL00 = 0
				KL10 = 0
				KL01 = 0
				KL11 = 0

				L0 = 0
				L1 = 0

				for itr in range(simulation_no):
					if itr not in died:
						timestart = timesteps - 200
						for t in range (timestart, timesteps):
							tDelayed = t - h  
							#all the permutations of the first probability term of TE formula
							if activitydata[itr][t][j] == 0 and activitydata[itr][tDelayed][j] == 0 and activitydata[itr][tDelayed][i] == 0:
								KLM000 +=1
							if activitydata[itr][t][j] >= 1 and activitydata[itr][tDelayed][j] == 0 and activitydata[itr][tDelayed][i] == 0:
								KLM100 +=1
							if activitydata[itr][t][j] == 0 and activitydata[itr][tDelayed][j] >= 1 and activitydata[itr][tDelayed][i] == 0:
								KLM010 +=1
							if activitydata[itr][t][j] == 0 and activitydata[itr][tDelayed][j] == 0 and activitydata[itr][tDelayed][i] >= 1:
								KLM001 +=1
							if activitydata[itr][t][j] >= 1 and activitydata[itr][tDelayed][j] >= 1 and activitydata[itr][tDelayed][i] == 0:
								KLM110 +=1
							if activitydata[itr][t][j] >= 1 and activitydata[itr][tDelayed][j] == 0 and activitydata[itr][tDelayed][i] >= 1:
								KLM101 +=1
							if activitydata[itr][t][j] == 0 and activitydata[itr][tDelayed][j] >= 1 and activitydata[itr][tDelayed][i] >= 1:
								KLM011 +=1
							if activitydata[itr][t][j] >= 1 and activitydata[itr][tDelayed][j] >= 1 and activitydata[itr][tDelayed][i] >= 1:
								KLM111 +=1

							#all the permutations of the conditional probability term components,lm,of TE formula	
							if activitydata[itr][tDelayed][j] == 0 and activitydata[itr][tDelayed][i] == 0:
								LM00 +=1
							if activitydata[itr][tDelayed][j] >= 1 and activitydata[itr][tDelayed][i] == 0:
								LM10 +=1
							if activitydata[itr][tDelayed][j] == 0 and activitydata[itr][tDelayed][i] >= 1:
								LM01 +=1
							if activitydata[itr][tDelayed][j] >= 1 and activitydata[itr][tDelayed][i] >= 1:
								LM11 +=1
							#all the permutations of the conditional probability term components,kl,of TE formula	
							if activitydata[itr][t][j] == 0 and activitydata[itr][tDelayed][j] == 0:
								KL00 +=1
							if activitydata[itr][t][j] >= 1 and activitydata[itr][tDelayed][j] == 0:
								KL10 +=1
							if activitydata[itr][t][j] == 0 and activitydata[itr][tDelayed][j] >= 1:
								KL01 +=1
							if activitydata[itr][t][j] >= 1 and activitydata[itr][tDelayed][j] >= 1:
								KL11 +=1
							#all the permutations of the conditional probability term components,l,of TE formula		
							if activitydata[itr][tDelayed][j] == 0:
								L0 +=1
							if activitydata[itr][tDelayed][j] >= 1:
								L1 +=1
				#all the permutations of the conditional probability term components,klm,of TE formula	
				probabilities[h][i][j]['000'] = float(KLM000)/(timestart)/simulation_no
				probabilities[h][i][j]['100'] = float(KLM100)/(timestart)/simulation_no
				probabilities[h][i][j]['010'] = float(KLM010)/(timestart)/simulation_no
				probabilities[h][i][j]['001'] = float(KLM001)/(timestart)/simulation_no
				probabilities[h][i][j]['110'] = float(KLM110)/(timestart)/simulation_no
				probabilities[h][i][j]['101'] = float(KLM101)/(timestart)/simulation_no
				probabilities[h][i][j]['011'] = float(KLM011)/(timestart)/simulation_no
				probabilities[h][i][j]['111'] = float(KLM111)/(timestart)/simulation_no
				#all the permutations of the conditional probability term components,lm,of TE formula	
				probabilities[h][i][j]['lm00'] = float(LM00)/(timestart)/simulation_no
				probabilities[h][i][j]['lm10'] = float(LM10)/(timestart)/simulation_no
				probabilities[h][i][j]['lm01'] = float(LM01)/(timestart)/simulation_no
				probabilities[h][i][j]['lm11'] = float(LM11)/(timestart)/simulation_no
				#all the permutations of the conditional probability term components,kl,of TE formula	
				probabilities[h][i][j]['kl00'] = float(KL00)/(timestart)/simulation_no
				probabilities[h][i][j]['kl10'] = float(KL10)/(timestart)/simulation_no
				probabilities[h][i][j]['kl01'] = float(KL01)/(timestart)/simulation_no
				probabilities[h][i][j]['kl11'] = float(KL11)/(timestart)/simulation_no
				#all the permutations of the conditional probability term components,l,of TE formula		
				probabilities[h][i][j]['l0'] = float(L0)/(timestart)/simulation_no
				probabilities[h][i][j]['l1'] = float(L1)/(timestart)/simulation_no

	return probabilities

'''			

def probability_calulation(G,activitydata, simulation_no, timesteps, timedelay):

	#create an nested list for each pair of nodes in each iteration in each simulation and each number of time delay	
	#TE = [[[[[[]for j in range(G.number_of_nodes())]for i in range(G.number_of_nodes())]for z in range(timedelay)]for y in range(timesteps)]for x in range(simulation_no)]
	TE = {} 

	#initialise case counters
	KLM000 = 0
	KLM010 = 0
	KLM001 = 0
	KLM110 = 0
	KLM101 = 0
	KLM011 = 0
	KLM111 = 0

	LM00 = 0
	LM10 = 0
	LM01 = 0
	LM11 = 0

	KL00 = 0
	KL10 = 0
	KL01 = 0
	KL11 = 0

	L0 = 0
	L1 = 0


	permutations = ([0,0,0],[1,0,0],[0,1,0],[0,0,1],[1,1,0],[1,0,1],[0,1,1],[1,1,1])
	probabilities = {}
	for h in range(timedelay):
		probabilities[h] = {}
		for i,nbrs1 in G.adjacency_iter():
			probabilities[h][i] = {}
			for j,nbrs2 in G.adjacency_iter():
				probabilities[h][i][j] = {}
				for t in range (timesteps):
					tDelayed = t - h  

					if tDelayed < 0:
						pass
					#all the permutations of the first probability term of TE formula
					elif activitydata[t][j] == 0 and activitydata[tDelayed][j] == 0 and activitydata[tDelayed][i] == 0:
						KLM000 +=1
					elif activitydata[t][j] ==1 and activitydata[tDelayed][j] == 0 and activitydata[tDelayed][i] == 0:
						KLM100 +=1
					elif activitydata[t][j] == 0 and activitydata[tDelayed][j] == 1 and activitydata[tDelayed][i] == 0:
						KLM010 +=1
					elif activitydata[t][j] == 0 and activitydata[tDelayed][j]  == 0 and activitydata[tDelayed][i] == 1:
						KLM001 +=1
					elif activitydata[t][j] == 1 and activitydata[tDelayed][j] == 1 and activitydata[tDelayed][i] == 0:
						KLM110 +=1
					elif activitydata[t][j] == 1 and activitydata[tDelayed][j] == 0 and activitydata[tDelayed][i] == 1:
						KLM101 +=1
					elif activitydata[t][j] == 0 and activitydata[tDelayed][j] == 1 and activitydata[tDelayed][i] == 1:
						KLM011 +=1
					elif activitydata[t][j] == 1 and activitydata[tDelayed][j] == 1 and activitydata[tDelayed][i] == 1:
						KLM111 +=1

					#all the permutations of the conditional probability term components,lm,of TE formula	
					elif activitydata[tDelayed][j] == 0 and activitydata[tDelayed][i] == 0:
						LM00 +=1
					elif activitydata[tDelayed][j] == 1 and activitydata[tDelayed][i] == 0:
						LM10 +=1
					elif activitydata[tDelayed][j] == 0 and activitydata[tDelayed][i] == 1:
						LM01 +=1
					elif activitydata[tDelayed][j] == 1 and activitydata[tDelayed][i] == 1:
						LM11 +=1
					#all the permutations of the conditional probability term components,kl,of TE formula	
					elif activitydata[t][j] == 0 and activitydata[tDelayed][j] ==  0:
						KL00 +=1
					elif activitydata[t][j] == 1 and activitydata[tDelayed][j] == 0:
						KL10 +=1
					elif activitydata[t][j] == 0 and activitydata[tDelayed][j] == 1:
						KL01 +=1
					elif activitydata[t][j] == 1 and activitydata[tDelayed][j] == 1:
						KL11 +=1
					#all the permutations of the conditional probability term components,l,of TE formula		
					elif activitydata[tDelayed][j] == 0:
						L0 +=1
					elif activitydata[tDelayed][j] == 1:
						L1 +=1
				#all the permutations of the conditional probability term components,klm,of TE formula	
				probabilities[h][i][j]['000'] = KLM000/timesteps
				probabilities[h][i][j]['100'] = KLM100/timesteps
				probabilities[h][i][j]['010'] = KLM010/timesteps
				probabilities[h][i][j]['001'] = KLM001/timesteps
				probabilities[h][i][j]['110'] = KLM110/timesteps
				probabilities[h][i][j]['101'] = KLM101/timesteps
				probabilities[h][i][j]['011'] = KLM011/timesteps
				probabilities[h][i][j]['111'] = KLM111/timesteps
				#all the permutations of the conditional probability term components,lm,of TE formula	
				probabilities[h][i][j]['lm00'] = LM00/timesteps
				probabilities[h][i][j]['lm10'] = LM10/timesteps
				probabilities[h][i][j]['lm01'] = LM01/timesteps
				probabilities[h][i][j]['lm11'] = LM11/timesteps
				#all the permutations of the conditional probability term components,kl,of TE formula	
				probabilities[h][i][j]['kl00'] = KL00/timesteps
				probabilities[h][i][j]['kl10'] = KL10/timesteps
				probabilities[h][i][j]['kl01'] = KL01/timesteps
				probabilities[h][i][j]['kl11'] = KL11/timesteps
				#all the permutations of the conditional probability term components,l,of TE formula		
				probabilities[h][i][j]['l0'] = L0/timesteps
				probabilities[h][i][j]['l1'] = L1/timesteps

	return probabilities			
'''
					 
def TE(G, probabilities, timedelay):
	TE = {}
	for h in range(1,timedelay):
		TE[h] = {}
		for i,nbrs1 in G.adjacency_iter():
			TE[h][i] = {}
			for j,nbrs2 in G.adjacency_iter():
				#upp = probabilities[h]['000']/probabilities[h]['lm00']
				#low = probabilities[h]['kl00']/probabilities[h]['l0']
				TEadd = 0
				if probabilities[h][i][j]['000'] > 0:
					TEadd = TEadd + probabilities[h][i][j]['000']*log((probabilities[h][i][j]['000']/probabilities[h][i][j]['lm00'])/(probabilities[h][i][j]['kl00']/probabilities[h][i][j]['l0']),2)									
				if probabilities[h][i][j]['100'] > 0:
					TEadd = TEadd + probabilities[h][i][j]['100']*log((probabilities[h][i][j]['100']/probabilities[h][i][j]['lm00'])/(probabilities[h][i][j]['kl10']/probabilities[h][i][j]['l0']),2)									
				if probabilities[h][i][j]['010'] > 0:
					TEadd = TEadd + probabilities[h][i][j]['010']*log((probabilities[h][i][j]['010']/probabilities[h][i][j]['lm10'])/(probabilities[h][i][j]['kl01']/probabilities[h][i][j]['l1']),2)				
				if probabilities[h][i][j]['001'] > 0:
					TEadd = TEadd + probabilities[h][i][j]['001']*log((probabilities[h][i][j]['001']/probabilities[h][i][j]['lm01'])/(probabilities[h][i][j]['kl00']/probabilities[h][i][j]['l0']),2)
				if probabilities[h][i][j]['110'] > 0:
					TEadd = TEadd + probabilities[h][i][j]['110']*log((probabilities[h][i][j]['110']/probabilities[h][i][j]['lm10'])/(probabilities[h][i][j]['kl11']/probabilities[h][i][j]['l1']),2)
				if probabilities[h][i][j]['101'] > 0:
					TEadd = TEadd + probabilities[h][i][j]['101']*log((probabilities[h][i][j]['101']/probabilities[h][i][j]['lm01'])/(probabilities[h][i][j]['kl10']/probabilities[h][i][j]['l0']),2)
				if probabilities[h][i][j]['011'] > 0:
					TEadd = TEadd + probabilities[h][i][j]['011']*log((probabilities[h][i][j]['011']/probabilities[h][i][j]['lm11'])/(probabilities[h][i][j]['kl01']/probabilities[h][i][j]['l1']),2)
				if probabilities[h][i][j]['111'] > 0:
					TEadd = TEadd + probabilities[h][i][j]['111']*log((probabilities[h][i][j]['111']/probabilities[h][i][j]['lm11'])/(probabilities[h][i][j]['kl11']/probabilities[h][i][j]['l1']),2)
				
				TE[h][i][j] = TEadd
	return TE

def dTE(G, TE, timedelay):
	dTE = {}
	for h in range(1, timedelay):
		dTE[h] = {}
		for i,nbrs1 in G.adjacency_iter():
			dTE[h][i] = {}
			for j,nbrs2 in G.adjacency_iter():
				val = TE[h][i][j] + TE[h][j][i]
				if val != 0:
					dTE[h][i][j] = TE[h][i][j] / (TE[h][i][j] + TE[h][j][i])
	return dTE

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
					if G.node[i]['degreeClass'] == 'High':
						dTE_high += dTE[h][i][j] 
						highCount += 1
					else:
						dTE_low += dTE[h][i][j]  
						lowCount += 1
		average_dTE_high = float(dTE_high)/float(highCount)
		average_dTE_low = float(dTE_low)/float(lowCount)
		HL[h] = average_dTE_high - average_dTE_low
	return HL