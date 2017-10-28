#assign neurodata type to nodes
def assign_neuro_type():
	colnames = ['NAME', 'GROUP', 'TYPE']
	data = pandas.read_csv('data/neurogroup.csv', names=colnames)
	names = data.NAME.tolist()
	groups = data.GROUP.tolist()
	types = data.TYPE.tolist()
	for n,nbrs in G.adjacency_iter():
		for i in range(len(names)):
			if G.node[n]['cell_name'] == names[i]:
				G.node[n]['cell_type'] = types[i]

#determine if a neuron is excitory or inhibitory
def exin():
	i = 1.0 
	for n,nbrs in G.adjacency_iter():
		NT_types = ['Ach', 'DA', 'GABA', '5-HT']
		if G.node[n]['neurotransmitters'] == NT_types[0]:
			G.node[n]['exin'] = 1
		elif G.node[n]['neurotransmitters'] == NT_types[1]:
			G.node[n]['exin'] = 1
		elif G.node[n]['neurotransmitters'] == NT_types[2]:
			#10% are labeled inhibitory
			G.node[n]['exin'] = -1
			i += 1.0
		else:
			#63% of nodes are unassigned
			G.node[n]['exin'] = 1
	ratio_inhibitory = i / 270.0		
	return ratio_inhibitory

#initialise all hub nodes to have the parameter activity
def init_activity_perimeter():
	init_active_nodes = 0
	for n,nbrs in G.adjacency_iter():
		#node is inactive when degree of node is smaller than 10
		if len(nbrs) < 10:
			G.node[n]['activity'] = 0
		#node is set to initialise as active when the degree of node is greater or equals to 10	
		else:
			G.node[n]['activity'] = 100
			init_active_nodes += 1
		#calculate the percentage of active nodes	
		percentage_init_active = init_active_nodes / G.number_of_nodes()
		print(percentage_init_active)

#initialise all perimeter nodes to have the parameter activity
def init_activity_perimeter():
	init_active_nodes = 0
	for n,nbrs in G.adjacency_iter():
		#node is inactive when degree of node is smaller than 10
		if len(nbrs) > 10:
			G.node[n]['activity'] = 0
		#node is set to initialise as active when the degree of node is greater or equals to 10	
		else:
			G.node[n]['activity'] = 100
			init_active_nodes += 1
		#calculate the percentage of active nodes	
		percentage_init_active = init_active_nodes / G.number_of_nodes()
		print(percentage_init_active)		

#initialise all nodes to have the parameter activity
def init_activity_random():
	init_active_nodes = 0
	for n,nbrs in G.adjacency_iter():
		#randomly activate roughly 30% of nodes
		if random.random() > 0.20:
			G.node[n]['activity'] = 0
		else:
			G.node[n]['activity'] = 100
			init_active_nodes += 1
		#calculate the percentage of active nodes	
	percentage_init_active = init_active_nodes / G.number_of_nodes()
	print(percentage_init_active)

#initialize refractory period, all 
def init_refractory():
	init_active_nodes = 0
	for n,nbrs in G.adjacency_iter():
		G.node[n]['refractory'] = 0

#pull function to get the current activity of nodes used to visualize color of nodes in graph
def get_activity():
	activity_array = {}
	for n,nbrs in G.adjacency_iter():
		activity_array[n] = G.node[n]['activity']
	return activity_array

#create an array of the degree of each node which can be used in visualization for node size
def node_size_map():
	size_array = [0] * G.number_of_nodes()
	i = 0
	for n,nbrs in G.adjacency_iter():
		size_array[i] = G.degree(n) * 5
		i += 1
	return size_array

def normalize_synapse_weight():
	#find maximum weights for each types of synapses
	max_e_weight = 1
	max_c_weight = 1

	for n,nbrs in G.adjacency_iter():
		for nbr,eattr in nbrs.items():
			for attr, data in eattr.items():
				if data['synapse_type'] == 'E':
					if data['weight'] > max_e_weight:
						max_e_weight = data['weight']
				if data['synapse_type'] == 'C':	
					if data['weight'] > max_c_weight:
						max_c_weight = data['weight']
	#normalize for each synapse
	for n,nbrs in G.adjacency_iter():
		for nbr,eattr in nbrs.items():
			for attr, data in eattr.items():
				if data['synapse_type'] == 'E':
					data['normal_weight'] = data['weight'] / max_e_weight
				if data['synapse_type'] == 'C':
					data['normal_weight'] = data['weight'] / max_c_weight	

#interate over all nodes to propogate neural activity
def single_time_step(node_sizes):
	integral= [0] * G.number_of_nodes()
	m = 0
	for n,nbrs in G.adjacency_iter():
		#check if the node is active

		#decay of activity of activated neuron in 2 time steps
		if G.node[n]['activity'] > 0:
			#an activated node will be activated for 2 timesteps
			G.node[n]['activity'] -= 50 

			current_activity = G.node[n]['activity']
			#set refractory period if activity of the node just ended
			if current_activity == 0:
				#refactory period takes 3 time steps to end
				G.node[n]['refractory'] = 1

		#if the node is in the refactory period reduce its count	
		elif G.node[n]['refractory'] > 0:
			G.node[n]['refractory'] -= 1
		
		#else determine the sum of all activities of its neighbouring nodes and decide if the integral is sufficient for firing	
		
		else:
			#initialize integral
			for nbr,eattr in nbrs.items():
				for attr, data in eattr.items():
					#'E' for electrical synapse
					if data['synapse_type'] == 'E':
						#summing the activity input into a node and store integral into a list
						integral[m] +=  G.node[nbr]['exin'] * G.node[nbr]['activity'] * data['normal_weight']
			#this threshold activation limit is chosen based on the proportion of neuron action potential			
			if integral[m] > 2:
				G.node[n]['activity'] = 100
				activationCount[n] += 1
		#for tracking the integral list		
		m += 1

		return activationCount

	"""
	#print current activities and integral
	print get_activity()
	print integral
	"""

#main function for time iteration that contain all smaller functions
def time_itr(time,iteration):
	assign_neuro_type()
	exin()
	init_activity_random()
	init_refractory()
	normalize_synapse_weight()
	node_sizes = node_size_map()

	for i in range(time):
		#figure perimeter set up
		plt.figure(figsize=(12,12))
		#pos=nx.spring_layout(G,iterations=100,scale=2.0)
		#n_colors=range(279)
		#e_colors=range(3225)
		#draw graphs so propogation can be seen in real time
		#nx.draw_spectral(G)
		nx.draw(G,pos, node_color=get_activity(), node_size=node_sizes, width=1, style='dotted', arrows=False, cmap=plt.cm.Blues)
		plt.savefig("img/test/" + str(iteration) + "step_n1_"  + str(i) + ".png")

		#nx.draw_circular(G, node_color=get_activity(), node_size=node_sizes, width=1, style='dotted', arrows=False, cmap=plt.cm.Blues)

		#plt.savefig("img/step_cr3_" + str(i) + ".png")
		#plt.show()
		single_time_step(node_sizes,iteration)
		activitydata[iteration][time] = get_activity()

		if sum(get_activity()) == 0
			dieDownTime[iteration] = i
			break

#importing the wormNet data from graphml file