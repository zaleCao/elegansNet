import networkx as nx
import matplotlib.pyplot as plt



#import worm connectome as a nx graph object



G=nx.Graph()

G.add_node(1)
G.add_nodes_from([2,3])

H=nx.path_graph(10)
G.add_nodes_from(H)

G.add_edge(1,2)
e=(2,3)
G.add_edge(*e)


G.number_of_nodes()

G.number_of_edges()

G.neighbors(2)

G.add_node(11, role='motor')

G.add_edge(10, 11, weight=4.7)
G.edge[10][11]
nx.shortest_path(G,1,3)

lollipop=nx.lollipop_graph(10,20)

plt.figure(figsize=(8,8))
nx.draw(lollipop)
plt.show()

G = nx.read_graphml("data/c.elegans.herm_pharynx_1.graphml")
