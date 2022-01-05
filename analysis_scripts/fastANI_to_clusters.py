#! /usr/bin python


### script <fastANI output> <threshold as number between 0 - 100>
import sys
import matplotlib.pyplot as plt
import networkx as nx

network_file = sys.argv[1]
threshold = int(sys.argv[2])
g =  nx.Graph() # create new empty undirected network

cnt = 0

new_network = []

with open(network_file) as f: # open the file
	for line in f: # read file line by line
		toks = line.strip().split(",") # split the line according to sep
		node1 = toks[0]
		node2 = toks[1]
		weight = float(toks[2])

		## add the nodes
		g.add_node(node1) #even if node1 or node2 are already in the network it's fine
		g.add_node(node2)

		if weight >= threshold:
			g.add_edge(node1, node2, weight = weight) # add an edge only if the weight bigger than threshold
			new_network.append(','.join(line.split(' ')[:3]))
		cnt+=1
		#if cnt == 1000:
		#	break


out = open("clusters.csv", "w")
out.write("Name, Cluster\n")

with open("network_output_threshold_"+str(threshold)+".csv", "w") as output:
	output.write("node1,node2,weight\n")

with open("network_output_threshold_"+str(threshold)+".csv", "a") as output:
	for line in new_network:
		output.write(line+"\n")



cc = nx.connected_components(g)  ## get clusters from graph (connected components)

cluster_num = 0
for c in cc: # go over each CC
	c = list(c)
	cluster_num += 1
	for node in c:
		out.write(node + "," + str(cluster_num) + "\n")

out.close()
