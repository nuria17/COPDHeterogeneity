## CODE TO CREATE THE RANDOM GRAPHS AND RETRIEVE THEIR OPTIMAL COMMUNITIES THROUGH MOLTI SOFTWARE
## Developed by Núria Olvera Ocaña
## Institut d'Investigacions Biomèdiques August Pi i Sunyer (IDIBAPS) and Barcelona Supercomputing Center (BSC)
## Inflammation and repair in respiratory diseases group at IDIBAPS and Computational Biology group at BSC
## Contact: olvera@clinic.cat



####LOADING PACKAGES####
import networkx as nx
import itertools
import subprocess
import matplotlib.pyplot as plt
import os
import statistics
import numpy as np
from scipy import stats
from statsmodels.distributions.empirical_distribution import ECDF
import pandas as pd
from collections import Counter
from scipy.spatial import distance
import re
import copy

####FUNCTIONS####
def ncol_format(corr_matrix, random_graph):
    
    '''
    Returns a dictionary of lists where the keys represent nodes and the lists 
    represent a node and the weigth (in this case 0 or 1). The file1 input corresponds
    to a correlation matrix and the random_graph is a dictionary of dictionaries, representing a graph.
    '''
    
    file= open(corr_matrix, 'r')    
    Dict= {}
    i= 0
    j= 0
    names= list()   
    for line in file:
        Line= line.rstrip().split(' ')          
        if line.startswith(' '): #First line corresponds to nodes' names
            names= Line
            names.pop(0)
        else:
            j= j + 1 #rows
            Dict[names[j-1]]= []     
            for element in Line:                
                i= i + 1 #columns                             
                if element != Line[0] and j != i and names[(i-2)] not in Dict.keys(): #not adding repeated edges in the Dict since the output will be an undirected graph             
                    if j-1 in random_graph.keys(): 
                        if (i-2) in random_graph[j-1]:
                            Dict[names[j-1]].append([names[i-2], 1])
                        elif (i-2) not in random_graph[j-1]:
                            Dict[names[j-1]].append([names[i-2],0])                    
            i=0    
                                                
    return Dict  

def sequence(seq):
    '''Creates a degree sequence from a degree distribution'''
    s=open(seq, 'r') #Degree distribution
    sequence=[]
    for edges in s:
        Edges= edges.rstrip().split()
        if sequence == []:
            sequence= [i for i in itertools.repeat(int(Edges[0]), int(Edges[1]))]
        else:
            sequence= sequence + [i for i in itertools.repeat(int(Edges[0]), int(Edges[1]))]
    return sequence

###CALLING FUNCTIONS###
Degree_file= 'Data/Degree_file'
Corr_matrix= 'Data/Corr_matrix'
path_random_graph= 'Data/Random_graphs/'


#Loop to create 10000 random graphs with the same degree sequence as the original network and finding their optimal communities with MolTi software.
for i in range(0,10000):        
    deg_seq= sequence(Degree_file)
    G=nx.random_degree_sequence_graph(deg_seq, tries= 10) #Return a simple random graph with the given degree sequence.
    Adj= nx.to_dict_of_dicts(G, edge_data=1) #Return adjacency representation of graph as a dictionary of dictionaries.
    fun= ncol_format(corr_matrix= Corr_matrix, random_graph= Adj) 
    #Creating the files in the suitable NCOL format to give the input to MolTi software.
    t=''
    for key in fun:
        for key2 in fun[key]:
            if t == '':
                t= str(key) + '\t' + str(key2[0]) + '\t' + str(key2[1]) + '\n' 
            else:
                t+= str(key) + '\t' + str(key2[0]) + '\t' + str(key2[1]) + '\n'
    f3= open(path_random_graph + 'Random_graph_micro_'+ str(i), 'w')
    f3.write(t.rstrip())
    f3.close()
    # Calling MolTi software through subprocess module
    process2 = subprocess.Popen(['molti-console','-p',str(1),'-o','Data/Random_clusters_' + str(i) , 'Data/Random_graph_file_' + str(i)],
                    stdout=subprocess.PIPE, 
                    stderr=subprocess.PIPE)
    stdout, stderr = process2.communicate()
    print(stdout.rstrip())