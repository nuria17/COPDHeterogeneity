## CODE FOR CALCULATING THE PROBABILITY OF FINDING A COMMUNITY OF THE SAME SIZE, WITH THE NODES AND WITH THE SAME PAIR OF NODES AS THE ORIGINAL NETWORK BY CHANCE.
## Developed by Núria Olvera Ocaña
## Institut d'Investigacions Biomèdiques August Pi i Sunyer (IDIBAPS) and Barcelona Supercomputing Center (BSC)
## Inflammation and repair in respiratory diseases group at IDIBAPS and Computational Biology group at BSC
## Contact: olvera@clinic.cat


### LOAD PACKAGES ###
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

### FUNCTIONS ###

def size_clusters(f):

    '''
     Parameters: a file returned by the MolTi software where the clusters are represented.
     Returns: a dictionary where keys represent clusters and values represent the cluster size.
    '''
    
    Dict_size={}
    file= open(f, 'r')
    counter= 0
    name= ''
    for line in file:
        Line= line.rstrip()     
        if Line.startswith('ClusterID') and name in Dict_size.keys():
            Dict_size[name]= counter -1
            counter= 0
            name= Line
            Dict_size[Line]= 0
        elif Line.startswith('ClusterID') and Dict_size == {}:
            name= Line
            Dict_size[name]= 0
        elif Line.startswith('#'):
            pass
        else:             
            counter= counter + 1
    Dict_size[name]= counter -1
    return Dict_size


def SizeProbability(path_random_graphs):
    
    '''
    Parameters: directory where the random graphs are located.
    Returns: a list with the probabilities of community size in the random networks.
    
    '''
    
    Files= [f for f in os.listdir(path_random_graphs) if f.endswith('effectif.csv')]
    size_cluster=[]
    for f in Files:
        F= open(path_random_graphs + f, 'r')
        for line in F:
            if line.startswith('Cluster'):
                Line=line.rstrip().split('\t')
                size_cluster.append(int(Line[1]))    

    Size_cluster= np.array(size_cluster)
    s = pd.Series(Size_cluster)

    vector = (s.groupby(s).transform('count') / len(s)).values
    freq_dict = Counter(Size_cluster)

    prob_dict = dict((k,round(val/len(Size_cluster),7)) for k, val in freq_dict.items())
    prob_dict_sorted= dict(sorted(prob_dict.items()))
    return_list= [vector,Size_cluster,prob_dict_sorted]
    return return_list

def id_to_cluster(molti_file):
    
    '''
    Parameters: file with communities returned by MolTi software.
    Returns: dictionary where keys represent community name and values represent nodes
    
    '''
    
    File_molti= open(molti_file, 'r')
    llista_ids=[]
    Dict_clust={}
    Dict_zeros={}
    name=''
    ids=["1","10","101","102","103","104","106","107","108","109","110","111","112","113","114","115","116","117","118","119","12","120","122","13","14","15","16","17","2","20","200","201","202","204","206","207","208","209","21","22","24","25","27","28","29","3","30","301","302","304","305","306","308","309","31","310","312","313","314","315","317","318","32","320","321","322","323","325","326","33","34","35","37","38","39","4","40","41","42","43","44","45","47","48","49","5","50","51","52","53","54","55","56","57","58","59","6","60","61","62","63","64","65","66","67","68","71","73","74","75","78","79","8","80","81","82","83","84","85","86","87","9","R100","R103","R113","R119","R126","R20","R22","R27","R33","R34","R36","R38","R39","R44","R46","R48","R50","R52","R53","R56","R63","R66","R67","R77","R78","R80","R84","R86","R87","R93","R94","R99"]
    c=0
    for line in File_molti:
        Line= line.rstrip()
        if Line.startswith('ClusterID'):
            name=Line
            Dict_clust[name]=[]
        elif name != '' and Line != '':
            Dict_clust[name].append(Line)
        else:
            pass
    return Dict_clust

def CommunityProbability(original_network, random_network):
    
    '''
    Parameters: a directory where the random graphs are located and the file with the original communities returned by MolTi.
    Returns: dictionary where keys represent communities and values the probabilities of getting the same community in the random networks. 
    
    '''
    
    Dict_com={}
    for f in os.listdir(random_network):
        fun= id_to_cluster(molti_file= random_network + f)
        fun2= id_to_cluster(molti_file= original_network)
        for key,value in fun2.items():
            if key not in Dict_com.keys():
                if value in fun.values():                
                    Dict_com[key]=1
                else:
                    Dict_com[key]= 0
            else: 
                if value in fun.values():
                    Dict_com[key]+= 1
    return Dict_com

def pair_nodes(file):
    
    '''
    Parameters: file with communities returned by MolTi software.
    Returns: a dictionary of lists where the keys represent nodes and the lists represent all the nodes that are in the same community.
    '''
    
    Dict_p={}
    value2=[]
    fun=id_to_cluster(file)
    for key,value in fun.items():
        for val in value:
            Dict_p[val]=[]
            if value2 == []:
                value2= copy.deepcopy(value)
                value2.remove(val)
                Dict_p[val]= value2
                value3= value2       
            else:
                value4= copy.deepcopy(value3)
                if value4 != []:
                    value4.remove(val)                
                    Dict_p[val]=value4
            value3= Dict_p[val]
        value2=[]
            
    return Dict_p

def PairsProbability(path_random_graphs, original_network):
    
    '''
    Parameters: path to the random graphs directory and file with the communities returned by MolTi software.
    Returns: dictionary of lists where keys represent the communinity name and the list contain the probability of finding every 
    pair of the community in the random networks.
    
    '''
    
    Dict_n={}
    num_comm=0
    for f in os.listdir(path_random_graphs):
        if re.match('Random_clusters_all_\d{4}$', f):        
            fun= pair_nodes(original_network)
            fun2= pair_nodes(path_random_graphs + f)
            a= len(id_to_cluster(molti_file= '/home/nuria/Desktop/CV_filt/Backbones_molti/Random_graphs_all/'+ f).keys())
            num_comm= num_comm + a
            for key,value in fun.items():
                if value != []:
                    for val in value:
                        name= str(val) + '_' + str(key)
                        name_i= str(key) + '_' + str(val)
                        if name not in Dict_n.keys() and name_i not in Dict_n.keys():
                            if val in fun2[key] or key in fun2[val]:
                                Dict_n[name] = 1
                            else:
                                Dict_n[name] = 0
                        elif name in Dict_n.keys() or name_i in Dict_n.keys():
                            if val in fun2[key] or key in fun2[val]:
                                Dict_n[name] = Dict_n[name] + 1
    fun3= id_to_cluster(molti_file='/home/nuria/Desktop/CV_filt/Backbones_molti/Backbone_closure_ultra_all_good')
    Dict_pc={}
    for key, val in Dict_n.items():
        Key= key.split('_')
        p1=Key[0]
        p2=Key[1]
        freq= val/num_comm
        for key2,val2 in fun3.items():    
            if key2 in Dict_pc.keys():
                if p1 in fun3[key2] and p2 in fun3[key2]:
                    Dict_pc[key2].append(freq)
            else:
                if p1 in fun3[key2] and p2 in fun3[key2]:
                    Dict_pc[key2]= [freq]
                    break
    return Dict_pc


