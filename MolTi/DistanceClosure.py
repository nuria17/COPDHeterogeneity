## CODE TO RETRIEVE THE BACKBONE EDGES BASED ON THE DISTANCE CLOSURE 
## Developed by Núria Olvera Ocaña
## Institut d'Investigacions Biomèdiques August Pi i Sunyer (IDIBAPS) and Barcelona Supercomputing Center (BSC)
## Inflammation and repair in respiratory diseases group at IDIBAPS and Computational Biology group at BSC
## Contact: olvera@clinic.cat


###LOADING PACKAGES###
from distanceclosure.utils import prox2dist, dist2prox
from distanceclosure.distance import pairwise_proximity
from distanceclosure.closure import transitive_closure
from distanceclosure.backbone import backbone
import numpy as np

###FUNCTIONS###
def from_corr_matrix_to_array(file):
    
    '''Parameters: File representing a correlation matrix.
       Returns: a numpy array of arrays representing the correlation matrix'''    

    file= open(file, 'r')
    Array= []
    for line in file:    
        if not line.startswith(' '):    
            Line= line.rstrip().split()   
            Array.append(Line[1:])
        X= np.array(Array).astype(np.float)
    file.close()
    return X

def map_index_to_ID(file):

    '''Parameters: File representing a correlation matrix.
       Returns: list of lists where each list contains connected nodes.'''

    llista=[]
    File= open(file, 'r')
    first_line=File.readlines()
    F= first_line[0].rstrip().split()
    for index in np.arange(len(rows)):    
        Index_r= rows[index]
        Index_c= cols[index]
        if [F[Index_c], F[Index_r]] not in llista:
            llista.append([F[Index_r], F[Index_c]])
    File.close()
    return llista

# First, the DEanalysis.r script needs to be run to get the individual x individual correlation matrices.
# Convert to Distance
corr_matrix= 'Data/Correlation_matrix.txt'
corr_matrix_array= from_corr_matrix_to_array(corr_matrix)
D=prox2dist(corr_matrix_array)


# Calculate transitive closure using the metric and ultra-metric (also known maximum flow) measure
Cm = transitive_closure(D, kind='metric')
Cu = transitive_closure(D, kind='ultrametric')

# Retrieve the backbone edges
Bm = backbone(D, Cm)
Bu = backbone(D, Cu)

# Backbone edges on Bm and Bu can be accessed, where edges with a 1 are metric.
# We used the backbone edges based on the ultra-metric transitive closure.
rows,cols = np.where(Bu==1)

# We get a two columns file where only the nodes that are connected through a backbone edge are showed.
t=""
for element in map_index_to_ID(corr_matrix):    
    if t== "":
        t= element[0] + ',' + element[1] + '\n'
    else:
        t += element[0] + ',' + element[1]+ '\n'

fd=open('Data/Backbone_edges.txt', 'w')
fd.write(t.rstrip())
fd.close()
