
## CODE FOR TRANSFORMING A EXPRESSION MATRIX TO A NCOL FORMAT NETWORK THROUGH A CUTOFF
## Developed by Núria Olvera Ocaña
## Institut d'Investigacions Biomèdiques August Pi i Sunyer (IDIBAPS) and Barcelona Supercomputing Center (BSC)
## Inflammation and repair in respiratory diseases group at IDIBAPS and Computational Biology group at BSC
## Contact: olvera@clinic.cat

# load packages
import sys
import argparse
import fileinput
import os
cwd = os.getcwd()
file2= open(cwd +'/'+ sys.argv[1], 'r')


# function
def ncol_format(file, cutoff=0):
    ''' 
    Parameters: file containing a individual x individual expression matrix. Cutoff to decide whether 
    two nodes should be connected or not.
    Returns: a dictionary of lists where the keys represent nodes and the lists contain a node and 0/1 depending 
    on whether they are connected or not.
    '''
    Dict= {}
    i= 0
    j= 0
    names= list()    
    for line in file:
        j= j + 1   
        Line= line.rstrip().split(' ')
        Dict[Line[0]]= []
        for element in Line:
            i= i + 1
            if line.startswith(' '):
                names= Line               
            else:                
                if element != Line[0] and j != i and names[i-1] not in Dict.keys():                   
                   
                    if float(element) > float(cutoff):
                        Dict[Line[0]].append([names[i-1], 1])
                    elif float(element) < float(cutoff):
                        Dict[Line[0]].append([names[i-1],0])
                
        i=0
    
                                              
    return Dict        

# call function
fun= ncol_format(file= file2, cutoff= sys.argv[2])

