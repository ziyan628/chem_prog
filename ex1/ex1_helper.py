import numpy as np
from collections import deque
from sage.all import *


"""
helper function, contains all the helper for the main program in ex1_main.py
"""

a = 0
b = -1

"""
Function to generate a huckel matrix for linear and cyclic polyene.
Takes argument n for number of molecules.
Optionl argument for cyclic polyene
"""

def get_huckel_matrix_polyene(n, cyclic = False):
    matrix = np.zeros((n,n))
    np.fill_diagonal(matrix, a)
    np.fill_diagonal(matrix[1:], b)
    np.fill_diagonal(matrix[:,1:], b)
    if cyclic:
        matrix[0][n-1] = b
        matrix[n-1][0] = b
    return matrix

"""
get_eval function takes a matrix and returns its eigenvalues
eig_values are rounded to 10 d.p. in order to aviod errors in comparing floats
"""

def get_evals(matrix):
    eig_values = np.linalg.eig(matrix)[0]
    return np.round(eig_values,10)

"""
degen_solver takes a list of eigenvalues and prints each eigenvalues, along with its degeneracy
"""

def degen_solver(evals):
    degen=dict()
    for eval in np.unique(evals):
        degen[eval] = int(np.count_nonzero(evals == eval))
    for key in degen:
        print("Eigenvalue, ",round(key,2),", has a degeneracy of,", degen[key])
    return 



"""generate platonic solids, where i is the key of the platonic solid to be generated
this uses the sagemath library
"""

def get_huckel_matrix_platonic(i):
   
    d = {'a':graphs.TetrahedralGraph, 'b':graphs.HexahedralGraph, 
        'c':graphs.OctahedralGraph, 'd':graphs.DodecahedralGraph, 'e':graphs.IcosahedralGraph}
   
    return d[i]().adjacency_matrix()*-1
 


"""
generate_cyclic_face takes in matrix and a list of carbon positions in cyclic order 
updates the huckel matrix accordingly
"""

def generate_cyclic_face(matrix, carbon_list):
    carbons_to_the_left = deque(carbon_list)
    carbons_to_the_left.rotate(1)

    carbons_to_the_right = deque(carbon_list)
    carbons_to_the_right.rotate(-1)

    for carbon in carbon_list:

        matrix[carbon][carbons_to_the_left.popleft()] = b
        matrix[carbon][carbons_to_the_right.popleft()] = b
    return matrix

"""
make_a_bucky_ball makes a list of list called segments (i used a 2d "orange peel" representation to label the carbons)
each "segment" is a list of carbon positions in cyclic order
the segements is then iterated over the generate_cyclic_face function to return a huckle matrix for buckminsterfullerene
"""

def make_a_bucky_ball():
    segments = [[0,1,2,3,4],
               [0,1,8,7,6,5],
               [1,2,11,10,9,8],
               [2,3,14,13,12,11],
               [3,4,17,16,15,14],
               [4,0,5,19,18,17],
               [5,6,20,39,19],
               [6,7,23,22,21,20],
               [7,8,9,24,23],
               [9,10,27,26,25,24],
               [10,11,12,28,27],
               [12,13,31,30,29,28],
               [13,14,15,32,31],
               [15,16,35,34,33,32],
               [16,17,18,36,35],
               [18,19,39,38,37,36],
               [38,39,20,21,40,54],
               [21,22,42,41,40],
               [22,23,24,25,43,42],
               [25,26,45,44,43],
               [26,27,28,29,46,45],
               [29,30,48,47,46],
               [30,31,32,33,49,48],
               [33,34,51,50,49],
               [34,35,36,37,52,51],
               [37,38,54,53,52],
               [53,54,40,41,55,59],
               [41,42,43,44,56,55],
               [44,45,46,47,57,56],
               [47,48,49,50,58,57],
               [50,51,52,53,59,58],
               [55,56,57,58,59]         
                ]
    matrix = np.zeros((60,60))
    np.fill_diagonal(matrix, a)
    for segment in segments:
        matrix = generate_cyclic_face(matrix,segment)
    return matrix