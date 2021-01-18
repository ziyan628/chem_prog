import numpy as np
from collections import deque
import ex1_helper as helper

"""
Specifications:
For an arbitrary molecule, specified by its connectivity between adjacent atoms, 
calculate and print the Hu ̈ckel energies and degeneracies of its π-system. 

Tasks:
1. A linear polyene with n carbons;
2. A cyclic polyene with n carbons;
3. The sp2-hybridized Platonic solids (i.e. 3D structures with faces of a single regular polygon, with carbons at each vertex and only 3 edges per vertex): tetrahedron, cube, and dodecahe- dron. (If you wish to suspend your chemical disbelief you might also try the octahedron and icosahedron );
4. (Optional) Buckminsterfullerene.
"""

tasks = ["Linear polyene", "Cyclic polyene", "Platonic solids", "Buckminsterfullerene"]

print("Welcome to this Ziyan's Huckel solver")

print("Select method:")

for method in tasks:
    print("FOR ", method, ", SELECT", tasks.index(method))

method = int(input())

if method < 3:
    print("Input number of molecules")
    n = int(input())
    if method == 0:
        matrix =  helper.get_huckel_matrix_polyene(n)
    elif method == 1:
        matrix =  helper.get_huckel_matrix_polyene(n, cyclic = False)
    elif method == 2:
        print("Input number of edges per vertex")
        j = int(input())
        matrix = helper.get_huckel_matrix_platonic(n,j)

if method == 3:
    matrix = helper.make_a_bucky_ball()

helper.degen_solver(helper.get_evals(matrix))



