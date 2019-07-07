"""rips.py by Jernej Puc"""

##############################################################################################
# MODULES

import numpy as np
from time import perf_counter


##############################################################################################
# MAIN

def cliques(VG, EG):
    """Finds all cliques in a graph, given as lists of vertices and edges."""
    
    # Redefine graph as a hash-map of neighbours
    G = {}

    for v in VG:
        G[v] = []
    
    for e in EG:
        G[e[0]].append(e[1])
        G[e[1]].append(e[0])
    
    # Call modified Bron-Kerbosch algorithm for all possible degrees of cliques
    C = {}

    for i in range(len(VG)):
        Ci = bronkerbosch([],VG,[],G,i)
        
        if Ci:
            C[i] = Ci
    
    return C


def VR(S, epsilon):
    """Returns the simplices in the Vietoris-Rips complex in their corresponding dimension."""
    
    # Redefine points as a hash-map of vertices
    V = {}

    for i in range(len(S)):
        V[i] = np.array(S[i])
    
    # Find valid edges between pairs of vertices
    EG = []

    for i in range(len(S)-1):
        for j in range(i+1,len(S)):
            if np.linalg.norm(V[i] - V[j]) <= epsilon:
                EG.append((i,j))
    
    # Get associated cliques
    return cliques(sorted(V.keys()),EG)


##############################################################################################
# AUXILIARY

def bronkerbosch(R,P_,X_,G,s):
    """Meant to find maximal cliques in a graph. Modified to find cliques of specific degree."""
    
    # Found clique of degree s
    if len(R)-1 == s:
        return [tuple(R)]
    
    # Found maximal clique, but not of degree s
    elif not P_ and not X_:
        return []
    
    # Recursive step
    else:
        Cs = []
        P = P_.copy()
        X = X_.copy()

        for v in P_:
            Cs += bronkerbosch(R+[v],
                               [p for p in P if p in G[v]],
                               [x for x in X if x in G[v]],
                               G,s)
            P.remove(v)
            X.append(v)
        
        return Cs


##############################################################################################
# Test

if __name__ == '__main__':

    # VR on a reference example
    S = [(0,0), (1,1), (2,3), (-1,2), (3,-1), (4,2)]
    epsilon = 3
    print(VR(S,epsilon),'\n')
    # Prints out:
    # {0: [(0,), (1,), (2,), (3,), (4,), (5,)], 1: [(0, 1), (0, 3), (1, 2), (1, 3), (1, 4), (2, 5)], 2: [(0, 1, 3)]}

    # Example 1 (cliques)
    VG = [0,1,2,3]
    EG = [(0,1),(0,3),(1,2),(1,3),(2,3)]
    print(cliques(VG,EG),'\n')
    # Prints out:
    # {0: [(0,), (1,), (2,), (3,)], 1: [(0, 1), (0, 3), (1, 2), (1, 3), (2, 3)], 2: [(0, 1, 3), (1, 2, 3)]}
    
    # Example 2 (cliques)
    VG = [0,1,2,3]
    EG = [(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]
    print(cliques(VG,EG),'\n')
    # Prints out:
    # {0: [(0,), (1,), (2,), (3,)], 1: [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)], 2: [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)], 3: [(0, 1, 2, 3)]}
    
    # Example 3 (cliques)
    VG = [0,1,2,3,4,5,6,7]
    EG = [(0,1),(0,5),(1,2),(1,6),(2,3),(2,6),(2,7),(3,4),(3,6),(3,7),(4,5),(4,7),(6,7)]
    print(cliques(VG,EG),'\n')
    # Prints out:
    # {0: [(0,), (1,), (2,), (3,), (4,), (5,), (6,), (7,)], 1: [(0, 1), (0, 5), (1, 2), (1, 6), (2, 3), (2, 6), (2, 7), (3, 4), (3, 6), (3, 7), (4, 5), (4, 7), (6, 7)], 2: [(1, 2, 6), (2, 3, 6), (2, 3, 7), (2, 6, 7), (3, 4, 7), (3, 6, 7)], 3: [(2, 3, 6, 7)]}

    # Example 4 (VR)
    S = [(-2,0),(0,-1),(2,0),(0,1)]
    epsilon = np.sqrt(5)
    print(VR(S,epsilon),'\n')
    # Prints out:
    # {0: [(0,), (1,), (2,), (3,)], 1: [(0, 1), (0, 3), (1, 2), (1, 3), (2, 3)], 2: [(0, 1, 3), (1, 2, 3)]}
    
    # Example 5 (VR)
    S = [(-1,0),(0,-1),(1,0),(0,1)]
    epsilon = 2
    print(VR(S,epsilon),'\n')
    # Prints out:
    # {0: [(0,), (1,), (2,), (3,)], 1: [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)], 2: [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)], 3: [(0, 1, 2, 3)]}
    
    # Example 6 (VR)
    S = [(-1,2),(-2,-0.5),(-1,-3),(1,-3),(2,-0.5),(1,2),(-1,-1),(1,-1)]
    epsilon = np.sqrt(8)
    print(VR(S,epsilon),'\n')
    # Prints out:
    # {0: [(0,), (1,), (2,), (3,), (4,), (5,), (6,), (7,)], 1: [(0, 1), (0, 5), (1, 2), (1, 6), (2, 3), (2, 6), (2, 7), (3, 4), (3, 6), (3, 7), (4, 5), (4, 7), (6, 7)], 2: [(1, 2, 6), (2, 3, 6), (2, 3, 7), (2, 6, 7), (3, 4, 7), (3, 6, 7)], 3: [(2, 3, 6, 7)]}

    # Cliques of a fully-connected graph
    n = 10  # < 1 sec
    # n = 11  # < 3 sec
    # n = 13  # ~ 50 sec
    VG = [i for i in range(n)]
    EG = []

    for i in range(n-1):
        for j in range(1,n):
            EG.append((i,j))

    ctr = perf_counter()
    _ = cliques(VG,EG)
    print(perf_counter() - ctr)
