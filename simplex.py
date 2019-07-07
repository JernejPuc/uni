"""simplex.py v19.05.28 by Jernej Puc"""


#####################################################################################################
# Modules

import numpy as np


#####################################################################################################
# Main

def simplex(c, A, b, J, K):
    """
    Implementation of the revised simplex method for linear programming.
    For reference: https://en.wikipedia.org/wiki/Revised_simplex_method

    Primary LP definition:
        min c'*x
        A*x == b
        x >= 0

    Input:
        c: vector of objective coefficients
        A: matrix of constraining coefficients
        b: vector of constraints
        J: vector of base variable indices
        K: vector of remaining indices

    Output:
        x: vector of non-negative variables, which solve the primary LP
    """

    x = np.zeros_like(c)

    while True:
        # Current base variables
        Aj = A[:,J]
        xj = np.linalg.solve(Aj, b)
        
        # Associated objective coefficients
        cj = c[J,0]

        # Current objective value
        f = np.dot(cj, xj)
        
        # Dual vector y
        y = np.linalg.solve(Aj.T, cj)
        
        # If ccT >= 0, the optimum has been reached
        ccT = c[K,0].T - np.dot(y.T, A[:,K])

        if np.all(ccT >= 0.):
            x[J,0] = xj.T
            break
        
        # Otherwise, an appropriate variable (Bland's rule) is chosen to enter the base
        s = np.min(K[ccT < 0.])
        
        # If aa <= 0, the problem is unbounded
        aa = np.linalg.solve(Aj, A[:,s])

        if np.all(aa <= 0.):
            raise Exception('LP is unbounded.')
        
        # Otherwise, an appropriate variable (Bland's rule) is chosen to exit the base
        aa = np.where(aa == 0., 1e-12, aa)
        v = np.divide(xj.T, aa)[0]
        r = np.min(J[v == np.min(v[v > 0.])])
        
        # Update of base variable indices
        J[J == r] = s
        K[K == s] = r
    
    return f[0], x, J, K


def simplex2(c, A, b):
    """
    Two-phase execution of the simplex method, asserting feasibility.
    """

    # Initial base variable indices
    m,n = A.shape
    J = np.arange(n,n+m)
    K = np.arange(n)

    # Phase 1
    c_ = np.vstack((np.zeros([n,1]), np.ones([m,1])))
    A_ = np.hstack((A, np.eye(m)))

    f1, _, J1, K1 = simplex(c_, A_, b, J, K)

    # Assert feasibility
    if not np.isclose(f1, 0):
        raise Exception('LP is infeasible.')

    # Feasible base
    J = J1[J1 < n]
    K = K1[K1 < n]

    # If a feasible solution includes base variables set to zero,
    # extended variables can remain and cause an issue
    mj = len(J)
    
    if mj != m:
        J = np.concatenate((J, K[mj-m:]))
        K = K[:mj-m]
    
    # Phase 2
    f2, x2, _, _ = simplex(c, A, b, J, K)

    # Final solution
    return f2, x2


#####################################################################################################
# Test

if __name__ == '__main__':

    # NOTE: c,A,b should already include slack/surplus variables
    
    # Problem 1
    c = np.array([[-3.],
                  [-1.],
                  [-3.],
                  [0.],
                  [0.],
                  [0.]])

    A = np.array([[2., 1., 1., 1., 0., 0.],
                  [1., 2., 3., 0., 1., 0.],
                  [2., 2., 1., 0., 0., 1.]])

    b = np.array([[2.],
                  [5.],
                  [6.]])

    print('\nAttempting problem 1:')
    
    f,x = simplex2(c,A,b)

    print('f:', round(f,6))
    print('x:\n', np.around(x,6), '\n')


    # Problem 2
    c = np.array([[-1.],
                  [-1.],
                  [0.],
                  [0.],
                  [0.]])

    A = np.array([[1.,0.,2.,1.,0.],
                  [0.,1.,-1.,0.,1.],
                  [1.,1.,1.,0.,0.]])

    b = np.array([[1.],
                  [1.],
                  [2.]])

    print('Attempting problem 2:')
    
    f,x = simplex2(c,A,b)

    print('f:', round(f,6))
    print('x:\n', np.around(x,6), '\n')
    
