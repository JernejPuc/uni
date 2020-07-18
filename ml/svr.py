"""
Support vector regression
"""

import numpy as np
import cvxopt
import unittest
from unittest.mock import patch


# Disable verbosity
cvxopt.solvers.options['show_progress'] = False


class Linear:
    """Linear kernel."""

    def __init__(self):
        pass

    def __call__(self, A, B):
        return A.dot(B.T)


class Polynomial:
    """
    Polynomial kernel.
    """

    def __init__(self, M=2):
        self.M = M

    def __call__(self, A, B):
        return np.power(np.dot(A, B.T) + 1., self.M)


class RBF:
    """
    Radial basis function kernel.
    """

    def __init__(self, sigma=0.5):
        self.coeff = -1. / (2. * sigma**2)

    def __call__(self, A, B):
        AdotB = np.dot(A, B.T)
        a_shape = None
        b_shape = None

        if len(A.shape) > 1:
            a_shape = (AdotB.shape[0],1) if len(AdotB.shape) > 1 else AdotB.shape
        
        if len(B.shape) > 1:
            b_shape = (1,AdotB.shape[1]) if len(AdotB.shape) > 1 else AdotB.shape

        return np.exp(self.coeff * (np.sum(A*A, axis=-1).reshape(a_shape) - 2.*AdotB + np.sum(B*B, axis=-1).reshape(b_shape)))


class SVRModel:
    """
    Support vector regression -- fitted model.
    """

    def __init__(self, alpha, b, kernel, X, sub=0., den=1.):
        self.alpha = alpha
        self.b = b
        self.kernel = kernel
        self.X = X
        self.sub = sub
        self.den = den

    def predict(self, X):
        return np.dot(self.alpha[:,0]-self.alpha[:,1], self.kernel((X-self.sub)/self.den, self.X).T) + self.b
    
    def get_alpha(self):
        return self.alpha
    
    def get_b(self):
        return self.b


class SVR:
    """
    Support vector regression -- fitting.

    Matrices (and vectors) P, q, G, h, A, b defined as per
    https://cvxopt.org/userguide/coneprog.html#quadratic-programming

    Note: Smola and Scholkopf state the optimisation as a maximisation problem,
    while cvxopt.solvers.qp aims to minimise the objective function.
    P and q thus differ from the original formulation by a factor of -1.
    """

    def __init__(self, kernel=RBF(sigma=0.5), lambda_=1e-4, epsilon=0.1):
        self.kernel = kernel
        self.C = 1./lambda_
        self.epsilon = epsilon
    
    def fit(self, X, y, data_prep=None, tol_0=1e-4, tol_C=1e-3):
        n = y.size

        # Data preprocessing
        if data_prep is not None:
            if data_prep == 'normalise':
                sub = np.min(X, axis=0)
                den = np.max(X - sub, axis=0)
            elif data_prep == 'standardise':
                sub = np.mean(X, axis=0)
                den = np.std(X, axis=0)
            else:
                sub, den = 0., 1.
        else:
            sub, den = 0., 1.
            
        X = (X - sub) / den
        
        # Def. P, q, G, h, A, b
        Pij = np.array([[1., -1.],
                        [-1., 1.]])
        P = np.tile(Pij, (n,n)) * np.repeat(np.repeat(self.kernel(X,X), 2, axis=0), 2, axis=1)

        q = np.empty((2*n,1))
        q[0::2,0] = -1.*y + self.epsilon
        q[1::2,0] = y + self.epsilon

        G = np.vstack((-1.*np.eye(2*n), np.eye(2*n)))
        h = np.concatenate((np.zeros(2*n), self.C*np.ones(2*n))).reshape(4*n,1)

        A = np.empty((1,2*n))
        A[0,0::2] = np.ones(n)
        A[0,1::2] = -1.*np.ones(n)
        b = 0.

        # Convert objects from numpy.ndarray to cvxopt.matrix
        P = cvxopt.matrix(P, tc='d')
        q = cvxopt.matrix(q, tc='d')
        G = cvxopt.matrix(G, tc='d')
        h = cvxopt.matrix(h, tc='d')
        A = cvxopt.matrix(A, tc='d')
        b = cvxopt.matrix(b, tc='d')

        # Run optimisation
        sol = cvxopt.solvers.qp(P, q, G, h, A, b)

        # Set alpha
        alpha = np.array(sol['x']).reshape(n,2)

        # Set b
        w_dot_X = np.dot(alpha[:,0] - alpha[:,1], self.kernel(X,X))

        cond_lower = (alpha[:,1] > tol_0) | (alpha[:,0] < (self.C * (1. - tol_C)))
        cond_upper = (alpha[:,0] > tol_0) | (alpha[:,1] < (self.C * (1. - tol_C)))
        
        b_lower = np.max((-1.*self.epsilon + y - w_dot_X)[cond_lower])
        b_upper = np.min((self.epsilon + y - w_dot_X)[cond_upper])

        b = np.mean((b_lower, b_upper))
        
        # Verify b condition collapse
        b_diff = np.round(b_upper - b_lower, 5)

        if np.abs(b_diff) > 0.15:
            # Rarely, it still goes over to about 0.2 or 0.3 or something
            # print('Houston, we have a problem...', b_diff, 'with C:', self.C, 'and eps:', self.epsilon)
            pass
        
        return SVRModel(alpha, b, self.kernel, X, sub, den)



class SVRTests(unittest.TestCase):

    def setUp(self):
        self.X = np.array([[0., 0],
                           [0, 1],
                           [1, 0],
                           [1, 1]])
        self.y = np.array([0., 1, 2, 3])

    def test_polynomial(self):
        fitter = SVR(kernel=Polynomial(M=2), lambda_=0.0001, epsilon=0.1)
        m = fitter.fit(self.X, self.y)
        pred = m.predict(self.X)
        np.testing.assert_allclose(pred, self.y, atol=0.11)

    def test_rbf(self):
        fitter = SVR(kernel=RBF(sigma=0.5), lambda_=0.0001, epsilon=0.1)
        m = fitter.fit(self.X, self.y)
        pred = m.predict(self.X)
        np.testing.assert_allclose(pred, self.y, atol=0.11)

    def test_predictor_get_info(self):
        fitter = SVR(kernel=Polynomial(M=2), lambda_=0.0001, epsilon=0.1)
        m = fitter.fit(self.X, self.y)

        alpha = m.get_alpha()
        np.testing.assert_equal(alpha.shape, (4, 2))  # two alpha for each sample
        # one value in row should be much bigger than the other (the other should be zero)
        np.testing.assert_allclose(alpha.sum(axis=1), alpha.max(axis=1), rtol=1e-5)

        b = m.get_b()
        float(b)

    @patch("cvxopt.solvers.qp")
    def test_enforce_quadratic_programming_order_of_x(self, qp):
        fitter = SVR(kernel=Linear(), lambda_=0.0001, epsilon=0)
        try:
            _ = fitter.fit(self.X, self.y)
        except:
            pass
        usedq = np.array(qp.call_args[0][1]).flatten()
        np.testing.assert_equal(usedq, [0,  0, -1, 1, -2, 2, -3, 3])

    def test_kernel(self):
        for kernel in [Linear(), Polynomial(M=3), RBF(sigma=0.2)]:
            number = kernel(self.X[0], self.X[1])
            float(number)
            a1d = kernel(self.X[0], self.X)
            self.assertTrue(len(a1d.shape) == 1)
            a1d = kernel(self.X, self.X[0])
            self.assertTrue(len(a1d.shape) == 1)
            a2d = kernel(self.X, self.X)
            self.assertTrue(len(a2d.shape) == 2)


if __name__ == '__main__':
    unittest.main()

