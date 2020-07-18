"""
Kernelised ridge regression
"""

import numpy as np
import unittest


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

    def __init__(self, M=3):
        self.M = M

    def __call__(self, A, B):
        return np.power(np.dot(A, B.T) + 1., self.M)


class RBF:
    """
    Radial basis function kernel.
    """

    def __init__(self, sigma=0.2):
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


class KRRModel:
    """
    Kernelised ridge regression -- fitted model.
    """

    def __init__(self, alpha, kernel, X):
        self.alpha = alpha
        self.kernel = kernel
        self.X = X

    def predict(self, X):
        return np.dot(self.alpha, self.kernel(X, self.X).T)


class KernelisedRidgeRegression:
    """
    Kernelised ridge regression -- fitting.
    """

    def __init__(self, kernel=RBF(sigma=0.2), lambda_=1e-4):
        self.kernel = kernel
        self.lambda_ = lambda_
    
    def fit(self, X, y):
        alpha = np.dot(np.linalg.inv(self.kernel(X, X) + self.lambda_ * np.eye(X.shape[0])), y)
        
        return KRRModel(alpha, self.kernel, X)


class KRRTests(unittest.TestCase):

    def setUp(self):
        self.X = np.array([[0, 0],
                           [0, 1],
                           [1, 0],
                           [1, 1]])
        self.y = np.array([0, 0, 1, 1])
    
    def test_polynomial(self):
        fitter = KernelisedRidgeRegression(kernel=Polynomial(M=2), lambda_=0.0001)
        m = fitter.fit(self.X, self.y)
        pred = m.predict(self.X)
        np.testing.assert_almost_equal(pred, self.y, decimal=3)

    def test_rbf(self):
        fitter = KernelisedRidgeRegression(kernel=RBF(sigma=0.5), lambda_=0.0001)
        m = fitter.fit(self.X, self.y)
        pred = m.predict(self.X)
        np.testing.assert_almost_equal(pred, self.y, decimal=3)

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

