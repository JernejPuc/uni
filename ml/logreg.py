"""
Logistic regression
"""

import numpy as np
from scipy.optimize import fmin_l_bfgs_b as fmin
import unittest


def sigmoid(X):
    return np.divide(np.ones_like(X), 1 + np.exp(-X))


def logreg_log_likelihood(beta,X,y,eps=1e-16):
    p = sigmoid(np.dot(beta,X.T))
    logp = np.log(p + eps).T
    log1_p = np.log(1-p + eps).T

    return np.dot(y, logp) + np.dot(1-y, log1_p)


def logreg_ll_grad(beta,X,y):
    """
    Log likelihood gradient.
    """

    return np.dot(y - sigmoid(np.dot(beta,X.T)), X)


def logreg_mle(X,y,approx_grad=False):
    """
    Maximum likelihood estimation.
    """

    n,m = X.shape
    
    if approx_grad:
        fprime = None
    else:
        fprime = lambda beta: (-1./n) * logreg_ll_grad(beta,X,y)
    
    beta_opt, _, _ = fmin(func=lambda beta: (-1./n) * logreg_log_likelihood(beta,X,y),
                          x0=(1./m) * np.ones(m),
                          fprime=fprime,
                          approx_grad=approx_grad
    )

    return beta_opt


def logreg_predict(beta,X):
    return (sigmoid(np.dot(beta,X.T)) > 0.5).astype(int)


class LogRegTests(unittest.TestCase):

    def setUp(self):
        self.X = np.array([[0, 0],
                      [0, 1],
                      [1, 0],
                      [1, 1]])
        self.y = np.array([0, 1, 1, 1])

    def test_logreg_log_likelihood(self):
        ll = logreg_log_likelihood(np.array([0.5, 0.1]), self.X, self.y)
        self.assertIsInstance(ll, float)
        self.assertLess(ll, 0)

    def test_logreg_mle(self):
        mle = logreg_mle(self.X, self.y)
        self.assertIsInstance(mle, np.ndarray)
        self.assertEqual(mle.shape, (2,))

    def test_logreg_predict(self):
        pred = logreg_predict(np.array([10, 10]), self.X)
        np.testing.assert_equal(pred, self.y)


if __name__ == '__main__':
    unittest.main()

