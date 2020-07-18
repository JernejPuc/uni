"""
Artificial neural networks
"""

import numpy as np
from scipy.optimize import fmin_l_bfgs_b
import unittest


#########################################################################################################
# ANN

class ANN:
    """
    Multi-layer fully-connected artificial neural networks for classification and regression.

    Modes:
        - Multi-class classification uses categorical cross-entropy as the loss and softmax
          as the activation function for the multiple (> 2) nodes in the final layer.
        - Binary classification uses binary cross-entropy as the loss and sigmoid as the
          activation function for the sole node in the final layer.
        - Regression uses mean squared error as the loss and identity as the activation
          function for the sole node in the final layer.

    Notes:
        - All nodes in the hidden layers use sigmoid as their activation function.
        - L-BFGS is used to optimise the weights.
    """

    def __init__(self, mode, units=[], lambda_=1e-4, seed=None):
        self.mode = mode
        self.units = [None] + units + [None]
        self.lambda_ = lambda_
        self.rng = np.random.default_rng(seed=seed)

        self.activation_fn = [self.sigmoid]*len(units) + [None]
        self.w = None
        self.b = None

        self.nodes_after_input = 0
        self.w_delim = None
        self.b_delim = None

        self.mean = 0.
        self.std = 1.

        if mode == 'classification':
            self.loss = self.log_loss
        elif mode == 'regression':
            self.loss = self.mse
        else:
            raise ValueError('Unrecognised mode.')

    @staticmethod
    def sigmoid(X):
        return 1. / (1. + num_exp(-1. * X))
    
    @staticmethod
    def softmax(X):
        exp_X = num_exp(X)
        return exp_X / np.sum(exp_X, axis=1).reshape(-1, 1)
    
    @staticmethod
    def identity(X):
        return X
    
    @staticmethod
    def sigmoid_derivative(X):
        return X * (1. - X)

    @staticmethod
    def log_loss(y, y_m):
        if len(y.shape) == 1:
            y_m = y_m.flatten()
            loss_val = np.dot(y, num_log(y_m)) + np.dot(1.-y, num_log(1.-y_m))
        else:
            loss_val = np.sum(y * num_log(y_m))

        return loss_val * (-1./y.shape[0])
    
    @staticmethod
    def mse(y, y_m):
        # 1/2 factor is for convenience of the derivative
        return 0.5 * np.mean(np.power(y - y_m.flatten(), 2))
    
    def params_to_wb(self, params):
        w_flat, b_flat = params[:-self.nodes_after_input], params[-self.nodes_after_input:]

        for i in range(len(self.b)):
            self.w[i] = (w_flat[self.w_delim[i]:self.w_delim[i+1]]).reshape(self.units[i], self.units[i+1])
            self.b[i] = b_flat[self.b_delim[i]:self.b_delim[i+1]]

        return None
    
    def wb_to_params(self, w=None, b=None):
        if w is None:
            w = self.w
        if b is None:
            b = self.b
        
        return np.concatenate([arr.flatten() for arr in w + b])

    def loss_grad(self, a_i, delta, w_i, n_samples):
        w_grad = (np.dot(a_i.T, delta) + self.lambda_ * w_i) / n_samples
        b_grad = np.mean(delta, axis=0)

        return w_grad, b_grad

    def forprop(self, X, return_all=True):
        a = [X]

        for i in range(len(self.b)):
            a.append(self.activation_fn[i](np.dot(a[-1], self.w[i]) + self.b[i]))
        
        # Prediction calls only need the last layer returned
        if return_all:
            return a
        else:
            return a[-1]

    def backprop(self, X, y):
        n_samples = X.shape[0]

        # Loss
        activations = self.forprop(X)
        loss = self.loss(y, activations[-1])

        # Regularisation
        loss += 0.5 * self.lambda_ * np.sum([np.dot(w_i.flatten(), w_i.flatten()) for w_i in self.w]) / n_samples

        # Gradients
        w_grad = [None]*(len(self.units)-1)
        b_grad = [None]*(len(self.units)-1)
        deltas = [None]*(len(self.units)-1)

        deltas[-1] = activations[-1] - y.reshape(activations[-1].shape)
        w_grad[-1], b_grad[-1] = self.loss_grad(activations[-2], deltas[-1], self.w[-1], n_samples)

        for i in range(len(self.units)-2, 0, -1):
            deltas[i-1] = np.dot(deltas[i], self.w[i].T) * self.sigmoid_derivative(activations[i])
            w_grad[i-1], b_grad[i-1] = self.loss_grad(activations[i-1], deltas[i-1], self.w[i-1], n_samples)

        return loss, w_grad, b_grad

    def loss_grad_lbfgs(self, params, X, y):
        self.params_to_wb(params)

        loss, w_grad, b_grad = self.backprop(X, y)

        grad = self.wb_to_params(w_grad, b_grad)

        return loss, grad

    def init_weights(self):
        """
        Init. such that activations are in the linear part of the sigmoid, i.e. z is within [-1,1]
        for z = sum(w) + b.

        In the worst case, all prior activations are either -1 or 1. Therefore, the sum of all
        absolute contributions, sum(abs(w)) + abs(b), should be less than 1.

        After initialising the weights within [-0.5, 0.5] range, this is achieved by setting:
            - w /= sum(abs(w)) + abs(b)
            - b /= sum(abs(w)) + abs(b)
        """

        w = [-0.5 + self.rng.random((self.units[i], self.units[i+1])) for i in range(len(self.units)-1)]
        b = [-0.5 + self.rng.random(self.units[i+1]) for i in range(len(self.units)-1)]

        for i in range(len(b)):
            for j in range(b[i].size):
                wb_ij = np.concatenate((w[i][:,j], (b[i][j],)))
                wb_ij_den = np.sum(np.abs(wb_ij))

                w[i][:,j] /= wb_ij_den
                b[i][j] /= wb_ij_den
        
        self.w = w
        self.b = b

        # Delimiters for convenience when unpacking params
        self.w_delim = np.cumsum([0] + [w_i.size for w_i in w])
        self.b_delim = np.cumsum([0] + [b_i.size for b_i in b])

        return None

    def weights(self):
        return [np.vstack(wb_i) for wb_i in zip(self.w, self.b)]

    def standardise(self, X):
        """
        Sets internal mean and std, so that further preprocessing at inference time is not needed.
        Returns preprocessed so that other models can use it as well.
        """

        self.mean = np.mean(X, axis=0)
        self.std = np.std(X, axis=0)

        return (X - self.mean) / self.std
    
    def fit(self, X, y):
        # Input size
        self.units[0] = X.shape[1]

        # Output size and activation (and y to proba)
        if self.mode == 'classification':
            n_labels = np.unique(y).size

            if n_labels > 2:
                self.units[-1] = n_labels
                self.activation_fn[-1] = self.softmax
                
                y_proba = np.zeros((y.size, n_labels))

                for idx in range(n_labels):
                    y_proba[y == idx, idx] = 1.
                
                y = y_proba
            
            else:
                self.units[-1] = 1
                self.activation_fn[-1] = self.sigmoid

                y = y.astype(np.float64)

        else:
            self.units[-1] = 1
            self.activation_fn[-1] = self.identity
        
        # For convenience when unpacking params
        self.nodes_after_input = sum(self.units[1:])

        # Random initialisation
        self.init_weights()

        # Optimisation
        opt_params, _, _ = fmin_l_bfgs_b(
            self.loss_grad_lbfgs,
            self.wb_to_params(),
            args=(X,y)
        )

        self.params_to_wb(opt_params)
        
        return self


class ANNClassification(ANN):
    """
    Sets the mode to classification and extends the parent ANN class
    with the corresponding prediction method.
    """

    def __init__(self, units=[], lambda_=1e-4, seed=None):
        ANN.__init__(self, mode='classification', units=units, lambda_=lambda_, seed=seed)

    # def predict_proba(self, X):
    def predict(self, X):
        p = self.forprop((X-self.mean)/self.std, return_all=False)

        if self.units[-1] == 1:
            p = np.hstack((1.-p, p))
        
        return p
    
    # def predict(self, X):
    def predict_label(self, X):
        # p = self.predict_proba(X)
        p = self.predict(X)
        return np.argmax(p, axis=1)


class ANNRegression(ANN):
    """
    Sets the mode to regression and extends the parent ANN class
    with the corresponding prediction method.
    """

    def __init__(self, units=[], lambda_=1e-4, seed=None):
        ANN.__init__(self, mode='regression', units=units, lambda_=lambda_, seed=seed)

    def predict(self, X):
        return self.forprop((X-self.mean)/self.std, return_all=False).flatten()


#########################################################################################################
# Util.

def num_log(X, eps=1e-16):
    """
    Enforces valid range.
    """

    return np.log(np.clip(X, eps, 1.))


def num_exp(X, lim=37):
    """"
    Prevents overflow.
    """

    return np.exp(np.clip(X, -lim, lim))


#########################################################################################################
# Test

class ANNTests(unittest.TestCase):

    def setUp(self):
        self.X = np.array([[0, 0],
                           [0, 1],
                           [1, 0],
                           [1, 1]])
        self.y = np.array([0, 1, 2, 3])
        self.hard_y = np.array([0, 1, 1, 0])

    def test_ann_classification_no_hidden_layer(self):
        fitter = ANNClassification(units=[], lambda_=0.0001)
        m = fitter.fit(self.X, self.y)
        pred = m.predict(self.X)
        self.assertEqual(pred.shape, (4, 4))
        np.testing.assert_allclose(pred, np.identity(4), atol=0.01)

    def test_ann_classification_no_hidden_layer_hard(self):
        fitter = ANNClassification(units=[], lambda_=0.0001)
        m = fitter.fit(self.X, self.hard_y)
        pred = m.predict(self.X)
        self.assertEqual(pred.shape, (4, 2))
        np.testing.assert_allclose(pred, 0.5, atol=0.01)

    def test_ann_classification_hidden_layer_hard(self):
        fitter = ANNClassification(units=[10], lambda_=0.0001)
        m = fitter.fit(self.X, self.hard_y)
        pred = m.predict(self.X)
        self.assertEqual(pred.shape, (4, 2))
        np.testing.assert_allclose(pred, [[1, 0], [0, 1], [0, 1], [1, 0]], atol=0.01)

    def test_ann_classification_hidden_layers_hard(self):
        fitter = ANNClassification(units=[10, 20], lambda_=0.0001)
        m = fitter.fit(self.X, self.hard_y)
        pred = m.predict(self.X)
        self.assertEqual(pred.shape, (4, 2))
        np.testing.assert_allclose(pred, [[1, 0], [0, 1], [0, 1], [1, 0]], atol=0.01)

    def test_ann_regression_no_hidden_layer(self):
        fitter = ANNRegression(units=[], lambda_=0.0001)
        m = fitter.fit(self.X, self.y)
        pred = m.predict(self.X)
        self.assertEqual(pred.shape, (4,))
        np.testing.assert_allclose(pred, self.y, atol=0.01)

    def test_ann_regression_no_hidden_layer_hard(self):
        fitter = ANNRegression(units=[], lambda_=0.0001)
        m = fitter.fit(self.X, self.hard_y)
        pred = m.predict(self.X)
        self.assertEqual(pred.shape, (4,))
        np.testing.assert_allclose(pred, 0.5, atol=0.01)

    def test_ann_regression_hidden_layer_hard(self):
        fitter = ANNRegression(units=[10], lambda_=0.0001)
        m = fitter.fit(self.X, self.hard_y)
        pred = m.predict(self.X)
        self.assertEqual(pred.shape, (4,))
        np.testing.assert_allclose(pred, self.hard_y, atol=0.01)

    def test_ann_regression_hidden_layers_hard(self):
        fitter = ANNRegression(units=[10,10], lambda_=0.0001)
        m = fitter.fit(self.X, self.hard_y)
        pred = m.predict(self.X)
        self.assertEqual(pred.shape, (4,))
        np.testing.assert_allclose(pred, self.hard_y, atol=0.01)

    def test_predictor_get_info(self):
        fitter = ANNRegression(units=[10, 5], lambda_=0.0001)
        m = fitter.fit(self.X, self.y)
        lw = m.weights()

        self.assertEqual(len(lw), 3)

        self.assertEqual(lw[0].shape, (3, 10))
        self.assertEqual(lw[1].shape, (11, 5))
        self.assertEqual(lw[2].shape, (6, 1))


if __name__ == '__main__':
    unittest.main()

