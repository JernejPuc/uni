"""
Decision trees and random forests
"""

import numpy as np
import random
import unittest


def random_feature(X, rand):
    return [rand.choice(list(range(X.shape[1])))]


class Node:
    def __init__(self, majority_class):
        self.majority_class = majority_class
        self.feature_idx = None
        self.threshold = None
        self.child_l = None
        self.child_r = None


class RootNode(Node):
    def __init__(self, majority_class):
        Node.__init__(self, majority_class)

    def predict(self, X):
        return np.array([self.find_leaf(x) for x in X])
    
    def find_leaf(self, x):
        curr_node = self

        while curr_node.child_l:
            if x[curr_node.feature_idx] < curr_node.threshold:
                curr_node = curr_node.child_l
            else:
                curr_node = curr_node.child_r
        
        return curr_node.majority_class


class Tree:
    def __init__(self,
                 rand=np.random.default_rng(),
                 get_candidate_columns=lambda x,_: np.arange(x.shape[1]),
                 min_samples=2
                ):

        self.rand = rand
        self.get_candidate_columns = get_candidate_columns
        self.min_samples = max(min_samples, 2)
    
    def build(self, X, y):
        self.n_classes = len(np.unique(y))

        return self.cart(X, y, root=True)

    def cart(self, X, y, root=False):
        """
        Implementation of the CART algorithm
        (CART - classification and regression trees).

        Should initially be called with root
        and without it on subsequent calls.
        """

        n_per_class = [np.sum(y == i) for i in range(self.n_classes)]
        majority_class = np.argmax(n_per_class)

        node = RootNode(majority_class) if root else Node(majority_class)

        if n_per_class[majority_class] < sum(n_per_class) and y.size >= self.min_samples \
        and len(np.unique(X, axis=0)) > 1:
            threshold = None

            # Threshold is None if all values are the same
            # Forces resample until some feature has different values
            while threshold is None:
                index_pool = self.get_candidate_columns(X, self.rand)
                feature_idx, threshold = self.split(X, y, n_per_class, index_pool)
            
            indices_l = X[:,feature_idx] < threshold
            indices_r = ~indices_l

            Xl, yl = X[indices_l], y[indices_l]
            Xr, yr = X[indices_r], y[indices_r]

            node.feature_idx = feature_idx
            node.threshold = threshold
            node.child_l = self.cart(Xl, yl)
            node.child_r = self.cart(Xr, yr)

        return node
    
    def split(self, X, y, n_per_class, index_pool):
        """
        The optimal split is found by:
            1) sorting the thresholds (values of X per each feature),
            2) looping over the thresholds,
            3) splitting the counts per class among the left and right splits for each threshold
               (performed incrementally to avoid unnecessary overhead),
            4) computing the associated gini index,
            5) returning the associated feature index
               and setting the threshold at the middle of the chosen interval
        """

        n = y.size

        gini = 1.
        feature_idx = None
        threshold = None

        for idx in index_pool:
            thresholds = X[:,idx]
            indices_t = thresholds.argsort()

            thresholds, classes = thresholds[indices_t], y[indices_t]

            n_pc_r = np.array(n_per_class)
            n_pc_l = np.zeros(n_pc_r.size)

            for i in range(1, n):
                n_pc_l[classes[i-1]] += 1
                n_pc_r[classes[i-1]] -= 1

                if thresholds[i] == thresholds[i-1]:
                    continue

                ratios_l = n_pc_l / i
                ratios_r = n_pc_r / (n-i)

                gini_l = 1. - np.sum(np.power(ratios_l, 2))
                gini_r = 1. - np.sum(np.power(ratios_r, 2))

                gini_w = (i/n)*gini_l + (1-i/n)*gini_r

                if gini_w == 0.:
                    return idx, (thresholds[i-1] + thresholds[i])/2
                
                elif gini_w < gini:
                    gini = gini_w
                    feature_idx = idx
                    threshold = (thresholds[i-1] + thresholds[i])/2

        return feature_idx, threshold


class Forest:
    """
    Predictions based on the consensus of the majority of individual decision trees.
    """

    def __init__(self):
        self.trees = []
    
    def predict_proba(self, X):
        predictions_per_sample = np.array([tree.predict(X) for tree in self.trees]).T

        p = np.zeros((X.shape[0], len(np.unique(predictions_per_sample))))

        for i in range(p.shape[1]):
            p[:,i] = np.sum(predictions_per_sample == i, axis=1) / predictions_per_sample.shape[1]
        
        return p

    # def predict(self, X):
    #     return np.argmax(self.predict_proba(X), axis=1)

    def predict(self, X):
        predictions_per_sample = np.array([tree.predict(X) for tree in self.trees]).T
        
        return np.array([self.majority_vote(*np.unique(prediction, return_counts=True)) \
                         for prediction in predictions_per_sample])
    
    @staticmethod
    def majority_vote(classes, counts):
        return classes[counts.argmax()]


class RandomForest:
    def __init__(
        self,
        rand=np.random.default_rng(),
        n=50,
        min_samples=2
        ):

        self.rand = rand
        self.n = n
        self.min_samples = min_samples
    
    @staticmethod
    def get_candidate_columns(X, rand):
        indices_f = np.arange(X.shape[1])
        rand.shuffle(indices_f)

        nf = int(np.sqrt(X.shape[1]))
        subsample_f = indices_f[:nf]

        return subsample_f

    def build(self, X, y):
        forest = Forest()
        tree_builder = Tree(rand=self.rand,
                            get_candidate_columns=self.get_candidate_columns,
                            min_samples=self.min_samples
                           )

        indices = np.arange(y.size)

        for _ in range(self.n):
            subsample = [self.rand.choice(indices) for _ in indices]
            
            Xs, ys = X[subsample,:], y[subsample]

            forest.trees.append(tree_builder.build(Xs, ys))

        return forest


class TreeTests(unittest.TestCase):

    def setUp(self):
        self.X = np.array([[0, 0],
                      [0, 1],
                      [1, 0],
                      [1, 1]])
        self.y = np.array([0, 0, 1, 1])
        self.train = self.X[:3], self.y[:3]
        self.test = self.X[3:], self.y[3:]

    def test_call_tree(self):
        t = Tree(rand=random.Random(1),
                 get_candidate_columns=random_feature,
                 min_samples=2)
        p = t.build(self.X, self.y)
        pred = p.predict(self.X)
        np.testing.assert_equal(pred, self.y)

    def test_call_randomforest(self):
        rf = RandomForest(rand=random.Random(0),
                          n=20,
                          min_samples=2)
        p = rf.build(self.X, self.y)
        pred = p.predict(self.X)
        np.testing.assert_equal(pred, self.y)


if __name__ == '__main__':
    unittest.main()
