"""Hidden Markov Models : Supervised, Viterbi & Baum-Welch learning."""
__author__ = 'JP'


########################################################################################################################
# Mods

import numpy as np


########################################################################################################################
# Core

def supervised(s, o, y, x):
    """Sets the transition and emission matrices by counting occurences in the given paths."""

    t = np.zeros([len(s), len(s)], dtype=int)
    e = np.zeros([len(s), len(o)], dtype=int)

    for i in range(len(y)):
        if i != len(y)-1:
            t[s.index(x[i]), s.index(x[i+1])] += 1

        e[s.index(x[i]), o.index(y[i])] += 1

    for i in range(len(s)):
        if np.sum(t[i, :]) == 0:
            t[i, :] = np.ones([1, len(s)], dtype=int)

        if np.sum(e[i, :]) == 0:
            e[i, :] = np.ones([1, len(o)], dtype=int)

    t = np.divide(t, np.sum(t, 1).reshape([len(s), 1]), dtype=float)
    e = np.divide(e, np.sum(e, 1).reshape([len(s), 1]), dtype=float)

    return t, e


def viterbi(s, o, y, pi, ti, ei, n):
    """Finds the unknown parameters of a HMM by using a supervised algorithm."""

    # Local defs
    p = np.array(pi)
    t = np.array(ti)
    e = np.array(ei)

    probability_table = np.zeros([len(s), len(y)])
    history = np.zeros([len(s), len(y)], dtype=int)

    for _ in range(n):

        # Find the "Viterbi path" (the most likely sequence of hidden states)
        for i in range(len(s)):
            probability_table[i, 0] = p[i] * e[i, o.index(y[0])]

        for i in range(1, len(y)):
            for j in range(len(s)):
                tmp = [probability_table[k, i-1] * t[k, j] * e[j, o.index(y[i])] for k in range(len(s))]

                probability_table[j, i] = max(tmp)
                history[j, i] = tmp.index(probability_table[j, i])

        final_probabilities = [probability_table[i, -1] for i in range(len(s))]
        pr_x = max(final_probabilities)
        s_index = final_probabilities.index(pr_x)

        # Hidden path
        sequence = s[s_index]

        for i in range(len(y)-1, 0, -1):
            s_index = history[s_index, i]
            sequence = s[s_index] + sequence
    
        t, e = supervised(s, o, y, sequence)

    return t, e


def forward(s, o, y, p, t, e):
    """Computes the forward probabilities."""

    # Initial distribution
    probability_table = np.zeros([len(s), len(y)])

    for i in range(len(s)):
        probability_table[i, 0] = p[i] * e[i, o.index(y[0])]

    # Formula
    for i in range(1, len(y)):
        for j in range(len(s)):
            probability_table[j,i] = sum(probability_table[k,i-1] * t[k,j] * e[j, o.index(y[i])] for k in range(len(s)))

    return probability_table, np.sum(probability_table[:, -1], 0)


def backward(s, o, y, p, t, e):
    """Computes the backward probabilities."""

    # Initial (final) distribution
    probability_table = np.zeros([len(s), len(y)])

    for i in range(len(s)):
        probability_table[i, -1] = 1.

    # Formula
    for i in range(len(y)-2, -1, -1):
        for j in range(len(s)):
            probability_table[j,i] = sum(probability_table[k,i+1] * t[j,k] * e[k, o.index(y[i+1])] for k in range(len(s)))

    return probability_table, np.sum(p * probability_table[:, 0] * e[:, o.index(y[0])])


def baum_welch(s, o, y, pi, ti, ei, n):
    """Finds the unknown parameters of a HMM by using a forward-backward algorithm."""

    # Local defs
    p = np.array(pi)
    t = np.array(ti)
    e = np.array(ei)

    # Observations indexed according to the observation space
    yoi = np.array([o.index(yi) for yi in y])

    for _ in range(n):
        # Forward-backward calls (fw = alpha, bw = beta in standard notation)
        fw, pr_fw = forward(s, o, y, p, t, e)
        bw, pr_bw = backward(s, o, y, p, t, e)

        # State probabilities (gamma in standard notation, pr_fw == pr_bw)
        state_probabilities = np.multiply(fw, bw) / pr_fw

        # Update initial probabilities
        p = state_probabilities[:,0]

        # Transient probabilities (xi in standard notation)
        transient_probabilities = np.zeros([len(y)-1, len(s), len(s)])

        for k in range(len(y)-1):
            for i in range(len(s)):
                for j in range(len(s)):
                    transient_probabilities[k, i, j] = t[i, j] * fw[i, k] * e[j, o.index(y[k+1])] * bw[j, k+1]

            transient_probabilities[k] /= np.sum(transient_probabilities[k])

        # Update transitions
        for i in range(len(s)):
            # Detail from https://en.wikipedia.org/wiki/Baum%E2%80%93Welch_algorithm#Update
            # gamma should only be summed to T-1
            divisor = np.sum(state_probabilities[i, :-1])
            # divisor = np.sum(transient_probabilities[:, i, :])

            for j in range(len(s)):
                t[i, j] = np.sum(transient_probabilities[:, i, j]) / divisor

        # Update emissions
        for i in range(len(s)):
            divisor = np.sum(state_probabilities[i, :])

            for j in range(len(o)):
                e[i, j] = np.sum(state_probabilities[i, yoi == j]) / divisor

    return t, e


########################################################################################################################
# Aux

def pretty(states, observations, transitions, emissions, alg=''):
    """Prints the matrices in a readable format."""

    out = ('\n' + alg + ':\n\n' if alg else '') + '  ' + ' '.join(states) + '\n'

    for i in range(len(states)):
        out += states[i] + ' ' + ' '.join('%.3f' % prob for prob in transitions[i, :]) + '\n'

    out += '\n  ' + ' '.join(observations) + '\n'

    for i in range(len(states)):
        out += states[i] + ' ' + ' '.join('%.3f' % prob for prob in emissions[i, :]) + '\n'

    print(out)


########################################################################################################################
# Test

if __name__ == '__main__':

    # Def
    state_space = ['A', 'B', 'C', 'D']
    observation_space = ['x', 'y', 'z']
    
    observations = list('yzzxzzzyxyxzyxzzyyzzxxzyzyyyyyyxzyxxyzzzyzxyxxxxyxzzzzyzxyxxzyyyyxyzyxzzyzyxyxzzxyxxxzxyxyyxzxxyzxzz')

    initial_transitions = [[0.056, 0.443, 0.334, 0.167],
                           [0.826, 0.052, 0.011, 0.111],
                           [0.022, 0.163, 0.811, 0.004],
                           [0.141, 0.696, 0.055, 0.108]]
    
    initial_emissions = [[0.340, 0.084, 0.576],
                         [0.196, 0.592, 0.212],
                         [0.344, 0.609, 0.047],
                         [0.320, 0.487, 0.193]]
    
    number_of_iterations = 100
    initial_probabilities = [1. / len(state_space)] * len(state_space)


    # Viterbi
    transitions, emissions = viterbi(state_space,
                                     observation_space,
                                     observations,
                                     initial_probabilities,
                                     initial_transitions,
                                     initial_emissions,
                                     number_of_iterations)

    pretty(state_space, observation_space, transitions, emissions, 'Viterbi')


    # Baum-Welch
    transitions, emissions = baum_welch(state_space,
                                        observation_space,
                                        observations,
                                        initial_probabilities,
                                        initial_transitions,
                                        initial_emissions,
                                        number_of_iterations)

    pretty(state_space, observation_space, transitions, emissions, 'Baum-Welch')
    
