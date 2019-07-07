"""Karger analysis v19.01.19 by JernejP"""

import sys
import numpy as np
from time import perf_counter
import matplotlib.pyplot as plt


###############################################################################
# Core

def karger(numv, edges):
    """Efficient Karger's alg. implementation (numpy magic)"""
    E = edges.copy()

    for _ in range(numv-2):
        u, v = E[np.random.randint(E.shape[0]),:]
        E = E[np.logical_and(np.logical_or(E[:,0] != u, E[:,1] != v),
                             np.logical_or(E[:,1] != u, E[:,0] != v))]
        E = np.where(E==v, u, E)

    return E.shape[0]


def cut_groups(numv, edges):
    """Yield vertex groups of a random cut (mincut or close)"""
    V = {}
    E = edges.copy()

    for i in range(1, numv+1):
        V[i] = [i]
    
    for _ in range(numv-2):
        u, v = E[np.random.randint(E.shape[0]),:]
        V[u] += V[v]
        del V[v]
        E = E[np.logical_and(np.logical_or(E[:,0] != u, E[:,1] != v),
                             np.logical_or(E[:,1] != u, E[:,0] != v))]
        E = np.where(E==v, u, E)
    
    groups = list(V.values())
    return sorted(groups[0]), sorted(groups[1])


def mincut_groups(numv, edges, opt):
    """Index until mincut found, then unpack into vertex groups"""
    H = np.empty([1,numv-2], dtype=int)
    cut = 0

    while cut != opt:
        E = edges.copy()

        for t in range(numv-2):
            i = np.random.randint(E.shape[0])
            u, v = E[i,:]
            E = E[np.logical_and(np.logical_or(E[:,0] != u, E[:,1] != v),
                                 np.logical_or(E[:,1] != u, E[:,0] != v))]
            E = np.where(E==v, u, E)
            H[0,t] = i
        
        cut = E.shape[0]
    
    V = {}
    E = edges.copy()

    for i in range(1, numv+1):
        V[i] = [i]
    
    for t in range(numv-2):
        u, v = E[H[0,t],:]
        V[u] += V[v]
        del V[v]
        E = E[np.logical_and(np.logical_or(E[:,0] != u, E[:,1] != v),
                             np.logical_or(E[:,1] != u, E[:,0] != v))]
        E = np.where(E==v, u, E)
    
    groups = list(V.values())
    return sorted(groups[0]), sorted(groups[1])


###############################################################################
# Analytics

def single_run(numv, edges, opt):
    """Produce typical run data for later analysis"""
    cuts = []
    bestcuts = []
    bestcut = len(edges)
    j = 0

    while bestcut != opt:
        cut = karger(numv, edges)
        cuts.append(cut)

        if cut < bestcut:
            bestcut = cut
        
        bestcuts.append(bestcut)

        j += 1
        print(j, cut)
    
    return cuts, bestcuts


def multi_run(numv, edges, opt, runs):
    """Produce distribution data for later analysis"""
    run_data = np.empty([1,runs], dtype=int)

    for i in range(runs):
        cut = 0
        j = 0

        while cut != opt:
            cut = karger(numv, edges)
            j += 1
        
        run_data[0,i] = j
        print(i, j)

    return np.bincount(run_data[0,:])


def time_run(numv, edges):
    """Measure single pass performance"""
    t = perf_counter()
    cut = karger(numv, edges)

    return perf_counter() - t, cut


###############################################################################
# Plotting

def plot_run(*args, filename=None):
    """Plot cut alternation and convergence of solution"""
    if filename is None:
        cuts = args[0]
        bestcuts = args[1]
    else:
        with open(filename, 'r') as infile:
            line = infile.readline()[:-1].split(' ')
            cuts = list(map(lambda x: int(x), line))
            
            line = infile.readline()[:-1].split(' ')
            bestcuts = list(map(lambda x: int(x), line))

    runs = [i for i in range(len(cuts))]

    fig = plt.figure(figsize=(6.4, 4.8), dpi=150)
    ax = fig.gca()
    plt.plot(runs, cuts, runs, bestcuts)
    ax.set_xlabel('Runs')
    ax.set_ylabel('Cuts')
    ax.legend(['Current cut', 'Best cut'])
    ax.set_title('Typical run')
    fig.tight_layout()
    plt.show()


def histogram(*args, filename=None):
    """Plot histogram from distribution data"""
    if filename is None:
        distribution = args[0]
    else:
        with open(filename, 'r') as infile:
            line = infile.readline()[:-1].split(' ')
            distribution = list(map(lambda x: int(x), line))

    run_data = []
    total_runs = len(distribution)

    for i in range(total_runs):
        runs = distribution[i]
        run_data += [i]*runs

    sigma = np.std(list(map(lambda x: -x, run_data[1::-1])) + run_data)
    gauss = ((1/(np.sqrt(2*np.pi) * sigma)) *
            np.exp(-0.5 * (1/sigma * np.array(run_data))**2)) * 2*len(run_data)

    fig = plt.figure(figsize=(6.4, 4.8), dpi=150)
    ax = fig.gca()
    ax.hist(run_data, total_runs)
    ax.plot(run_data, gauss)
    ax.set_xlabel('Runs before optimum')
    ax.set_ylabel('Run density')
    ax.set_title('Run distribution')
    fig.tight_layout()
    plt.show()


###############################################################################
# Auxiliary

def read_graph(filename, dtype='uint16'):
    """Read graph data in DIMACS format"""
    with open(filename, 'r') as infile:
        numv, nume = infile.readline()[:-1].split(' ')[2:]
        numv, nume = (int(numv), int(nume))

        edges = np.empty([nume, 2], dtype=dtype)

        for edge in range(nume):
            v1, v2 = infile.readline()[2:-1].split(' ')
            v1, v2 = (int(v1), int(v2))
            edges[edge] = np.array([v1, v2])
    
    return numv, edges


def tlog(cuts, bestcuts, filename='ka_tlog.txt'):
    """Log typical run data"""
    lines = ' '.join(list(map(lambda x: str(x), cuts))) + '\n'
    lines += ' '.join(list(map(lambda x: str(x), bestcuts))) + '\n'

    with open(filename, 'w') as outfile:
        outfile.write(lines)


def mlog(distribution, filename='ka_mlog.txt'):
    """Log multi run data"""
    line = ' '.join(list(map(lambda x: str(x), distribution))) + '\n'

    with open(filename, 'w') as outfile:
        outfile.write(line)


def merge_mlogs(filenames, filename='ka_mlog.txt'):
    """Combine separate distribution data"""
    run_data = []

    for datafile in filenames:
        with open(datafile, 'r') as infile:
            line = infile.readline()[:-1].split(' ')
            run_data.append(list(map(lambda x: int(x), line)))

    maxlen = max([len(distribution) for distribution in run_data])
    distributions = np.array([rd + [0]*(maxlen-len(rd)) for rd in run_data])
    
    new_distribution = np.sum(distributions, axis=0)
    new_data = ' '.join(list(map(lambda x: str(x), new_distribution))) + '\n'
    
    with open(filename, 'w') as outfile:
        outfile.write(new_data)


#############################################################################
# CLI

if __name__ == '__main__':
    mode = sys.argv[1]

    if mode == 'single':
        datafile = sys.argv[2]
        numv, edges = read_graph(datafile)
        opt = int(sys.argv[3])
        cuts, bestcuts = single_run(numv, edges, opt)
        
        if len(sys.argv) > 4:
            tlog(cuts, bestcuts, filename=sys.argv[4])
        else:
            plot_run(cuts, bestcuts)
    
    elif mode == 'plot':
        datafile = sys.argv[2]
        plot_run(filename=datafile)

    elif mode == 'multi':
        datafile = sys.argv[2]
        numv, edges = read_graph(datafile)
        opt = int(sys.argv[3])
        runs = int(sys.argv[4])
        distribution = multi_run(numv, edges, opt, runs)

        if len(sys.argv) > 5:
            mlog(distribution, filename=sys.argv[5])
        else:
            histogram(distribution)
    
    elif mode == 'merge':
        datafile = sys.argv[2]
        datafiles = sys.argv[3:]
        merge_mlogs(datafiles, filename=datafile)

    elif sys.argv[1] == 'hist':
        datafile = sys.argv[2]
        histogram(filename=datafile)
    
    elif mode == 'time':
        datafile = sys.argv[2]
        numv, edges = read_graph(datafile)
        dt, _ = time_run(numv, edges)
        print(dt)
    
    elif mode == 'cut':
        datafile = sys.argv[2]
        numv, edges = read_graph(datafile)
        a, b = cut_groups(numv, edges)
        print(a)
        print(b)
    
    elif mode == 'mincut':
        datafile = sys.argv[2]
        opt = int(sys.argv[3])
        numv, edges = read_graph(datafile)
        a, b = mincut_groups(numv, edges, opt)
        print(a)
        print(b)
    
    else:
        print("Mode '" + mode + "' not implemented.")
