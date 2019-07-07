def reduce_kg_sat(k, v, e, edges):
    # Confirm feasibility
    if k < 2:
        return None

    # Preamble
    varnum = k*v
    clauses = v*(1 + k*(k-1)//2) + k*e
    dimacs = 'p cnf ' + str(varnum) + ' ' + str(clauses) + '\n'

    # Clauses
    # Vertices
    for i in range(1, varnum+1, k):
        dimacs += ' '.join([str(i+j) for j in range(k)]) + ' 0\n'

        for m in range(i, i+k-1):
            for n in range(m+1, i+k):
                dimacs += str(-m) + ' ' + str(-n) + ' 0\n'

    # Edges
    kv = [[0, 0]]

    for i in range(1, k*v+1, k):
        kv.append(list(range(i, i+k)))

    for i in range(e):
        for j in range(k):
            dimacs += str(-kv[edges[i][0]][j]) + ' ' + str(-kv[edges[i][1]][j]) + ' 0\n'

    return dimacs


if __name__ == '__main__':
    import sys

    with open(sys.argv[1], 'r') as in_file:
        (v, e) = in_file.readline()[:-1].split(' ')[2:]
        (v, e) = (int(v), int(e))

        edges = [(0, 0)] * e

        for edge in range(e):
            (v1, v2) = in_file.readline()[2:-1].split(' ')
            edges[edge] = (int(v1), int(v2))

    sat_reduction = reduce_kg_sat(int(sys.argv[2]), v, e, edges)
    # print(sat_reduction)

    with open(sys.argv[3], 'w') as outfile:
        outfile.write(sat_reduction)
