def reduce_nq_sat(n):
    # Confirm feasibility
    if n < 4:
        return None

    # Preamble
    varnum = n ** 2
    clauses = 2*n + (n-1)*n*(n+1)
    m = n-1

    while m > 1:
        clauses += 2 * m*(m-1)
        m -= 1

    dimacs = 'p cnf ' + str(varnum) + ' ' + str(clauses) + '\n'

    # Clauses
    # Rows
    for i in range(1, varnum+1, n):
        dimacs += ' '.join([str(i+j) for j in range(n)]) + ' 0\n'

        for j in range(i, i+n-1):
            for k in range(j+1, i+n):
                dimacs += str(-j) + ' ' + str(-k) + ' 0\n'

    # Columns
    for i in range(1, n+1):
        dimacs += ' '.join([str(i+j) for j in range(0, varnum, n)]) + ' 0\n'

        for j in range(i, varnum, n):
            for k in range(j+n, varnum+1, n):
                dimacs += str(-j) + ' ' + str(-k) + ' 0\n'

    # Upper LR diagonals
    for i in range(n-1, 0, -1):
        for j in range(i, n*(n-i+1), n+1):
            for k in range(j+n+1, n*(n-i+1)+1, n+1):
                dimacs += str(-j) + ' ' + str(-k) + ' 0\n'

    # Lower LR diagonals
    for i in range(n+1, n*(n-2)+2, n):
        for j in range(i, varnum - i//n, n+1):
            for k in range(j+n+1, varnum - i//n + 1, n+1):
                dimacs += str(-j) + ' ' + str(-k) + ' 0\n'

    # Upper RL diagonals
    for i in range(2, n+1):
        for j in range(i, i + (n-1)*(i-1), n-1):
            for k in range(j+n-1, i + (n-1)*(i-1) + 1, n-1):
                dimacs += str(-j) + ' ' + str(-k) + ' 0\n'

    # Lower RL diagonals
    for i in range(2*n, varnum, n):
        for j in range(i, varnum - n + i//n, n-1):
            for k in range(j+n-1, varnum - n + i//n + 1, n-1):
                dimacs += str(-j) + ' ' + str(-k) + ' 0\n'

    return dimacs


# Test
if __name__ == '__main__':
    import sys

    sat_reduction = reduce_nq_sat(int(sys.argv[1]))
    # print(sat_reduction)

    with open(sys.argv[2], 'w') as outfile:
        outfile.write(sat_reduction)
