"""Ham sandwich cut application for a special case of the "ghosts vs busters" problem by Jernej Puc"""


######################################################################################################
# Mods

from time import perf_counter
import matplotlib.pyplot as plt
from random import uniform, triangular, choice, random


######################################################################################################
# Main

def hscut(P,Q):
    """Does the thing in O(n*log(n)**2). Probably."""

    # Termination
    if len(P) < 2:
        return P,Q

    # Edges of the initial interval
    # maxpx = max(P, key=lambda p: p[1])
    # minqx = min(Q, key=lambda q: q[1])
    # maxpy = max(P, key=lambda p: p[2])
    # minqy = min(Q, key=lambda q: q[2])
    # minpy = min(P, key=lambda p: p[2])
    # maxqy = max(Q, key=lambda q: q[2])
    # a = (maxpy[2] - minqy[2])/(maxpx[1] - minqx[1]) - 1.
    # b = (maxqy[2] - minpy[2])/(minqx[1] - maxpx[1]) + 1.

    # If possible, one can begin with sufficiently large constants instead
    a = -10000.
    b = 10000.

    # Determines the interval, at which we first filter the lists of candidate lines
    dba = 10.

    # High median levels on the edges of the initial interval
    m = len(P)//2
    mpa = sorted(p[1]*a - p[2] for p in P)[m]
    mqa = sorted(q[1]*a - q[2] for q in Q)[m]
    mpb = sorted(p[1]*b - p[2] for p in P)[m]
    mqb = sorted(q[1]*b - q[2] for q in Q)[m]

    # Initial sets of candidate lines
    Pi = P.copy()
    Qi = Q.copy()

    # Find median level intersection
    while len(Pi) > 1 or len(Qi) > 1:
        
        # Reduce the interval
        while True:
            pr = choice(Pi)
            qr = choice(Qi)
            c = (pr[2]-qr[2])/(pr[1]-qr[1])
            
            if a < c < b:
                c += (random() - 0.5)*1e-8
                break

        mpc = sorted(p[1]*c - p[2] for p in P)[m]
        mqc = sorted(q[1]*c - q[2] for q in Q)[m]

        if (mpa-mqa)*(mpc-mqc) < 0.:
            b = c
            mpb = mpc
            mqb = mqc
        else:
            a = c
            mpa = mpc
            mqa = mqc
        
        # Reduce the sets of candidate lines
        if len(Pi) < 4 or b-a < dba:
            dba = (b-a)/4.
            
            k = (mqb - mpa)/(b - a)
            n = mpa - k*a

            Pi = [p for p in Pi if mpa >= p[1]*a - p[2] >= mqa or mqb >= p[1]*b - p[2] >= mpb or a < (n+p[2])/(p[1]-k) < b]
            Qi = [q for q in Qi if mpa >= q[1]*a - q[2] >= mqa or mqb >= q[1]*b - q[2] >= mpb or a < (n+q[2])/(q[1]-k) < b]
    
    # Bisector
    k = (Qi[0][2] - Pi[0][2])/(Qi[0][1] - Pi[0][1])
    n = Pi[0][2] - k*Pi[0][1] + 1e-12
    
    # Divide point clouds
    Pu = [p for p in P if p[2] > k*p[1] + n]
    Qu = [q for q in Q if q[2] > k*q[1] + n]
    Pl = [p for p in P if p[2] <= k*p[1] + n]
    Ql = [q for q in Q if q[2] <= k*q[1] + n]

    Pu, Qu = hscut(Pu, Qu)
    Pl, Ql = hscut(Pl, Ql)

    # Merge solutions
    return Pu+Pl, Qu+Ql


######################################################################################################
# Backup

def chull(R):
    """Does the thing in O(n**2), but works well for a general case (collinearity included)."""

    # Data structures
    P = []
    Q = []
    upper = []
    lower = []
    R.sort(key=lambda r: r[1])

    # Overhead    
    upper_clear = upper.clear
    lower_clear = lower.clear
    upper_pop = upper.pop
    lower_pop = lower.pop
    upper_append = upper.append
    lower_append = lower.append
    R_remove = R.remove
    P_append = P.append
    Q_append = Q.append

    # (Re)construct the convex hull by using an optimised Andrew's monotone chain algorithm
    while len(R):
        # Endpoints
        r0 = R[0]
        rn = R[-1]

        # Initial heuristic line
        ku = (rn[2] - r0[2])/(rn[1] - r0[1] + 1e-12)
        nu = r0[2] - ku*r0[1]
        kl = ku
        nl = nu

        # Add left endpoint
        upper_clear()
        lower_clear()
        upper_append(r0)
        lower_append(r0)
        lenu = 1
        lenl = 1

        # Construct the middle of the upper and lower hull sections
        for j in range(1,len(R)-1):
            r = R[j]

            if r[2] > ku*r[1] + nu:
                while lenu > 1 and (upper[-1][1] - upper[-2][1])*(r[2] - upper[-2][2]) > (r[1] - upper[-2][1])*(upper[-1][2] - upper[-2][2]):
                    upper_pop()
                    lenu-=1
                upper_append(r)
                lenu+=1

                # Update upper heuristic line
                ku = (rn[2] - r[2])/(rn[1] - r[1] + 1e-12)
                nu = r[2] - ku*r[1]

            elif r[2] < kl*r[1] + nl:
                while lenl > 1 and (lower[-1][1] - lower[-2][1])*(r[2] - lower[-2][2]) < (r[1] - lower[-2][1])*(lower[-1][2] - lower[-2][2]):
                    lower_pop()
                    lenl-=1
                lower_append(r)
                lenl+=1

                # Update lower heuristic line
                kl = (rn[2] - r[2])/(rn[1] - r[1] + 1e-12)
                nl = r[2] - kl*r[1]
        
        # Add right endpoint
        while lenu > 1 and (upper[-1][1] - upper[-2][1])*(rn[2] - upper[-2][2]) > (rn[1] - upper[-2][1])*(upper[-1][2] - upper[-2][2]):
            upper_pop()
            lenu-=1
        
        while lenl > 1 and (lower[-1][1] - lower[-2][1])*(rn[2] - lower[-2][2]) < (rn[1] - lower[-2][1])*(lower[-1][2] - lower[-2][2]):
            lower_pop()
            lenl-=1
        
        upper_append(rn)
        lower_append(rn)
        lenu+=1
        lenl+=1
        
        # Find the vertices of both ordinate crossing edges        
        u1,u2 = [(upper[j],upper[j+1]) for j in range(lenu) if upper[j][1] < 0. and upper[j+1][1] > 0.][0]
        v1,v2 = [(lower[j],lower[j+1]) for j in range(lenl) if lower[j][1] < 0. and lower[j+1][1] > 0.][0]

        # Remove them from R and add to output
        R_remove(u1)
        R_remove(u2)
        P_append(u1)
        Q_append(u2)

        # Prevent certain edge cases
        if u1 != v1 and v2 != u2:
            R_remove(v1)
            R_remove(v2)
            P_append(v1)
            Q_append(v2)
    
    # Report
    return P,Q


######################################################################################################
# Aux

def randpoints(n, hw=1000., hh=1000.):
    """Returns lists of n (total) randomly generated points within the given planar region."""
    P = [(i, uniform(-hw,-hw/4), uniform(-hh,hh)) for i in range(n//2)]
    Q = [(i, uniform(hw/4, hw), uniform(-hh,hh)) for i in range(n//2,n//1)]
    
    return P,Q


######################################################################################################
# Test

if __name__ == '__main__':

    # Init
    n = 16
    P,Q = randpoints(n)

    # Do the thing
    # ctr = perf_counter()
    
    if n <= 1000:
        P,Q = chull(P+Q)    # Works faster for smaller sizes and hand-made examples
    else:
        P,Q = hscut(P,Q)    # Can process more than to 30000 points in less than 5 seconds
    
    # print(perf_counter() - ctr)

    # NOTE: comment/uncomment plotting/timer for high n

    fig,ax = plt.subplots()

    for i in range(len(P)):
        ax.plot((P[i][1], Q[i][1]), (P[i][2], Q[i][2]), 'k-')
        ax.plot((P[i][1]), (P[i][2]), 'ro')
        ax.plot((Q[i][1]), (Q[i][2]), 'bo')

    plt.show()
