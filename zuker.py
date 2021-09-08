from collections import deque
import numpy as np
from tkinter import *
import nussinov
import visuals
import argparse

"""
Compute a RNA according to Zuker algorithm (min free energy)

Arguments:
    -s: sequence to read in (pasted in terminal)
Outputs:
    list of base pairings in RNA (0 indexed)
    dot bracket notation
    arc diagram
Example Usage:
    python zuker.py -s <RNA sequence>
DOES NOT HANDLE NON RNA SEQUENCES (no ambigious placeholders)
"""

R = 1.987*(10**(-3))  # gas constant kcal/
T = 310  # temp in kelvin (body temp)

# all data from Turner (2004) unless otherwise specified

# stacking energy of base rings
stack_energy = {  # (i,j) (i+1,j-1) stack
    "AU": {
        "AU": -0.9, "CG": -2.2, "GC": -2.1, "GU": -0.6, "UA": -1.1, "UG": -1.4
    }, "CG": {
        "AU": -2.1, "CG": -3.3, "GC": -2.4, "GU": -1.4, "UA": -2.1, "UG": -2.1
    }, "GC": {
        "AU": -2.4, "CG": -3.4, "GC": -3.3, "GU": -1.5, "UA": -2.2, "UG": -2.5
    }, "GU": {
        "AU": -1.3, "CG": -2.5, "GC": -2.1, "GU": -0.5, "UA": -1.4, "UG": 1.3
    }, "UA": {
        "AU": -1.3, "CG": -2.4, "GC": -2.1, "GU": -1.0, "UA": -0.9, "UG": -1.3
    }, "UG": {
        "AU": -1.0, "CG": -1.5, "GC": -1.4, "GU": 0.3, "UA": -0.6, "UG": -0.5
    }
}

# loop energy initialization (destabilizing)
# [internal, bulge, hairpin]
loop_energy = [  # each row is the size of the loop
    [0, 0, 0],  # 0
    [0, 3.8, 0],  # 1 (null values replaced with 0)
    [0, 2.8, 0],  # 2
    [0, 3.2, 5.4],  # 3
    [1.1, 3.6, 5.6],  # 4
    [2.0, 4.0, 5.7],  # 5
    [2.0, 4.4, 5.4],  # 6
    [2.1, 4.6, 6.0],  # 7
    [2.3, 4.7, 5.5],  # 8
    [2.4, 4.8, 6.4],  # 9
    [2.5, 4.9, 6.5],  # 10
    [2.6, 5.0, 6.6],  # 11
    [2.7, 5.1, 6.7],  # 12
    [2.8, 5.2, 6.8],  # 13
    [2.9, 5.3, 6.9],  # 14
    [2.9, 5.4, 6.9],  # 15
    [3.0, 5.4, 7.0],  # 16
    [3.1, 5.5, 7.1],  # 17
    [3.1, 5.5, 7.1],  # 18
    [3.2, 5.6, 7.2],  # 19
    [3.3, 5.7, 7.2],  # 20
    [3.3, 5.7, 7.3],  # 21
    [3.4, 5.8, 7.3],  # 22
    [3.4, 5.8, 7.4],  # 23
    [3.5, 5.8, 7.4],  # 24
    [3.5, 5.9, 7.5],  # 25
    [3.5, 5.9, 7.5],  # 26
    [3.6, 6.0, 7.5],  # 27
    [3.6, 6.0, 7.6],  # 28
    [3.7, 6.0, 7.6],  # 29
    [3.7, 6.1, 7.7]  # 30
]

# terminal mismatches are non-canonical pairs next to *helix* ends
terminal_mismatch = {
    "AU": {
        "AA": -0.8, "AC": -1.0, "AG": -0.8, "AU": -1.0,
        "CA": -0.6, "CC": -0.7, "CG": -0.6, "CU": -0.7,
        "GA": -0.8, "GC": -1.0, "GG": -0.8, "GU": -1.0,
        "UA": -0.6, "UC": -0.8, "UG": -0.6, "UU": -0.8,
    }, "CG": {
        "AA": -1.5, "AC": -1.5, "AG": -1.4, "AU": -1.5,
        "CA": -1.0, "CC": -1.1, "CG": -1.0, "CU": -0.8,
        "GA": -1.4, "GC": -1.5, "GG": -1.6, "GU": -1.5,
        "UA": -1.0, "UC": -1.4, "UG": -1.0, "UU": -1.2,
    }, "GC": {
        "AA": -1.1, "AC": -1.5, "AG": -1.3, "AU": -1.5,
        "CA": -1.1, "CC": -0.7, "CG": -1.1, "CU": -0.5,
        "GA": -1.6, "GC": -1.5, "GG": -1.4, "GU": -1.5,
        "UA": -1.1, "UC": -1.0, "UG": -1.1, "UU": -0.7,
    }, "GU": {
        "AA": -0.3, "AC": -1.0, "AG": -0.8, "AU": -1.0,
        "CA": -0.6, "CC": -0.7, "CG": -0.6, "CU": -0.7,
        "GA": -0.6, "GC": -1.0, "GG": -0.8, "GU": -1.0,
        "UA": -0.6, "UC": -0.8, "UG": -0.6, "UU": -0.6,
    }, "UA": {
        "AA": -1.0, "AC": -0.8, "AG": -1.1, "AU": -0.8,
        "CA": -0.7, "CC": -0.6, "CG": -0.7, "CU": -0.5,
        "GA": -1.1, "GC": -0.8, "GG": -1.2, "GU": -0.8,
        "UA": -0.7, "UC": -0.6, "UG": -0.7, "UU": -0.5,
    }, "UG": {
        "AA": -1.0, "AC": -0.8, "AG": -1.1, "AU": -0.8,
        "CA": -0.7, "CC": -0.6, "CG": -0.7, "CU": -0.5,
        "GA": -0.5, "GC": -0.8, "GG": -0.8, "GU": -0.8,
        "UA": -0.7, "UC": -0.6, "UG": -0.7, "UU": -0.5,
    }
}

# special hairpins
sp_hps = {
    "CAACG": 6.8, "GUUAC": 6.9, "CUACGG": 2.8, "CUCCGG": 2.7, "CUUCGG": 3.7,
    "CUUUGG": 3.7, "CCAAGG": 3.3, "CCCAGG": 3.4, "CCGAGG": 3.5, "CCUAGG": 3.7,
    "CCACGG": 3.7, "CCGCGG": 3.6, "CCUCGG": 2.5, "CUAAGG": 3.6, "CUCAGG": 3.7,
    "CUUAGG": 3.5, "CUGCGG": 2.8, "CAACGG": 5.5, "ACAGUGCU": 2.9,
    "ACAGUGAU": 3.6, "ACAGUGUU": 1.8, "ACAGUACU": 2.8
}

tetraloops = {'GGGGAC': -3.0, 'GGUGAC': -3.0, 'CGAAAG': -3.0, 'GGAGAC': -3.0,
              'CGCAAG': -3.0, 'GGAAAC': -3.0, 'CGGAAG': -3.0, 'CUUCGG': -3.0,
              'CGUGAG': -3.0, 'CGAAGG': -2.5, 'CUACGG': -2.5, 'GGCAAC': -2.5,
              'CGCGAG': -2.5, 'UGAGAG': -2.5, 'CGAGAG': -2.0, 'AGAAAU': -2.0,
              'CGUAAG': -2.0, 'CUAACG': -2.0, 'UGAAAG': -2.0, 'GGAAGC': -1.5,
              'GGGAAC': -1.5, 'UGAAAA': -1.5, 'AGCAAU': -1.5, 'AGUAAU': -1.5,
              'CGGGAG': -1.5, 'AGUGAU': -1.5, 'GGCGAC': -1.5, 'GGGAGC': -1.5,
              'GUGAAC': -1.5, 'UGGAAA': -1.5}

# free energy of stacking bp, depends on all four bases in a stack


def eS(i, j, sequence):
    # ignore symmetry penalty
    length = abs(i-j)-1
    if length < 2:
        return float("INF")
    b1 = sequence[i]
    b2 = sequence[j]
    b3 = sequence[i+1]
    b4 = sequence[j-1]
    p1 = b1+b2
    p2 = b3+b4
    matches = {'AU', 'CG', 'GC', 'GU', 'UA', 'UG'}

    if p1 not in matches or p2 not in matches:
        return 0

    # dG helix inititation = 4.09
    return stack_energy[p1][p2]

# free energy of the hairpin closed by (i,j) bp
# depends on loopsize,unpaired bases adjacent to i,j in loop


def eH(i, j, sequence):
    # ignore penalty for C loop
    length = abs(i-j)-1
    if length < 3:
        return float("INF")

    b1 = sequence[i]
    b2 = sequence[j]
    b3 = sequence[i+1]
    b4 = sequence[j-1]
    # UU/GA first mismatch (not AG) added stability, GG first mismatch added stability
    uu_gaMM = 0
    ggMM = 0
    if (b3 == 'U' and b4 == 'U') or (b3 == 'G' and b4 == 'A'):
        uu_gaMM = -0.9
    if (b3 == 'G' and b4 == 'G'):
        ggMM = -0.8
    # GU-GG closure gets added stability
    gu_closure = 0
    if (b1 == 'G' and b2 == 'U') or (b2 == 'G' and b1 == 'U'):
        if (b3 == 'G' and b4 == 'G'):
            gu_closure = -2.2

    # terminal mismatch
    p1 = b1+b2  # pair 1
    p2 = b3+b4  # pair 2
    tm = terminal_mismatch[p1][p2]

    special = False
    if length == 3 or length == 4 or length == 6:
        hp = sequence[i:j+1]
        if hp in sp_hps.keys():
            lp_energy = sp_hps[hp]
            special = True

    tl_bonus = 0
    ''' #tends to overfit for these tetraloops; data from Turner 1999 (not 2004)
    if length == 4:
        tl = sequence[i:j+1]
        if tl in tetraloops.keys():
            tl_bonus = tetraloops[tl]
    '''
    if special:
        return lp_energy+gu_closure+uu_gaMM+ggMM+tm+tl_bonus
    elif length <= 30:
        return loop_energy[length][2]+gu_closure+uu_gaMM+ggMM+tm+tl_bonus
    # dG(n>6)=dG(6)+1.75RTln(n/6)
    return loop_energy[9][2]+1.75*R*T*np.log(length/9)+gu_closure+uu_gaMM+ggMM+tm+tl_bonus

# bulge/internal loop energy


def eL(i, j, p, q, sequence):
    i_bulge = p-i+1
    j_bulge = j-q+1
    if i_bulge == 0 or j_bulge == 0:  # bulge
        if i_bulge == 0:
            bulge = j_bulge
        else:
            bulge = i_bulge
        if bulge <= 30:
            return loop_energy[bulge][1]
        # dG(n>6)=dG(6)+1.75RTln(n/6)
        return loop_energy[6][1]+1.75*R*T*np.log(bulge/6)
    else:  # internal loop
        b1 = sequence[i]
        b2 = sequence[j]
        b3 = sequence[p]
        b4 = sequence[q]
        p1 = b1+b2  # pair 1
        p2 = b3+b4  # pair 2
        penalty = 0
        # AU/GU closure penalty
        if p1 == 'AU' or p1 == 'UA' or p1 == 'GU' or p1 == 'UG':
            penalty += 0.65
        if p2 == 'AU' or p2 == 'UA' or p2 == 'GU' or p2 == 'UG':
            penalty += 0.65

        int_loop = i_bulge+j_bulge
        asymmetry = abs(i_bulge-j_bulge)
        asym_penalty = 0.48*asymmetry
        # print(int_loop)
        if int_loop <= 30:
            return loop_energy[int_loop][0]+asym_penalty

        # dG(n>6)= dG(6)+1.08×ln(n/6)
        return loop_energy[6][0]+1.08*np.log(int_loop/6)+asym_penalty


def traceback(W, V, sequence):
    n = len(sequence)
    stack = deque()
    stack.append((0, n-1))
    pairs = []
    while stack:
        f = stack.pop()
        i, j = f
        if i >= j:
            continue
        elif W[i+1][j] == W[i][j]:
            stack.append((i+1, j))
        elif W[i][j-1] == W[i][j]:
            stack.append((i, j-1))
        elif V[i][j] == W[i][j]:
            pairs.append((i, j))
            stack.append((i+1, j-1))
        else:
            for k in range(i+1, j-1):
                if W[i][k]+W[k+1][j] == W[i][j]:
                    stack.append((k+1, j))
                    stack.append((i, k))
                    break
    return pairs


def zuker(sequence):
    n = len(sequence)
    W = [[0]*n for _ in range(n)]
    V = [[float("INF")]*n for _ in range(n)]
    for k in range(1, n):
        i = 0
        j = k
        while j < n:
            if not nussinov.delta(sequence, i, j):
                V[i][j] = float("INF")
            else:
                if j-i-1 > 2:  # need >2 "inside bases" for these recurrence
                    eM = min(W[i+1][k] + W[k+1][j-1]
                             for k in range(i+2, j-1))
                    # double loop O(n^4), limit size of internal/bulge loops
                    # to restrict runtime. bulges limited to 30, so internal
                    # loops limited to 60.
                    eL1 = min((eL(i, j, p, q, sequence) + V[p][q]) for p in range(i, min(j, i+30))
                              for q in range(max(p, j-30), j))

                else:
                    eM = float("INF")
                    eL1 = float("INF")

                V[i][j] = min(
                    eH(i, j, sequence),  # hairpin  (E1)
                    eS(i, j, sequence) + V[i+1][j-1],  # stacking energy (E2)
                    eL1,  # bulge/internal loop (E2)

                    eM  # bifurcation loop (E3)
                )
            W[i][j] = min(
                W[i+1][j],  # i unpaired
                W[i][j-1],  # j unpaired
                V[i][j],  # ij pair
                min((W[i][k]+W[k+1][j]) for k in range(i, j))  # bifurcation
            )
            i += 1
            j += 1
    return traceback(W, V, sequence)


def main():
    parser = argparse.ArgumentParser(
        description='Compute RNA fold using Zuker algorithm')
    parser.add_argument('-s', action="store", dest="s", type=str)
    args = parser.parse_args()
    sol = zuker(args.s)
    print(sol)
    print(visuals.vis_brak(sol, args.s))
    visuals.vis_arc(sol, args.s)


if __name__ == '__main__':
    main()
