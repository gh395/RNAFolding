from collections import deque
import visuals
import argparse

"""
Compute a RNA according to Nussinov algorithm (maximal base pair)

Arguments:
    -s: sequence to read in (pasted in terminal)
Outputs:
    list of base pairings in RNA (0 indexed)
    arc diagram
Example Usage:
    python nussinov.py -s <RNA sequence>
DOES NOT HANDLE NON RNA SEQUENCES (no ambigious placeholders)
"""


def delta(sequence, i, j):
    base1 = sequence[i]
    base2 = sequence[j]
    return (base1 == "A" and base2 == "U") or (base1 == "U" and base2 == "A") or \
        (base1 == "G" and base2 == "C") or (base1 == "C" and base2 == "G") or \
        (base1 == "G" and base2 == "U") or (base1 == "U" and base2 == "G")


def traceback(dp, sequence):
    n = len(sequence)
    stack = deque()
    stack.append((0, n-1))
    pairs = []
    while stack:
        f = stack.pop()
        i, j = f
        if i >= j:
            continue
        elif dp[i+1][j] == dp[i][j]:
            stack.append((i+1, j))
        elif dp[i][j-1] == dp[i][j]:
            stack.append((i, j-1))
        elif dp[i+1][j-1]+delta(sequence, i, j) == dp[i][j]:
            pairs.append((i, j))
            stack.append((i+1, j-1))
        else:
            for k in range(i+1, j-1):
                if dp[i][k]+dp[k+1][j] == dp[i][j]:
                    stack.append((k+1, j))
                    stack.append((i, k))
                    break
    return pairs


def nussinov(sequence):
    n = len(sequence)
    dp = [[None]*n for _ in range(n)]
    for i in range(0, n):
        dp[i][i] = 0
    for i in range(1, n):
        dp[i][i-1] = 0
    for k in range(1, n):
        i = 0
        j = k
        while j < n:
            x = max(
                dp[i+1][j],
                dp[i][j-1],
                dp[i+1][j-1]+delta(sequence, i, j),
                max((dp[i][k]+dp[k+1][j]) for k in range(i, j))
            )
            dp[i][j] = x
            i += 1
            j += 1
    return traceback(dp, sequence)


def main():
    parser = argparse.ArgumentParser(
        description='Compute RNA fold using Nussinov algorithm')
    parser.add_argument('-s', action="store", dest="s", type=str)
    args = parser.parse_args()
    sol = nussinov(args.s)
    print(sol)
    print(visuals.vis_brak(sol, args.s))
    visuals.vis_arc(sol, args.s)


if __name__ == '__main__':
    main()
