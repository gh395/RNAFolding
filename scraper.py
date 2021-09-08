import argparse
import visuals
from sortedcontainers import SortedList

"""
Gets the RNA sequence information of a .ct file from RNAStrAlign

Arguments:
    -f: a compatible .ct file 
Outputs:
    the RNA sequence
    the true base pairings in a fold
    dot bracket notation (does not distinguish between pseudoknots)
    WILL OUTPUT A WARNING WHEN:
    -there are non RNA bases e.g. N placeholder
    -there is a pseudoknot (for your convience/info)
Example Usage:
    python zuker.py -f <file_path_to_ct_file>
DOES NOT HANDLE NON RNA SEQUENCES (no ambigious placeholders)
"""


def read_rnafile(filename):
    with open(filename, "r") as f:
        s = []
        seq = ""
        for l in f.readlines()[1:]:
            s.append(l.strip().split()[-2])
            char = l.strip().split()[1]
            if char not in ['A', 'G', 'C', 'U']:
                print("WARNING: contains non RNA bases")
            seq += char

        s2 = []
        in_pair = set()
        unclosed = set()
        matchings = dict()
        pk = False
        for count, pair in enumerate(s):
            p = int(pair)-1
            if pair != "0":
                if count not in in_pair:
                    s2.append((count, p))
                    unclosed.add(count)
                    unclosed.add(p)
                    matchings[count] = p
                if p in in_pair:
                    unclosed.remove(count)
                    unclosed.remove(p)
                else:
                    potential_pk = []
                    for u in unclosed:
                        if count > u:
                            potential_pk.append(u)
                    for u in potential_pk:
                        if matchings[u] < p:
                            pk = True

                in_pair.add(count)

                in_pair.add(p)
        if pk:
            print("WARNING: THERE IS AT LEAST ONE PSEUDOKNOT")
        return s2, seq


def main():
    parser = argparse.ArgumentParser(
        description='Scrape RNA database file for fold')
    parser.add_argument('-f', action="store", dest="f", type=str)
    args = parser.parse_args()
    output, seq = read_rnafile(args.f)
    print(seq)
    print(output)
    vb = visuals.vis_brak(output, seq)
    print(vb)
    visuals.vis_arc(output, seq)


if __name__ == '__main__':
    main()
