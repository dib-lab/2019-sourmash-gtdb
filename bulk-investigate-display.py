#! /usr/bin/env python
"""
Summarize and display output from bulk-investigate.py
"""
from pickle import load
from sourmash.lca import lca_utils
import math
import pprint
from collections import defaultdict
import argparse

def main():
    p = argparse.ArgumentParser()
    p.add_argument('pickle_file', help='pickle file from bulk-investigate')
    args = p.parse_args()

    with open(args.pickle_file, 'rb') as fp:
        combo_counts = load(fp)

    # now figure out who the key class/etc players are: sort by number
    # of times they show up.
    combo_counts_items = list(combo_counts.items())
    combo_counts_items.sort(key = lambda x: -len(x[1]))

    total = 0
    for n, (blame_lineages, track_lineages) in enumerate(combo_counts_items):
        total += len(track_lineages)

    print('total lineages:', total)
        
    print('showing top 10 clusters:')
    print('')
    for n, (blame_lineages, track_lineages) in enumerate(combo_counts_items):
        if n == 10:
            break

        count = len(track_lineages)

        print('cluster {} showed up {} times'.format(n, count))

        x = []
        for signame, *rest in track_lineages:
            d = defaultdict(int)
            for lineages, count in rest:
                lineages = tuple(sorted(lineages))

                d[lineages] += count

            pairs = d.items()
            x.append((signame, pairs))

        for lin in blame_lineages:
            lintext = lca_utils.display_lineage(lin)
            print('   ', lintext)

        # print out some exemplary foo
        def calc_entropy(x):
            total = 0
            for lin, num in x[1:]:
                total += n

            H = 0.0
            if total:
                for lin, num in x[1:]:
                    H += -(num / total) * math.log(num / total, 2)

            return -H

        track_lineages.sort(key = lambda x: calc_entropy(x))
        for i in range(min(4, len(track_lineages))):
            print(' * lineage contributions #{}'.format(i + 1))
            signame = track_lineages[i][0]
            pairs = track_lineages[i][1:]
            print('    {} {:g}'.format(signame, calc_entropy(pairs)))
            for lin, count in pairs:
                print('     -', count, lca_utils.display_lineage(lin))


if __name__ == '__main__':
    main()
