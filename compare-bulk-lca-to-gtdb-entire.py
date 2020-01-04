#! /usr/bin/env python
import os
import csv
import sys
import argparse
import pprint
from collections import defaultdict

import sourmash
from sourmash.lca import lca_utils
from sourmash.lca.command_index import load_taxonomy_assignments


def main():
    p = argparse.ArgumentParser()
    p.add_argument('gtdb_lineages')
    p.add_argument('bulk_lineages')
    p.add_argument('-v', '--verbose', action='store_true')
    args = p.parse_args()

    gtdb_lineages, _ = load_taxonomy_assignments(args.gtdb_lineages, start_column=3)
    bulk_lineages, _ = load_taxonomy_assignments(args.bulk_lineages)

    print('gtdb only:', len(set(gtdb_lineages.keys()) - set(bulk_lineages.keys())))
    print('bulk only:', len(set(bulk_lineages.keys()) - set(gtdb_lineages.keys())))
    print('common:', len(set(gtdb_lineages.keys()).intersection(bulk_lineages.keys())))

    common = set(gtdb_lineages.keys()).intersection(bulk_lineages.keys())

    same = 0
    different = 0
    consistent = 0
    consistent_at = defaultdict(int)
    disagree = 0
    disagree_at = defaultdict(int)

    fp = open('superkingdom-and-phylum-list.txt', 'wt')
    
    for k in common:
        gtdb = gtdb_lineages[k]
        bulk = bulk_lineages[k]

        if gtdb == bulk:
            same += 1
        else:
            different += 1

            tree = lca_utils.build_tree([gtdb, bulk])
            lca, reason = lca_utils.find_lca(tree)
            if lca == gtdb:
                consistent += 1
                rank = 'root'
                if bulk:
                    rank = bulk[-1].rank
                consistent_at[rank] += 1

                if rank in ('superkingdom', 'phylum'):
                    fp.write('{}\n'.format(k))

            else:
                disagree += 1

                rank = 'root'
                if lca:
                    rank = lca[-1].rank
                disagree_at[rank] += 1
#                print(lca_utils.display_lineage(gtdb))
#                print(lca_utils.display_lineage(bulk))
#                print(lca_utils.display_lineage(lca))
                

        assert gtdb[-1].rank == 'species', gtdb[-1].rank


    print('same:', same)
    print('different:', different)
    print('')
    print('different but consistent:', consistent)
    for rank in lca_utils.taxlist(include_strain=False):
        print('   rank: {} / count: {}'.format(rank, consistent_at[rank]))
    print('')
    print('inconsistent:', disagree)

    for rank in lca_utils.taxlist(include_strain=False):
        print('   rank: {} / count: {}'.format(rank, disagree_at[rank]))


if __name__ == '__main__':
    main()
