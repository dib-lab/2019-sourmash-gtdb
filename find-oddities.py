#! /usr/bin/env python
"""
Look for compositional + taxonomic oddities in an LCA database.
"""
import sourmash
import sys
from collections import defaultdict
import argparse

from sourmash.logging import error, debug, set_quiet, notify
from sourmash.lca import lca_utils


def make_lca_counts(dblist, min_num=0):
    """
    Collect counts of all the LCAs in the list of databases.
    """

    # gather all hashvalue assignments from across all the databases
    assignments = defaultdict(set)
    for lca_db in dblist:
        for hashval, idx_list in lca_db.hashval_to_idx.items():
            if min_num and len(idx_list) < min_num:
                continue

            for idx in idx_list:
                lid = lca_db.idx_to_lid.get(idx)
                if lid is not None:
                    lineage = lca_db.lid_to_lineage[lid]
                    assignments[hashval].add(lineage)

    # now convert to trees -> do LCA & counts
    counts = defaultdict(int)
    mixdict = defaultdict(set)
    for hashval, lineages in assignments.items():

        # for each list of tuple_info [(rank, name), ...] build
        # a tree that lets us discover lowest-common-ancestor.
        debug(lineages)
        tree = lca_utils.build_tree(lineages)

        # now find either a leaf or the first node with multiple
        # children; that's our lowest-common-ancestor node.
        lca, reason = lca_utils.find_lca(tree)

        # find cross-superkingdom hashes, and record combinations of lineages
        # that have them.
        if not len(lca) or lca[-1].rank == 'superkingdom':
            xx = []
            for lineage in lineages:
                xx.append(tuple(lineage))
            xx = tuple(xx)

            mixdict[xx].add(hashval)

        counts[lca] += 1

    mixdict_items = list(mixdict.items())
    mixdict_items.sort(key = lambda x: -len(x[1]))

    confused_hashvals = set()

    # filter & display
    for n, (k, v) in enumerate(mixdict_items):
        # insist on more than 5 hash vals
        if len(v) > 5:
            # display:
            print('cluster {} has {} assignments for {} hashvals / {} bp'.format(n, len(k), len(v), dblist[0].scaled * len(v)))
            confused_hashvals.update(v)
            for lineage in k:
                print('* ', "; ".join(lca_utils.zip_lineage(lineage)))

                lids = dblist[0].lineage_to_lids[lineage]
                for lid in lids:
                    idxs = dblist[0].lid_to_idx[lid]
                    for idx in idxs:
                        ident = dblist[0].idx_to_ident[idx]
                        print('  ', ident)
            print('')

    return counts, confused_hashvals


def main(args):
    p = argparse.ArgumentParser()
    p.add_argument('db', nargs='+')
    p.add_argument('--scaled', type=float)
    p.add_argument('-q', '--quiet', action='store_true',
                   help='suppress non-error output')
    p.add_argument('-d', '--debug', action='store_true',
                   help='output debugging output')
    p.add_argument('--minimum-num', type=int, default=0,
                   help='Minimum number of different lineages a k-mer must be in to be counted')
    args = p.parse_args(args)

    if not args.db:
        error('Error! must specify at least one LCA database with --db')
        sys.exit(-1)

    set_quiet(args.quiet, args.debug)

    if args.scaled:
        args.scaled = int(args.scaled)

    # load all the databases
    dblist, ksize, scaled = lca_utils.load_databases(args.db, args.scaled)
    assert len(dblist) == 1

    # count all the LCAs across these databases
    counts, confused_hashvals = make_lca_counts(dblist, args.minimum_num)

    with open('confused_hashvals.txt', 'wt') as fp:
        fp.write("\n".join([ str(i) for i in confused_hashvals ]))


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
