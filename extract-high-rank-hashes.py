#! /usr/bin/env python
"""
Look for compositional + taxonomic oddities in an LCA database.
"""
import sourmash
import sys
from collections import defaultdict

from sourmash.logging import error, debug, set_quiet, notify
from sourmash.lca import lca_utils
from sourmash.sourmash_args import SourmashArgumentParser


def make_lca_counts(dblist, min_num=0):
    """
    Collect counts of all the LCAs in the list of databases.

    CTB this could usefully be converted to a generator function.
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
    crossdict = defaultdict(set)
    for hashval, lineages in assignments.items():

        # for each list of tuple_info [(rank, name), ...] build
        # a tree that lets us discover lowest-common-ancestor.
        debug(lineages)
        tree = lca_utils.build_tree(lineages)

        # now find either a leaf or the first node with multiple
        # children; that's our lowest-common-ancestor node.
        lca, reason = lca_utils.find_lca(tree)
        counts[lca] += 1
        if lca:
            crossdict[lca[-1].rank].add(hashval)
        else:
            crossdict['root'].add(hashval)

    return crossdict

def rankinfo_main(args):
    """
    rankinfo!
    """
    p = SourmashArgumentParser(prog="sourmash lca rankinfo")
    p.add_argument('db', nargs='+')
    p.add_argument('--scaled', type=float)
    p.add_argument('-q', '--quiet', action='store_true',
                   help='suppress non-error output')
    p.add_argument('-d', '--debug', action='store_true',
                   help='output debugging output')
    p.add_argument('--minimum-num', type=int, default=0,
                   help='Minimum number of different lineages a k-mer must be in to be counted')
    p.add_argument('--output', type=str)
    args = p.parse_args(args)

    if not args.db:
        error('Error! must specify at least one LCA database with --db')
        sys.exit(-1)

    set_quiet(args.quiet, args.debug)

    if args.scaled:
        args.scaled = int(args.scaled)

    # load all the databases
    dblist, ksize, scaled = lca_utils.load_databases(args.db, args.scaled)

    # count all the LCAs across these databases
    crossdict = make_lca_counts(dblist, args.minimum_num)

    for rank, v in crossdict.items():
        print(rank, len(v))

    n = 0
    if args.output:
        with open(args.output, 'wt') as fp:
            for rank in ('root', 'superkingdom', 'class', 'phylum', 'order'):
                for hashval in crossdict[rank]:
                    fp.write("{}\n".format(hashval))
                    n += 1

    total = sum([ len(v) for v in crossdict.values() ])
    print('wrote {} confused hashvals, of {} total'.format(n, total))


if __name__ == '__main__':
    sys.exit(rankinfo_main(sys.argv[1:]))
