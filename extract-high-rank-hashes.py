#! /usr/bin/env python
"""
Look for compositional + taxonomic oddities in an LCA database.
"""
import sourmash
import sys
from collections import defaultdict
import pprint
import argparse

from sourmash.logging import error, debug, set_quiet, notify
from sourmash.lca import lca_utils
from sourmash.sourmash_args import SourmashArgumentParser


def make_lca_counts(dblist):
    """
    Collect counts of all the LCAs in the list of databases.
    """

    # gather all hashvalue assignments from across all the databases
    assignments = defaultdict(set)
    for lca_db in dblist:
        for hashval, idx_list in lca_db.hashval_to_idx.items():
            for idx in idx_list:
                lid = lca_db.idx_to_lid.get(idx)
                if lid is not None:
                    lineage = lca_db.lid_to_lineage[lid]
                    assignments[hashval].add(lineage)

    # now convert to trees -> do LCA for each hashval
    counts = defaultdict(int)
    crossdict = defaultdict(set)
    for hashval, lineages in assignments.items():

        # for each list of tuple_info [(rank, name), ...] build
        # a tree that lets us discover lowest-common-ancestor.
        debug('lineages: {}', lineages)
        tree = lca_utils.build_tree(lineages)

        # now find either a leaf or the first node with multiple
        # children; that's our lowest-common-ancestor node.
        lca, reason = lca_utils.find_lca(tree)

        if lca:
            rank = lca[-1].rank
            crossdict[rank].add(hashval)
        else:
            crossdict['root'].add(hashval)

    return crossdict


def main(args):
    p = argparse.ArgumentParser(prog="sourmash lca rankinfo")
    p.add_argument('db', nargs='+')
    p.add_argument('--scaled', type=float)
    p.add_argument('-q', '--quiet', action='store_true',
                   help='suppress non-error output')
    p.add_argument('-d', '--debug', action='store_true',
                   help='output debugging output')
    p.add_argument('-o', '--output', type=str, help='output filename')
    p.add_argument('--lowest-rank', default='phylum')
    args = p.parse_args(args)

    if not args.db:
        error('Error! must specify at least one LCA database with --db')
        sys.exit(-1)

    set_quiet(args.quiet, args.debug)

    if args.scaled:
        args.scaled = int(args.scaled)

    if not args.output:
        error('Please specify an output filename')
        sys.exit(-1)

    keep_ranks = ['root']
    for rank in lca_utils.taxlist():
        keep_ranks.append(rank)
        if rank == args.lowest_rank:
            break

    print('outputting hashvals at following ranks:', keep_ranks)

    # load all the databases
    dblist, ksize, scaled = lca_utils.load_databases(args.db, args.scaled)
    assert len(dblist) == 1

    # count all the LCAs across these databases
    crossdict = make_lca_counts(dblist)

    # output basic stats
    for rank, v in crossdict.items():
        print(rank, len(v))

    n = 0
    all_hashvals = set()
    if args.output:
        with open(args.output, 'wt') as fp:
            for rank in keep_ranks:
                for hashval in crossdict[rank]:
                    fp.write("{}\n".format(hashval))
                    n += 1
                all_hashvals.update(crossdict[rank])
    else:
        assert 0

    total = sum([ len(v) for v in crossdict.values() ])
    print('wrote {} confused hashvals, of {} total'.format(n, total))


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
