#! /usr/bin/env python
"""
Look for compositional oddities that do not match ANI in an LCA.
(Ignore taxonomy except for reporting.)
"""
import sourmash
import sys
from collections import defaultdict

from sourmash.logging import error, debug, set_quiet, notify
from sourmash.lca import lca_utils
from sourmash.sourmash_args import SourmashArgumentParser


def make_assignment_counts(lca_db, min_num=0):
    """
    Collect counts of all the lineages in the database

    CTB this could usefully be converted to a generator function.
    """

    # gather all hashvalue assignments from across all the databases
    assignments = defaultdict(set)
    for hashval, idx_list in lca_db.hashval_to_idx.items():
        if min_num and len(idx_list) < min_num:
            continue

        for idx in idx_list:
            assignments[hashval].add(idx)

    return assignments


def main(args):
    """
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
    lca_db = dblist[0]

    # count all the LCAs across these databases
    counts = make_assignment_counts(lca_db, args.minimum_num)

    idx_groups = defaultdict(int)
    for hashval, idx_list in counts.items():
        if len(idx_list) >= 2:
            idx_list = list(idx_list)
            idx_list.sort()
            idx_list = tuple(idx_list)
            idx_groups[idx_list] += 1

    n = 0
    sigd = lca_db._signatures
    for idx_group, count in idx_groups.items():
        if count >= 5:
            keep = False
            for idx in idx_group:
                mh1 = sigd[idx]
                for idx2 in idx_group:
                    mh2 = sigd[idx2]

                    print(mh1.similarity(mh2))
                    if mh1.similarity(mh2) < 0.01:
                        keep = True

            if keep:
                print('group:')
                for idx in idx_group:
                    print('* ', lca_db.idx_to_ident[idx])
                    lid = lca_db.idx_to_lid[idx]
                    print('  ', ";".join(lca_utils.zip_lineage(lca_db.lid_to_lineage[lid])))
                print('')
                n += 1

    print(n)


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
