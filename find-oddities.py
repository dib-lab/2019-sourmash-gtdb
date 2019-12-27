#! /usr/bin/env python
"""
Look for compositional + taxonomic oddities in an LCA database.
"""
import sourmash
import sys
from collections import defaultdict
import argparse
import csv

from sourmash.logging import error, debug, set_quiet, notify
from sourmash.lca import lca_utils


def make_lca_counts(dblist, lowest_rank='phylum', min_num=0, min_hashes=5,
                    prefix='oddities'):
    """
    Collect counts of all the LCAs in the list of databases.
    """
    assert len(dblist) == 1

    keep_ranks = ['root']
    for rank in lca_utils.taxlist():
        keep_ranks.append(rank)
        if rank == lowest_rank:
            break
    print('keeping hashvals at following ranks:', keep_ranks)
    print('min number of lineages:', min_num)
    print('min number of shared hashes:', min_hashes)

    print('---')

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
        debug("{}", lineages)
        tree = lca_utils.build_tree(lineages)

        # now find either a leaf or the first node with multiple
        # children; that's our lowest-common-ancestor node.
        lca, reason = lca_utils.find_lca(tree)

        # find cross-superkingdom hashes, and record combinations of lineages
        # that have them.
        rank = 'root'
        if lca:
            rank = lca[-1].rank
        
        if rank in keep_ranks:
            xx = []
            for lineage in lineages:
                xx.append(tuple(lineage))
            xx = tuple(xx)

            mixdict[xx].add(hashval)

        counts[lca] += 1

    # sort on number of confused hash vals by combination of lineages.
    mixdict_items = list(mixdict.items())
    mixdict_items.sort(key = lambda x: -len(x[1]))

    confused_hashvals = set()

    fp = open(prefix + '.csv', 'wt')
    w = csv.writer(fp)
    w.writerow(['cluster', 'num_lineages', 'shared_kmers', 'ksize', 'rank',
                'lca', 'ident1', 'lineage1', 'ident2', 'lineage2'])

    # filter & display
    cluster_n = 0
    for n, (lineages, hashvals) in enumerate(mixdict_items):
        # insist on more than N hash vals
        if len(hashvals) < min_hashes:
            continue
        
        # display summary:
        print('cluster {} has {} assignments for {} hashvals / {} bp'.format(n, len(lineages), len(hashvals), dblist[0].scaled * len(hashvals)))
        confused_hashvals.update(hashvals)

        tree = lca_utils.build_tree(lineages)
        lca, reason = lca_utils.find_lca(tree)
        if lca:
            rank = lca[-1].rank
        else:
            rank = 'root'
        print('rank & lca:', rank, lca_utils.display_lineage(lca))

        all_idxs = []
        for lineage_n, lineage in enumerate(lineages):
            print('* ', lca_utils.display_lineage(lineage))

            lids = dblist[0].lineage_to_lids[lineage]
            for lid in lids:
                idxs = dblist[0].lid_to_idx[lid]
                all_idxs.extend(idxs)
                for idx in idxs:
                    ident = dblist[0].idx_to_ident[idx]
                    print('  ', ident)
        print('')

        for i in range(len(all_idxs)):
            idx1 = all_idxs[i]
            lid1 = dblist[0].idx_to_lid[idx1]
            lin1 = dblist[0].lid_to_lineage[lid1]
            for j in range(i):
                idx2 = all_idxs[j]
                lid2 = dblist[0].idx_to_lid[idx2]
                lin2 = dblist[0].lid_to_lineage[lid2]

                ident1 = dblist[0].idx_to_ident[idx1]
                ident2 = dblist[0].idx_to_ident[idx2]

                debug("{} x {}", ident1, ident2)

                this_tree = lca_utils.build_tree([lin1, lin2])
                this_lca, this_reason = lca_utils.find_lca(this_tree)

                if lca != this_lca:
                    #print('skip lca', lca, this_lca)
                    continue

                mh1 = dblist[0]._signatures[idx1]
                mh2 = dblist[0]._signatures[idx2]

                mins1 = set(mh1.get_mins())
                mins2 = set(mh2.get_mins())

                intersect_size = len(mins1.intersection(mins2))

                if intersect_size < min_hashes:
                    #print('skip', intersect_size, len(mh1), len(mh2))
                    continue

                w.writerow(['cluster{}'.format(cluster_n),
                            len(lineages),
                            intersect_size * dblist[0].scaled,
                            dblist[0].ksize,
                            rank,
                            lca_utils.display_lineage(lca),
                            ident1,
                            lca_utils.display_lineage(lin1),
                            ident2,
                            lca_utils.display_lineage(lin2)])

                cluster_n += 1

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
    p.add_argument('--minimum-hashes', type=int, default=5,
                   help='Minimum number of hashes lineages must share to be reported')
    p.add_argument('--lowest-rank', default='phylum')
    p.add_argument('--prefix', default=None, help='prefix for output files')
    args = p.parse_args(args)

    if not args.db:
        error('Error! must specify at least one LCA database with --db')
        sys.exit(-1)

    set_quiet(args.quiet, args.debug)

    if args.scaled:
        args.scaled = int(args.scaled)

    # load all the databases
    print('loading databases:', args.db)
    dblist, ksize, scaled = lca_utils.load_databases(args.db, args.scaled)
    assert len(dblist) == 1

    # count all the LCAs across these databases
    counts, confused_hashvals = make_lca_counts(dblist,
                                                lowest_rank=args.lowest_rank,
                                                min_num=args.minimum_num,
                                                min_hashes=args.minimum_hashes,
                                                prefix=args.prefix
    )

    with open('confused_hashvals.txt', 'wt') as fp:
        fp.write("\n".join([ str(i) for i in confused_hashvals ]))


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
