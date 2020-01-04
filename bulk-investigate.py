#! /usr/bin/env python
"""
Dig into results from bulk-classify-sbt-with-lca.py
"""
import os
import argparse
import sys
from collections import defaultdict
import pprint
import math
from pickle import dump

from sourmash import load_signatures
from sourmash.logging import error, notify, set_quiet
from sourmash.lca import lca_utils
from sourmash.lca.command_summarize import summarize
from sourmash.lca.command_classify import classify_signature
from sourmash import sourmash_args

DEFAULT_THRESHOLD=5


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--db', nargs='+', action='append')
    p.add_argument('--query', nargs='+', action='append')
    p.add_argument('--threshold', type=int, default=DEFAULT_THRESHOLD)
    p.add_argument('--traverse-directory', action='store_true',
                        help='load all signatures underneath directories.')
    p.add_argument('-o', '--output', help='output pickle file name')
    p.add_argument('--scaled', type=float)
    p.add_argument('-q', '--quiet', action='store_true',
                   help='suppress non-error output')
    p.add_argument('-d', '--debug', action='store_true',
                   help='output debugging output')

    args = p.parse_args()

    if not args.output:
        error("must supply -o/--output pickle file")
        sys.exit(-1)

    if not args.db:
        error('Error! must specify at least one LCA database with --db')
        sys.exit(-1)

    if not args.query:
        error('Error! must specify at least one query signature with --query')
        sys.exit(-1)

    set_quiet(args.quiet, args.debug)

    if args.scaled:
        args.scaled = int(args.scaled)

    # flatten --db and --query
    args.db = [item for sublist in args.db for item in sublist]
    args.query = [item for sublist in args.query for item in sublist]

    # load all the databases
    dblist, ksize, scaled = lca_utils.load_databases(args.db, args.scaled)

    # find all the queries
    notify('finding query signatures...')
    if args.traverse_directory:
        inp_files = list(sourmash_args.traverse_find_sigs(args.query))
    else:
        inp_files = list(args.query)

    combo_counts = defaultdict(list)

    # for each query, gather all the hashvals across databases
    total_count = 0
    n = 0
    total_n = len(inp_files)
    for n, query_filename in enumerate(inp_files):
        if n and n % 100 == 0:
            print('...', n)

        hashvals = defaultdict(int)
        n += 1
        for query_sig in load_signatures(query_filename, ksize=ksize):
            total_count += 1

            mh = query_sig.minhash.downsample_scaled(scaled)
            for hashval in mh.get_mins():
                hashvals[hashval] += 1

            # get the full counted list of lineage counts in this signature
            lineage_counts = summarize(hashvals, dblist, args.threshold)

            if not lineage_counts:
                continue

            # also, separately, classify the signature, to get the lca:
            lineage, status = classify_signature(query_sig, dblist,
                                                 args.threshold)

            # figure out the rank-after-classify => that's where it's confusing
            lca_rank = 'root'
            next_rank = 'superkingdom'
            if lineage:
                lca_rank = lineage[-1].rank
                if lca_rank == 'superkingdom':
                    next_rank = 'phylum'
                elif lca_rank == 'phylum':
                    next_rank = 'class'

            print('---\nassigned at {} -- {}'.format(lca_rank, query_filename))
            blame_lineages = []
            track_lineages = [query_filename]
            for (lineage, count) in lineage_counts.items():
                this_rank = 'root'
                if lineage:
                    this_rank = lineage[-1].rank

                if lca_rank in ('root', 'superkingdom', 'phylum') and \
                   next_rank == this_rank:
                    print(lca_utils.display_lineage(lineage), count)
                    blame_lineages.append(lineage)
                    track_lineages.append((lineage, count))
                    continue

            if blame_lineages:
                blame_lineages.sort()
                blame_lineages = tuple(blame_lineages)

                combo_counts[blame_lineages].append(track_lineages)
                


    with open(args.output, 'wb') as fp:
        dump(combo_counts, fp)

    sys.exit(0)
    

if __name__ == '__main__':
    main()
