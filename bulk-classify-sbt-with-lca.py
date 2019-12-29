#! /usr/bin/env python
"""
Try to classify everything in an SBT using an LCA database.

Uses `sourmash lca classify` algorithm on every signature in the SBT.

Outputs a spreadsheet of classifications. The first column indicates
the lowest taxonomic rank at which something is unambiguously classified,
and will be "MISSED" if no classification is available.

Also outputs a directory of unclassified signatures.
"""
import sourmash
import sys
from collections import defaultdict
import pprint
import csv
import os

from sourmash.logging import error, debug, set_quiet, notify
from sourmash.lca import lca_utils
from sourmash.lca.command_classify import classify_signature
import argparse


def main(args):
    """
    """
    p = argparse.ArgumentParser()
    p.add_argument('prefix')
    p.add_argument('lca_db', nargs='+')
    p.add_argument('sbt')
    p.add_argument('--scaled', type=float)
    p.add_argument('-q', '--quiet', action='store_true',
                   help='suppress non-error output')
    p.add_argument('-d', '--debug', action='store_true',
                   help='output debugging output')
    args = p.parse_args(args)

    dirname = '{}-unclassified-sigs'.format(args.prefix)
    try:
        print('making output sigs dir: {}'.format(dirname))
        os.mkdir(dirname)
    except:
        print('WARNING: {} already exists.'.format(dirname), file=sys.stderr)

    if not args.lca_db:
        error('Error! must specify at least one LCA database with')
        sys.exit(-1)

    set_quiet(args.quiet, args.debug)

    if args.scaled:
        args.scaled = int(args.scaled)

    # load all the databases
    dblist, ksize, scaled = lca_utils.load_databases(args.lca_db, args.scaled)

    print(ksize, scaled)

    sbt_db = sourmash.load_sbt_index(args.sbt)

    counts = defaultdict(int)
    n_missed = 0
    fp = open('{}-bulk-classify.csv'.format(args.prefix), 'wt')
    w = csv.writer(fp)
    w.writerow(["rank", "name", "filename", "md5sum", "lineage"])
    for n, sig in enumerate(sbt_db.signatures()):
        if n % 100 == 0:
            print('...', n)
            fp.flush()

        lineage = ''
        classified_as, why = classify_signature(sig, dblist, 5)

        if classified_as:
            rank = classified_as[-1].rank
            counts[rank] += 1
            lineage = lca_utils.display_lineage(classified_as)
        elif why == 'disagree':
            rank = 'root'
            lineage = ''
        else:
            rank = 'MISSED'
            n_missed += 1

        w.writerow([rank, sig.name(), sig.d['filename'], sig.md5sum(), lineage])

        if rank not in ('genus', 'species', 'family', 'order'):
            md5name = sig.md5sum()
            with open('{}/{}.sig'.format(dirname, md5name), 'wt') as fp2:
                sourmash.save_signatures([sig], fp2)

        if n % 1000 == 0 and n:
            print('at', n, 'genomes...')
            pprint.pprint(list(counts.items()))
            print('missed:', n_missed, 'of', n)
        
    pprint.pprint(list(counts.items()))
    print('missed:', n_missed, 'of', n)


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
