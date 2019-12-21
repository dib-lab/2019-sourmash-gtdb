#! /usr/bin/env python
"""
Dig into the bulk classification results from bulk-classify-sbt-with-lca.py.

Briefly, this script sorts results classed above order level into two bins
* those that contain multiple hashes belonging unambiguously to two different
  species, i.e. chimerae.
* other, e.g. things classed as a single lineage.
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
from sourmash.lca.command_summarize import summarize
import argparse

FILTER_AT='order'


def main(args):
    """
    """
    p = argparse.ArgumentParser()
    p.add_argument('prefix')
    p.add_argument('lca_db', nargs='+')
    p.add_argument('classify_csv')
    p.add_argument('--scaled', type=float)
    p.add_argument('--threshold', type=int, default=5)
    p.add_argument('-q', '--quiet', action='store_true',
                   help='suppress non-error output')
    p.add_argument('-d', '--debug', action='store_true',
                   help='output debugging output')
    args = p.parse_args(args)

    dirname = '{}-unclassified-sigs'.format(args.prefix)
    dirname2 = '{}-unclassified-sigs-chimera.info'.format(args.prefix)
    try:
        os.mkdir(dirname2)
    except:
        pass
    set_quiet(args.quiet, args.debug)

    if args.scaled:
        args.scaled = int(args.scaled)

    # load all the databases
    dblist, ksize, scaled = lca_utils.load_databases(args.lca_db, args.scaled)

    print(ksize, scaled)

    assert len(dblist) == 1
    lca_db = dblist[0]

    ###

    fp = open(args.classify_csv, 'rt')
    r = csv.DictReader(fp, fieldnames=['rank', 'name', 'filename', 'md5sum'])

    fp2 = open('{}-dig.csv'.format(args.prefix), 'wt')
    w = csv.writer(fp2)

    n = 0
    m = 0
    for row in r:
        if row['rank'] in ('MISSED', 'species', 'genus', 'family', 'order'):
            continue
        name = row['name']
        md5sum = row['md5sum']
        sig = sourmash.load_one_signature(os.path.join(dirname, md5sum) + '.sig')

        hashvals = defaultdict(int)
        for hashval in sig.minhash.get_mins():
            hashvals[hashval] += 1

        lineage_counts = summarize(hashvals, dblist, args.threshold)

        match_list = set()
        for lineage, count in lineage_counts.items():
            if lineage:
                pair = lineage[-1]
                if pair.rank == FILTER_AT:
                    match_list.add(lineage)

        if len(match_list) >= 2:
            print(name)
            for lineage in match_list:
                print('   ', ";".join(lca_utils.zip_lineage(lineage)))
            print('----\n')

            fp3 = open(os.path.join(dirname2, row['md5sum']) + '.txt', 'wt')
            for lineage, count in lineage_counts.items():
                printme = 0
                if lineage in match_list:
                    fp3.write("{} {}\n".format(count, ";".join(lca_utils.zip_lineage(lineage))))
            fp3.close()
                
            n += 1

            w.writerow(['chimera', row['name'], row['filename'], row['md5sum']])
        else:
            w.writerow(['other', row['name'], row['filename'], row['md5sum']])
            m += 1

    print(n, m)


if __name__ == '__main__':
    main(sys.argv[1:])
