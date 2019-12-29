#! /usr/bin/env python
import sys
import os
import argparse
from sourmash.lca import lca_utils


def main():
    p = argparse.ArgumentParser()
    p.add_argument('lca_db')
    p.add_argument('scrublist')
    p.add_argument('-o', '--output')
    args = p.parse_args()

    (lca_db, ksize, scaled) = lca_utils.load_single_database(args.lca_db)

    scrublist = set()
    with open(args.scrublist, 'rt') as fp:
        for line in fp:
            hashval = int(line)
            scrublist.add(hashval)

    print('loaded {} hashvals from scrublist {}'.format(len(scrublist), args.scrublist))

    for hashval in scrublist:
        del lca_db.hashval_to_idx[hashval]

    filename = 'scrub-{}'.format(args.lca_db)
    if args.output:
        filename = args.output
    print('saving scrubbed LCA db to {}'.format(filename))
    lca_db.save(filename)


if __name__ == '__main__':
    main()
    
