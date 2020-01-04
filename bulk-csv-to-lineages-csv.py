#! /usr/bin/env python
import os
import csv
import sys
import argparse
import pprint

import sourmash
from sourmash.lca import lca_utils


def main():
    p = argparse.ArgumentParser()
    p.add_argument('bulk_classify_csvs', nargs='+')
    p.add_argument('lineages_csv_out')
    p.add_argument('-v', '--verbose', action='store_true')
    args = p.parse_args()

    ident_to_tax = {}
    for filename in args.bulk_classify_csvs:
        with open(filename, 'rt') as fp:
            r = csv.DictReader(fp)
            n = 0
            for n, row in enumerate(r):
                #ident = row['name'].split(' ')[0].split('.')[0]
                ident = row['name']
                if ident in ident_to_tax:
                    print('*** WARNING *** ident {} occurs more than once!?'.format(ident))
                ident_to_tax[ident] = row['lineage'].split(';')

        print('loaded {} rows from bulk classify csv {}'.format(n + 1, filename))

    print('{} rows total'.format(len(ident_to_tax)))

    with open(args.lineages_csv_out, 'wt') as fp:
        w = csv.writer(fp)
        w.writerow('accession,superkingdom,phylum,class,order,family,genus,species'.split(','))
        for ident, tax in ident_to_tax.items():
            row = [ident] + tax
            w.writerow(row)

    print('done! {} rows written to {}'.format(len(ident_to_tax), args.lineages_csv_out))


if __name__ == '__main__':
    main()
