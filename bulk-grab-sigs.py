#! /usr/bin/env python
import csv
import argparse
import shutil
import os


def main():
    p = argparse.ArgumentParser()
    p.add_argument('prefixes', nargs='+')
    p.add_argument('--sigdir', help='directory for signatures')
    args = p.parse_args()

    assert args.sigdir, "must supply --sigdir"

    try:
        os.mkdir(args.sigdir)
    except FileExistsError:
        print('warning, sigdir {} already exists'.format(args.sigdir))
        print('continuing...')

    for prefix in args.prefixes:
        csvname = prefix + '-bulk-classify.csv'
        dirname = prefix + '-unclassified-sigs'

        with open(csvname, 'rt') as fp:
            n = 0
            r = csv.DictReader(fp)
            for m, row in enumerate(r):
                if m % 10000 == 0:
                    print(prefix, m, n)
                if row['rank'] in ('superkingdom', 'root'):
                    n += 1
                    sigfile = os.path.join(dirname, row['md5sum'] + '.sig')
                    outfile = os.path.join(args.sigdir, row['md5sum'] + '.sig')
                    shutil.copyfile(sigfile, outfile)

            print(prefix, n)
    

if __name__ == '__main__':
    main()
