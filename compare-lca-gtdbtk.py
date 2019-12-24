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
    p.add_argument('gtdbtk_dir')
    p.add_argument('lca_classify_out')
    p.add_argument('-v', '--verbose', action='store_true')
    args = p.parse_args()

    lca_d = {}
    with open(args.lca_classify_out, 'rt') as fp:
        r = csv.DictReader(fp)
        for row in r:
            lca_d[row['ID']] = row
    print('loaded {} rows from sourmash lca classify in {}'.format(len(lca_d), args.lca_classify_out))

    tk_d = {}
    bac_summary = os.path.join(args.gtdbtk_dir, 'gtdbtk.bac120.summary.tsv')
    with open(bac_summary, 'rt') as fp:
        r = csv.DictReader(fp, delimiter='\t')
        for row in r:
            tk_d[row['user_genome']] = row
    ar_summary = os.path.join(args.gtdbtk_dir, 'gtdbtk.ar122.summary.tsv')
    with open(ar_summary, 'rt') as fp:
        r = csv.DictReader(fp, delimiter='\t')
        for row in r:
            tk_d[row['user_genome']] = row

    print('loaded {} rows from gtdbtk classify_wf in dir {}'.format(len(tk_d), args.gtdbtk_dir))

    # make keys match
    lca_d_2 = {}
    for k, v in lca_d.items():
        k = os.path.basename(k)
        if k.endswith('.gz'):
            k = k[:-3]
        lca_d_2[k] = v
    lca_d = lca_d_2

    lca_keys = set(lca_d.keys())
    tk_keys = set(tk_d.keys())

    print('{} input genomes'.format(len(lca_d)))

    if lca_keys - tk_keys:
        print('gtbtk failed to output:')
        pprint.pprint(lca_keys - tk_keys)

    assert not (tk_keys - lca_keys)

    # remove the nomatch from lca classify
    n = 0
    m = 0
    for row in lca_d.values():
        if row['status'] == 'nomatch':
            n += 1
        if not row['species']:
            m += 1

    print('lca classify completely failed on: {}'.format(n))
    print('lca classify failed to assign species on: {}'.format(m))

    lca_lineages = {}
    for k, row in lca_d.items():
        if row['status'] != 'nomatch':
            lineage = []
            for rank in lca_utils.taxlist(include_strain=False):
                name = row[rank]
                lineage.append(lca_utils.LineagePair(rank=rank, name=name))
            lineage = tuple(lineage)
            lca_lineages[k] = lineage

    tk_lineages = {}
    for k, row in tk_d.items():
        classify = row['classification']
        classify = classify.split(';')

        idx = len(classify) - 1
        while classify[idx].endswith('__'):
            assert len(classify[idx]) == 3
            idx -= 1
        classify = classify[:idx + 1]
        lineage = []
        for rank, name in zip(lca_utils.taxlist(include_strain=False), classify):
            lineage.append(lca_utils.LineagePair(rank=rank, name=name))

        tk_lineages[k] = lineage

    both_lin = set(tk_lineages).union(set(lca_lineages))
    n_match1 = 0
    n_match2 = 0
    total = 0
    for k in both_lin:
        tk_lin = tk_lineages.get(k, ())
        lca_lin = lca_lineages.get(k, ())

        if not tk_lin or not lca_lin:
            if not tk_lin and not lca_lin:
                print('skipping comparison b/c absent from both gtdb-tk and sourmash:', k)
            elif not tk_lin:
                print('skipping comparison b/c absent from gtdb-tk:', k)
            elif not lca_lin:
                print('skipping comparison b/c absent from sourmash:', k)
            continue

        total += 1

        tree = lca_utils.build_tree((tk_lin, lca_lin))
        tree_lca, reason = lca_utils.find_lca(tree)

        lca_lin2 = list(lca_lin)
        while lca_lin2 and not lca_lin2[-1].name:
            lca_lin2.pop()

        display = False
        tree_lca = list(tree_lca)
        if tree_lca[:len(lca_lin2)] == lca_lin2:
            n_match1 += 1
        else:
            if args.verbose:
                print('** mismatch at sourmash lca classification level', k)
            display = True

        if tree_lca == list(lca_lin):
            n_match2 += 1
        else:
            display = True
            if args.verbose:
                print('** mismatch at full GTDB-Tk level', k)

        if display and args.verbose:
            print('tk', k, ';'.join(lca_utils.zip_lineage(tk_lin, include_strain=False)))
            print('lca', k,  ';'.join(lca_utils.zip_lineage(lca_lin, include_strain=False)))
            
#        pprint.pprint(tree_lca)
#        pprint.pprint(lca_lin)
#        print('------')

    print('match at full sourmash lca: {:.1f}'.format(n_match1 / total * 100))
    print('match at full GTDB-Tk lineage: {:.1f}'.format(n_match2 / total * 100))




if __name__ == '__main__':
    main()
