#! /usr/bin/env python
import argparse
import glob
import csv
import os
from pymummer import coords_file, alignment, nucmer
import gzip


def find_genome_filename(genomes_dir, ident):
    pattern = os.path.join(genomes_dir, ident + '*_genomic.fna.gz')
    matches = glob.glob(pattern)
    if len(matches) != 1:
        assert 0, "more than one match to {}; {}".format(ident, matches)

    return matches[0]


def main():
    p = argparse.ArgumentParser()
    p.add_argument('oddities_csv')
    p.add_argument('genomes_dir', help='fastani database dir')
    p.add_argument('--percent-threshold', type=float,
                   default=95.0)
    p.add_argument('--length-threshold', type=int, default=500)
    p.add_argument('-v', '--verbose', action='store_true')
    args = p.parse_args()

    print('loading', args.oddities_csv)
    print('getting genomes from:', args.genomes_dir)
    print('length threshold for alignments (bp):', args.length_threshold)
    print('lower cutoff for identity (%):', args.percent_threshold)

    fp = open(args.oddities_csv, 'rt')
    r = csv.DictReader(fp)

    try:
        os.mkdir('alignments')
    except:
        pass

    for row in r:
        cluster_name = row['cluster']
        ident1 = row['ident1']
        ident2 = row['ident2']

        fn1 = find_genome_filename(args.genomes_dir, ident1)
        fn2 = find_genome_filename(args.genomes_dir, ident2)

        if args.verbose:
            print(cluster_name, ident1, ident2)

        genome1 = os.path.join('alignments', ident1 + '.fa')
        with gzip.open(fn1) as fp1:
            with open(genome1, 'wb') as fp2:
                      fp2.write(fp1.read())
        genome2 = os.path.join('alignments', ident2 + '.fa')
        with gzip.open(fn2) as fp1:
            with open(genome2, 'wb') as fp2:
                      fp2.write(fp1.read())
            
        nucmer_output_name = os.path.join('alignments', cluster_name + '.a')

        if not os.path.exists(nucmer_output_name):
            print('running...')
            runner = nucmer.Runner(genome1, genome2, nucmer_output_name)
            runner.run()
            print('...done!')
        else:
            if args.verbose:
                print('using cached alignments file', nucmer_output_name)

        file_reader = coords_file.reader(nucmer_output_name)
        alignments = [coord for coord in file_reader if not coord.is_self_hit()]

        # alignment obj:
        # 'frame', 'hit_length_qry', 'hit_length_ref', 'intersects_variant', 'is_self_hit', 'on_same_strand', 'percent_identity', 'qry_coords', 'qry_coords_from_ref_coord', 'qry_end', 'qry_length', 'qry_name', 'qry_start', 'ref_coords', 'ref_coords_from_qry_coord', 'ref_end', 'ref_length', 'ref_name', 'ref_start', 'reverse_query', 'reverse_reference', 'to_msp_crunch']

        #alignments.sort(key = lambda x: -x.hit_length_qry)

        sum_bp = 0
        for alignment in alignments:
            if alignment.hit_length_qry >= args.length_threshold and \
               alignment.percent_identity >= args.percent_threshold:
                sum_bp += alignment.hit_length_qry
        lca_name = "(root)"
        if row['lca']:
            lca_name = row['lca']
        print(cluster_name, sum_bp, row['shared_bp'], lca_name)


if __name__ == '__main__':
    main()
