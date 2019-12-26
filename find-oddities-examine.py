#! /usr/bin/env python
"""
Dig into shared nt segments by doing nucmer alignments.
"""
import argparse
import glob
import csv
import os
from pymummer import coords_file, alignment, nucmer
import gzip
from collections import defaultdict
import screed


def find_genome_filename(genomes_dir, ident):
    pattern = os.path.join(genomes_dir, ident + '*_genomic.fna.gz')
    matches = glob.glob(pattern)
    if len(matches) != 1:
        assert 0, "more than one match to {}; {}".format(ident, matches)

    return matches[0]


def remove_contigs(ident, genomefile, keep_d, verbose=True):
    """
    remove contigs from 'genomefile' whose names are in keep_d keys.

    removed contigs go in genomefile + '.removed.fa'
    retained contigs go in genomefile + '.kept.fa'

    return the total number of bp removed.
    """
    bp_skipped = 0
    contigs_skipped = 0
    bp_total = 0
    contigs_total = 0

    kept_outfp = open(genomefile + '.kept.fa', 'wt')
    removed_outfp = open(genomefile + '.removed.fa', 'wt')
    
    for record in screed.open(genomefile):
        bp_total += len(record.sequence)
        contigs_total += 1

        name = record.name.split(' ')[0]
        if name in keep_d:
            # filter
            bp_skipped += len(record.sequence)
            contigs_skipped += 1
            removed_outfp.write('>{}\n{}\n'.format(record.name, record.sequence))
        else:
            kept_outfp.write('>{}\n{}\n'.format(record.name, record.sequence))

    if verbose:
        print('{}: removed {:.0f}kb of {:.0f}kb ({:.0f}%), {} of {} contigs'.format(ident, bp_skipped / 1000, bp_total / 1000, bp_skipped / bp_total * 100, contigs_skipped, contigs_total))

    return bp_skipped


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

    prefix = os.path.basename(args.oddities_csv)
    if prefix.endswith('.csv'):
        prefix = prefix[:-4]

    fp = open(args.oddities_csv, 'rt')
    r = csv.DictReader(fp)

    alignments_dir = prefix + '.alignments'
    print('putting alignments in:', alignments_dir)
    try:
        os.mkdir(alignments_dir)
    except FileExistsError:
        pass

    print('----')

    for row in r:
        cluster_name = row['cluster']
        ident1 = row['ident1']
        ident2 = row['ident2']

        fn1 = find_genome_filename(args.genomes_dir, ident1)
        fn2 = find_genome_filename(args.genomes_dir, ident2)

        if args.verbose:
            print(cluster_name, ident1, ident2)

        genome1 = os.path.join(alignments_dir, ident1 + '.fa')
        with gzip.open(fn1) as fp1:
            with open(genome1, 'wb') as fp2:
                      fp2.write(fp1.read())
        genome2 = os.path.join(alignments_dir, ident2 + '.fa')
        with gzip.open(fn2) as fp1:
            with open(genome2, 'wb') as fp2:
                      fp2.write(fp1.read())
            
        nucmer_output_name = os.path.join(alignments_dir, cluster_name + '.a')

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

        alignments.sort(key = lambda x: -x.hit_length_qry)
        keep_alignments = []

        aligned_bp = 0

        all_bp = 0
        weighted_percent_identity = 0.
        skipped_bp = 0
        skipped_aln = 0

        for alignment in alignments:
            weighted_percent_identity += alignment.percent_identity * alignment.hit_length_qry
            all_bp += alignment.hit_length_qry

            if alignment.hit_length_qry >= args.length_threshold and \
               alignment.percent_identity >= args.percent_threshold:
                aligned_bp += alignment.hit_length_qry
                keep_alignments.append(alignment)
            else:
                skipped_bp += alignment.hit_length_qry
                skipped_aln += 1

        lca_name = "(root)"
        if row['lca']:
            lca_name = row['lca']

        print('{}: {:.0f}kb aln ({:.0f}k {}-mers) across {}; longest contig: {:.0f} kb'.format(cluster_name, aligned_bp / 1000, int(row['shared_kmers']) / 1000, int(row['ksize']), lca_name, alignments[0].hit_length_qry / 1000))
        print('weighted percent identity across alignments: {:.1f}%'.format(weighted_percent_identity / all_bp))
        print('skipped {:.0f} kb of alignments in {} alignments (too short or low identity)'.format(skipped_bp / 1000, skipped_aln))

        ###

        keep_d = defaultdict(set)
        for aln in keep_alignments:
            keep_d[aln.qry_name].add(aln)

        bp_removed = remove_contigs(ident2, genome2, keep_d)

        if bp_removed > 1.5*aligned_bp:
            print('** FLAG, too much removed from', ident2)
            os.unlink(genome2 + '.removed.fa')
            os.unlink(genome2 + '.kept.fa')

        ###

        keep_d = defaultdict(set)
        for aln in keep_alignments:
            keep_d[aln.ref_name].add(aln)

        bp_removed = remove_contigs(ident1, genome1, keep_d)

        if bp_removed > 1.5*aligned_bp:
            print('** FLAG, too much removed from', ident1)
            os.unlink(genome1 + '.removed.fa')
            os.unlink(genome1 + '.kept.fa')

        print('')


if __name__ == '__main__':
    main()
