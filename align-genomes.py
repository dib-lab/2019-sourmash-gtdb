#! /usr/bin/env python
"""
Compare and summarize two genomes by doing nucmer alignments.
"""
import sys
import argparse
import glob
import csv
import os
from pymummer import coords_file, alignment, nucmer
import gzip
from collections import defaultdict
import screed
import shutil
import math


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
    p.add_argument('input_genome1')
    p.add_argument('input_genome2')
    p.add_argument('--percent-threshold', type=float,
                   default=95.0)
    p.add_argument('--length-threshold', type=int, default=0)
    p.add_argument('-v', '--verbose', action='store_true')
    args = p.parse_args()

    ident1 = os.path.basename(args.input_genome1)
    ident2 = os.path.basename(args.input_genome2)
    alignments_dir = 'alignments'

    try:
        os.mkdir(alignments_dir)
    except FileExistsError:
        pass
    
    genome1 = os.path.join(alignments_dir, '{}.fa'.format(ident1))
    xopen = open
    if args.input_genome1.endswith('.gz'): xopen = gzip.open
    with xopen(args.input_genome1, 'rb') as fp1:
        with open(genome1, 'wb') as fp2:
                  fp2.write(fp1.read())

    genome2 = os.path.join(alignments_dir, '{}.fa'.format(ident2))
    xopen = open
    if args.input_genome2.endswith('.gz'): xopen = gzip.open
    with xopen(args.input_genome2, 'rb') as fp1:
        with open(genome2, 'wb') as fp2:
                  fp2.write(fp1.read())

    nucmer_output_name = os.path.join(alignments_dir, ident1 + '.x.' + ident2)

    if not os.path.exists(nucmer_output_name):
        print('running {} alignments...'.format(nucmer_output_name))
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
            print(alignment.hit_length_qry)
            print(alignment.percent_identity)
            skipped_bp += alignment.hit_length_qry
            skipped_aln += 1

    if not keep_alignments:
        print('** FLAG: no kept alignments!')
        sys.exit(-1)

    print('{}.x.{}: {:.0f}kb aln; longest contig: {:.0f} kb'.format(ident1, ident2, aligned_bp / 1000, keep_alignments[0].hit_length_qry / 1000))
    print('weighted percent identity across alignments: {:.1f}%'.format(weighted_percent_identity / all_bp))
    print('skipped {:.0f} kb of alignments in {} alignments (< {} bp or < {:.0f}% identity)'.format(skipped_bp / 1000, skipped_aln, args.length_threshold, args.percent_threshold))

    ###

    keep_d = defaultdict(set)
    for aln in keep_alignments:
        keep_d[aln.qry_name].add(aln)

    bp_removed = remove_contigs(ident2, genome2, keep_d)

    flag_2 = 0
    if bp_removed > 2.5*aligned_bp:
        flag_2 = 1

        # reset to rm kept, and removed is empty.
        os.unlink(genome2 + '.kept.fa')
        with open(genome2 + '.removed.fa', 'wt') as fp:
            pass

    ###

    keep_d = defaultdict(set)
    for aln in keep_alignments:
        keep_d[aln.ref_name].add(aln)

    bp_removed = remove_contigs(ident1, genome1, keep_d)

    flag_1 = 0
    if bp_removed > 2.5*aligned_bp:
        flag_1 = 1
        # reset to rm kept, and removed is empty.
        os.unlink(genome1 + '.kept.fa')
        with open(genome1 + '.removed.fa', 'wt') as fp:
            pass

    if flag_1 and flag_2:
        print('** FLAGFLAG, too much removed from both!')
    elif flag_1 and not flag_2:
        print('** FLAG, {} is probably contaminated (too much rm from {})'.format(ident2, ident1))
    elif flag_2 and not flag_1:
        print('** FLAG, {} is probably contaminated (too much rm from {})'.format(ident1, ident2))

    print('')


if __name__ == '__main__':
    main()
