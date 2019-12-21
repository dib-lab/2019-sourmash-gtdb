#! /usr/bin/env python
import argparse
import sourmash
import screed


def main():
    p = argparse.ArgumentParser()
    p.add_argument('query_genome')
    p.add_argument('subject_genome')
    p.add_argument('-o', '--output')
    p.add_argument('-k', '--ksize', type=int, default=21)
    p.add_argument('--scaled', type=int, default=500)
    p.add_argument('--threshold', type=float, default=0.8)
    args = p.parse_args()

    minhash = sourmash.MinHash(n=0, ksize=args.ksize, scaled=args.scaled)

    assert args.output
    ofp = open(args.output, 'wt')

    db_minhash = minhash.copy_and_clear()
    for record in screed.open(args.subject_genome):
        db_minhash.add_sequence(record.sequence, False)

    for record in screed.open(args.query_genome):
        q_minhash = minhash.copy_and_clear()
        q_minhash.add_sequence(record.sequence)
        if q_minhash.contained_by(db_minhash) >= args.threshold:
            print('saving', record.name)
            ofp.write('>{}\n{}\n'.format(record.name, record.sequence))


if __name__ == '__main__':
    main()
    
