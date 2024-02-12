#!/usr/bin/env python


import sys
from Bio import SeqIO
from collections import namedtuple


Feature = namedtuple('Feature', ['id',
                                 'chromosome',
                                 'start',
                                 'end',
                                 'strand'])


def parse_gff(file_name, feature_types=None):
    if feature_types is None:
        feature_types = {'CDS'}

    # output dict
    # key: feature ID
    # value: Feature NamedTuple
    features = {}

    with open(file_name, 'r') as gff:
        for line in gff:
            if line.lstrip().startswith('##FASTA'):
                # start of FASTA entries, end of file
                break

            elif line.lstrip().startswith('#'):
                # comment, ignore
                continue

            # should be a valid GFF3 line
            entries = line.split('\t')

            try:
                ftype = entries[2]

                if ftype not in feature_types:
                    continue

                chrom = entries[0]
                start = int(entries[3])
                end = int(entries[4])
                strand = entries[6]

                # integer takes up less space
                if strand == '+':
                    strand = 1
                else:
                    strand = -1

                # fetch the feature ID from the last field
                ID = None
                for entry in entries[8].split(';'):
                    if entry.startswith('ID') and '=' in entry:
                        ID = entry.split('=')[1]

                # could not find it, skip this entry
                if ID is None:
                    continue

                # save the relevant details
                features[ID] = Feature(ID, chrom, start, end, strand)

            except Exception as e:
                # not distinguishing between exceptions
                # not great behaviour
                logger.warning(f'{e}, skipping line "{line.rstrip()}" from {file_name}')
                continue

    return features


if __name__ == "__main__":
    if len(sys.argv) < 3:
        sys.stderr.write('USAGE: python3 gff2faa.py INPUT.gff INPUT.fasta\n')
        sys.exit(1)

    gff, fasta = sys.argv[1: 3]

    d = {s.id: s
         for s in SeqIO.parse(fasta, 'fasta')}

    feats = parse_gff(gff)

    for fid, f in feats.items():
        s = d[f.chromosome]
        if f.strand > 0:
            p = s[f.start-1: f.end].translate().seq
        else:
            p = s[f.start-1: f.end].reverse_complement().translate().seq
        p = str(p)

        if not p.endswith('*') or '*' in p[:-1]:
            sys.stderr.write(f'Skipping {fid} for stop codon irregularities\n')
            continue

        p = p[:-1]

        print(f'>{fid} (f.chromosome) {f.start} {f.end} {f.strand}\n{p}')
