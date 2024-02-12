#!/usr/bin/env python


import sys
from Bio import SeqIO


if __name__ == '__main__':
    if len(sys.argv) < 2:
        sys.stderr.write('USAGE: python extract_pks.py INPUT.gbk\n')
        sys.exit(1)

    s = SeqIO.read(sys.argv[1], 'genbank')
    for f in s[2566: 52787].features:
        if f.type != 'CDS':
            continue
        g = f.qualifiers.get('gene', [None,])[0]
        p = f.qualifiers['protein_id'][0]
        t = f.qualifiers['translation'][0]
        if g is None:
            g = p
        print(f'>{g}\n{t}')
