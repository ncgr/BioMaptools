#!/usr/bin/env python3

import os
import sys
import argparse
from glob import glob
from maptools import markers_blastn
parser= argparse.ArgumentParser(description='''

Align marker footprints using BLASTn to provided reference

''', formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('--markers', metavar = '<markers.fasta>',
required=True,
help='''Tab delimited file of markers, col1:label, col2:sequence\n\n''')

parser.add_argument('--reference_db', metavar = '<reference.fasta>',
required=True,
help='''Reference FASTA\n\n''')

parser.add_argument('--blast_path', metavar = '</path/to/ncbi_blast+/bin>',
required=True,
help='''Path to blastn executable\n\n''')

parser._optionals.title = '''Program Options'''
args = parser.parse_args()


if __name__ == '__main__':
    reference = os.path.abspath(args.reference_db)
    blast_path = os.path.abspath(args.blast_path)
    markers = os.path.abspath(args.markers)
    nin = glob('{}*nin'.format(reference))
    if not nin:
        print('ERROR: could not find blast nin for {}'.format(reference))
        sys.exit(1)
    markers_blastn(markers, reference, blast_path)
