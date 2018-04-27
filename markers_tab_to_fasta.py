#!/usr/bin/env python3

import sys
import argparse
import subprocess
from glob import glob
from maptools import markers_tab_to_fasta
parser= argparse.ArgumentParser(description='''

Tool to create markers.fasta from markers.tab. 

If a reference is provided, align markers using BLASTn to provided reference

''', formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('--markers', metavar = '<markers.tab>',
required=True,
help='''Tab delimited file of markers, col1:label, col2:sequence\n\n''')

parser.add_argument('--reference_db', metavar = '<reference.fasta>',
help='''Reference FASTA\n\n''')

parser.add_argument('--blast_path', metavar = '</path/to/ncbi_blast+/bin>',
help='''Path to blastn executable\n\n''')

parser._optionals.title = '''Program Options'''
args = parser.parse_args()


if __name__ == '__main__':
    reference = args.reference_db
    blast_path = args.blast_path
    markers = args.markers
    fasta = markers_tab_to_fasta(markers)
    if reference:
        reference = os.path.abspath(reference)
        nin = glob('{}.*nin')
        if not nin:
            print('ERROR: could not find blast nin for {}'.format(reference))
            sys.exit(1)
        if blast_path:
            blast_path = os.path.abspath(blast_path)
#        blast_markers()
