#!/usr/bin/env python3

import sys
import argparse
from glob import glob
from maptools import markers_tab_to_fasta
parser= argparse.ArgumentParser(description='''

Tool to create markers.fasta from markers.tab. 

''', formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('--markers', metavar = '<markers.tab>',
required=True,
help='''Tab delimited file of markers, col1:label, col2:sequence\n\n''')

parser._optionals.title = '''Program Options'''
args = parser.parse_args()


if __name__ == '__main__':
    markers = args.markers
    fasta = markers_tab_to_fasta(markers)
