#!/usr/bin/env python3

import os
import sys
import argparse
import subprocess
from glob import glob
from maptools import filter_markers_blastn
parser= argparse.ArgumentParser(description='''

Tool to Filter Markers from Reference Alignmnet Using BLASTn format 6.

Produces a tsv File with Marker IDs, Scaffold, and Physical position

If a reference is provided, align markers using BLASTn to provided reference

''', formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('--alignments', metavar = '<markers_blast_fmt6.tab>',
required=True,
help='''File of Alignments blast fmt6 or fmt7 std.\n\n''')

parser.add_argument('--reference', metavar = '<reference.fasta>', 
required=True,
help='''Reference FASTA\n\n''')

parser.add_argument('--out_file', metavar = '</path/to/outfile.tab>',
default='./filtered_markers_by_scaffold.tab',
help='''Output file name (default:./filtered_markers_by_scaffold.tab)\n\n''')

parser.add_argument('--min_markers', metavar = '<INT>', default=2,
type=int,
help='''Minimum Markers Aligned to a Scffold to Consider this Scaffold (default:2)\n\n''')

parser.add_argument('--max_targets', metavar = '<INT>', default=5,
type=int,
help='''Maximum Marker Alignment Targets (default:5)

Marker will be assigned to the highest bit-score.  
If bit-scores are tied assigns marker to the longest sequence\n\n''')


parser._optionals.title = '''Program Options'''
args = parser.parse_args()


if __name__ == '__main__':
    alignments = os.path.abspath(args.alignments)
    outfile = os.path.abspath(args.out_file)
    reference = os.path.abspath(args.reference)
    filter_markers_blastn(alignments, reference, outfile, args.min_markers,
                          args.max_targets)

#        blast_markers()
