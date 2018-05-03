#!/usr/bin/env python3

import os
import sys
import argparse
from maptools import format_for_allmaps
parser= argparse.ArgumentParser(description='''

Tool to format marker and scaffold input files for allmaps genetic map input.

''', formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('--markers_gm', metavar = '<markers_gm.tab>',
required=True,
help='''Tab delimited file of markers, col1:label, col2:linkage group, col3:genetic distance\n\n''')

parser.add_argument('--markers_phys', metavar = '<markers_phys.tab>',
required=True,
help='''Tab file of markers with physical positions and scaffold assignments.  col1:marker_label, col2:scaffold_id, col3:physical position (bp)\n\n''')


parser._optionals.title = '''Program Options'''
args = parser.parse_args()


if __name__ == '__main__':
    markers_gm = os.path.abspath(args.markers_gm)
    markers_phys = os.path.abspath(args.markers_phys)
    format_for_allmaps(markers_gm, markers_phys)
