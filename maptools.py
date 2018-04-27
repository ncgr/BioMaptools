import os
import re
import sys
import subprocess
import gzip
import errno


def markers_tab_to_fasta(markers):
    '''Formats markers.tab into markers.fasta.

       Assumes Column 1 is the Label and Column 2 is the Sequence
    '''
    markers = os.path.abspath(markers)
    with open(markers) as fopen:
        for line in fopen:
            line = line.rstrip()
            if line.startswith('#') or line.startswith('VARIANT') or not line:
                continue  # continue if comment common header or blank
            fields = line.split('\t')  # get fields from line
            if len(fields) != 2:
                print('ERROR: line {}, does nto have 2 fields!'.format(line))
                sys.exit(1)
            label = fields[0]
            sequence = fields[1]
            parts = sequence.split('/')  # get variant parts
            if len(parts) < 2:  # Add log message here skip if not 2 or more
                continue
            new_sequence = []
            count = 0
            part1 = re.compile('[A-z]\[([A-z]+)')
            part2 = re.compile('[A-z]+\]')
            for p in parts:
                if part2.search(p):
                    p = part2.sub('', p)
                if part1.search(p):
                    allele = part1.search(p).groups()[0]
                    length = len(allele) + 1
                    p = '{}{}'.format(p[:-length], allele)
                new_sequence.append(p)
            record = '>{}\n{}'.format(label, ''.join(new_sequence))
            print(record)


def create_directories(dirpath):
    '''make directory path'''
    try:
        os.makedirs(dirpath)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise


def return_filehandle(open_me):
    '''get me a filehandle, common compression or text'''
    magic_dict = {
                  b'\x1f\x8b\x08': 'gz'
#                  '\x42\x5a\x68': 'bz2',
#                  '\x50\x4b\x03\x04': 'zip'
                 }
    max_bytes = max(len(t) for t in magic_dict)
    with open(open_me, 'rb') as f:
        s = f.read(max_bytes)
    for m in magic_dict:
        if s.startswith(m):
            t = magic_dict[m]
            if t == 'gz':
                return gzip.open(open_me, 'rt')
#            elif t == 'bz2':
#                return bz2.open(open_me)
#            elif t == 'zip':
#                return zipfile.open(open_me)
    return open(open_me)


if __name__ == '__main__':
    print('Please import!')
    sys.exit(0)
