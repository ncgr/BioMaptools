import os
import re
import sys
import subprocess
import gzip
import errno
from Bio import SeqIO


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


def get_seqio_record(seq_handle):
    '''Parses a filehandle seq_handle and yields the formatted records

       Generator for SeqIO record objects
    '''
    with seq_handle as sopen:
        for record in SeqIO.parse(sopen, 'fasta'):  # iterate with SeqIO
            yield record  # yield each record as it is iterated


def get_genetic_markers(markers):
    '''Accepts a tab file of genetic markers, their LGs and their genetic

       distance.

       Returns a dictionary with marker id as key and a dict as value.
    '''
    genetic_markers = {}
    fh = return_filehandle(markers)  # get filehandle for markers file
    with fh as mopen:
        for line in mopen:
            line = line.rstrip()
            if not line or line.startswith('#') or line.startswith('VARIANT'):
                continue  # skip comments and headers
            fields = line.split('\t')
            if len(fields) < 3:  # skip empty or unassigned
                continue
            genetic_markers[fields[0]] = {'lg': fields[1], 
                                          'distance': fields[2]}
    return genetic_markers


def format_for_allmaps(markers_gm, markers_phys):
    '''Accepts a tab file of genetic marker labels, their LG and their

       genetic distance.

       Loads this into genetic_markers dict.

       Also accepts a physical position file parsed from filter_marker_blastn

       these will be unioned and an AllMaps CSV file will be produced.
    '''
    markers_gm = os.path.abspath(markers_gm)
    markers_phys = os.path.abspath(markers_phys)
    genetic_markers = get_genetic_markers(markers_gm)
    fh = return_filehandle(markers_phys)
    with fh as popen:
        for line in popen:
            line = line.rstrip()
            if not line or line.startswith('#') or line.startswith('VARIANT'):
                continue  # skip comments and headers
            fields = line.split('\t')
            if len(fields) < 3:  # skip empty or unassigned
                continue
            if not genetic_markers.get(fields[0]):  # skip if no distance etc
                continue
            marker_id = fields[0]
            ref = fields[1]
            position = fields[2]
            lg = genetic_markers[marker_id]['lg']
            genetic_distance = genetic_markers[marker_id]['distance']
            output = '{},{},{},{}'.format(ref, position, lg, genetic_distance)
            print(output)  # print output line formatted above for allmaps


def get_sequence_lengths(fasta):
    '''Gets sequence lengths and returns dictionary 
    
       with ID key and Length value
    '''
    sequence_lengths = {}  # dictionary to return
    fh = return_filehandle(fasta)  # get a filehandle
    for record in get_seqio_record(fh):  # generator for SeqIO
        length = len(record.seq)
        seq_id = record.id
        sequence_lengths[seq_id] = length
    return sequence_lengths


def populate_marker_alignments(alignments, ref):
    '''Populates a dictionary with marker alignment targets,

       their bitscores, physical positions and reference sequence lengths
    '''
    ref_lengths = get_sequence_lengths(ref)  # dict to hold ref lengths
    fh = return_filehandle(alignments)
    marker_alignments = {} 
    with fh as bopen:  # iterate through alignments file
        for line in bopen:
            line = line.rstrip()
            if line.startswith('#') or not line:  # skip if comment or None
                continue
            fields = line.split('\t')
            query = fields[0]
            ref = fields[1]
            ref_length = int(ref_lengths[ref])
            ref_start = int(fields[8])
            bitscore = float(fields[-1])
            last = 0
            if not marker_alignments.get(query):
                marker_alignments[query] = {'alignments': {}, 'num_hits': 0}
            if not marker_alignments[query]['alignments'].get(ref):
                marker_alignments[query]['alignments'][ref] = []
            target = {'length': ref_length, 'start': ref_start,
                      'bitscore': bitscore, 'ref': ref}  # target object
            marker_alignments[query]['alignments'][ref].append(target)
            last = bitscore
            marker_alignments[query]['num_hits'] += 1
    return marker_alignments


def filter_markers_blastn(alignments, ref, outfile, min_markers, max_targets):
    '''Filters BLASTn format 6 output based on provided min_markers and

       max_targets.  min_markers filters scaffold with >= min_markers.

       max_targets filters out markers with >= max_target alignments.
    '''
    alignments = os.path.abspath(alignments)
    ref = os.path.abspath(ref)
    outfile = os.path.abspath(outfile)
    marker_alignments = populate_marker_alignments(alignments, ref)  # markers
    passing_markers = {}
    passing_refs = {}
    for query in marker_alignments:
        if marker_alignments[query]['num_hits'] > max_targets:  # >targets
            continue
        passing_markers[query] = {'length': 0, 'start': 0,
                                  'bitscore': 0, 'ref': ''}
        for ref in marker_alignments[query]['alignments']:
            target = {}
            skip_switch = 0
            bitscore = 0
            if len(marker_alignments[query]['alignments'][ref]) > 1:
                count = 0
                hold = {}
                for align in marker_alignments[query]['alignments'][ref]:
                    if align['bitscore'] > bitscore:
                        hold = align
                        bitscore = align['bitscore']
                    elif align['bitscore'] == bitscore:
                        skip_switch = 1  # cannot use multi target same ref
                        break  # break loop and continue reject this ref
                if not skip_switch:
                    target = hold  # assign marker this target
                if skip_switch:
                    continue  # don't use markers ambiguous on same scaffold
            else:
                target = marker_alignments[query]['alignments'][ref][0]
            assert target, 'Target must be defined'
            bitscore = target['bitscore']
            if bitscore > passing_markers[query]['bitscore']:
                passing_markers[query] = target  # set new target
            elif bitscore == passing_markers[query]['bitscore']:
                if target['length'] > passing_markers[query]['length']:
                    passing_markers[query] = target  # set to longest ref
        if not passing_refs.get(passing_markers[query]['ref']):
            passing_refs[passing_markers[query]['ref']] = 0
        passing_refs[passing_markers[query]['ref']] += 1
    for query in passing_markers:  # now iterate through passing and get refs
        ref = passing_markers[query]['ref']
        if not ref:  # if marker passed but ref was too ambiguous
            continue
        if passing_refs[ref] >= min_markers:  # check min markers
            position = passing_markers[query]['start']
            output_me = '{}\t{}\t{}'.format(query, ref, position)
            print(output_me)


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


if __name__ == '__main__':
    print('Please import!')
    sys.exit(0)
