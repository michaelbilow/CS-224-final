import os
import sys
from collections import defaultdict
from ast import literal_eval
import string

UPPERCASE = set(string.ascii_uppercase)


def read_reference(ref_fn):
    with open(ref_fn, 'r') as ref_file:
        genome_name = ref_file.readline()
        chrom_name = ref_file.readline()
        ref = ''.join([line.strip() for line in ref_file])
    return ref


def process_line(line):
    """
    :param line:
    :return:
    """
    start_point, read, cigar = line.split('\t')
    start_point = literal_eval(start_point)
    return start_point, read, cigar


def cigar_to_read_length(cigar):
    """
    returns the integer read length based on the cigar
    :param cigar:
    :return:
    """
    read_length = 0
    number_string = ''
    for letter in cigar:
        if letter in UPPERCASE:
            number = int(number_string)
            if letter in "MXD":  # If there's a match, mismatch, or deletion, we add this to the total distance
                                 # covered by the read in the reference
                read_length += number
            elif letter == 'I':  # An insertion does not change the distance covered in the reference by the read.
                pass
            else:  # Only allow MXID from extended cigar dictionary.
                raise ValueError
    return read_length


def process_pileup(pileup_fn, ref_fn):
    """
    :param pileup_fn:
    :param ref_fn:
    :return:
    """
    ref = read_reference(ref_fn)
    active_lines = []
    output_seq = ''
    insertions = []
    deletions = []
    SNPs = []
    with(open(pileup_fn, 'r')) as pileup_file:
        current_line = pileup_file.readline()
        start_point, read, cigar = process_line(current_line)
        print start_point, read, cigar
        for i in range(len(ref)):
            if start_point[1] <= i: # Add to the new set of active reads.
                while True:
                    active_lines.append([start_point, read, cigar, 0, cigar_to_read_length(cigar)])
                    current_line = pileup_file.readline()
                    if not current_line:
                        break
                    else:
                        start_point, read, cigar = process_line(current_line)
                    if start_point[1] <= i:
                        continue
                    else:
                        break

            if not active_lines:  # If there's no active lines, we have to the reference for the sequence.
                output_seq += ref[i]
                continue
            else:  # Examine the reads that are active at the current base,
                   # and determine if there are any SNPs, insertions, or deletions
                next_read =

    return



if __name__ == "__main__":
    input_folder = './W1_genome'
    pileup_fn = os.path.join(input_folder, 'pileup.txt')
    ref_fn = os.path.join(input_folder, 'ref', 'ref_genomeW1.txt')
    process_pileup(pileup_fn, ref_fn)