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


def cigar_to_instructions(cigar):
    """
    Converts a cigar string to a set of more easily-processable and modifiable set of instructions
    :param cigar:
    :return:
    """
    number_string = ''
    instruction_list = []
    for letter in cigar:
        if letter in UPPERCASE:
            number = int(number_string)
            instruction = letter
            instruction_list.append([number, instruction])
        else:
            number_string += letter
    return instruction_list


def grab_next_base(old_active_lines):
    first_bases = []  # This will be the first base plus any insertions.
    new_active_lines = []

    for line in old_active_lines:
        bases, instructions = line
        current_instruction = instructions[0]
        current_instruction_type = current_instruction[1]
        current_instruction_run = current_instruction[0]
        if current_instruction == 'D':  # If there's a deletion, don't chomp any bases,
            first_base, new_bases = '-', bases
        elif current_instruction == 'I':
            first_base, new_bases = bases[:current_instruction_run], bases[current_instruction_run:]
        elif current_instruction_type in 'MX':
            first_base, new_bases = bases[0], bases[1:]
        else:
            raise ValueError

        if current_instruction == I:

        if len(current_instruction_run) > 1:

            new_instructions = instructions
            new_instructions[0][0] -= 1
        else:


        new_active_lines.append([new_bases, new_instructions])


    return first_bases, new_active_lines


def chomp_cigar(cigar):
    number_str = ''
    for i in range(len(cigar))
        if cigar[i] in UPPERCASE:
            number_str =
        else:
            number_str +=



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
        for i in range(2*len(ref)):
            '''
            This proceeds in two phases. In the first, only matches/mismatches/deletions are
            allowed, and in the second only insertions are allowed.
            '''
            i = float(i)/2.0


            if start_point[1] == i:  # Add to the new set of active reads.
                while True:
                    active_lines.append([read, cigar_to_instructions(cigar)])
                    current_line = pileup_file.readline()
                    if not current_line:
                        break
                    else:
                        start_point, read, cigar = process_line(current_line)
                    if start_point[1] <= i:
                        continue
                    else:
                        break

            next_bases, active_lines = grab_next_base(active_lines)
            if not next_bases:  # If there's no bases, we have to the reference for the sequence.
                output_seq += ref[i]
                continue
            else:  # Examine the reads that are active at the current base,
                   # and determine if there are any SNPs, insertions, or deletions
                next_bases,active = [ for x in ]

    return



if __name__ == "__main__":
    input_folder = './W1_genome'
    pileup_fn = os.path.join(input_folder, 'pileup.txt')
    ref_fn = os.path.join(input_folder, 'ref', 'ref_genomeW1.txt')
    process_pileup(pileup_fn, ref_fn)