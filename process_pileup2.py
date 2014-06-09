__author__ = 'michael'

'''
FUCK SO STUPEEED
'''


import os
import sys
from collections import defaultdict
from ast import literal_eval
import string
import swalign
import zipfile
SCORING = swalign.NucleotideScoringMatrix()
ALIGNER = swalign.LocalAlignment(SCORING, globalalign=True, gap_penalty=-5)
UPPERCASE = set(string.ascii_uppercase)
AGREEMENT_THRESHOLD = .99
PRIOR_WEIGHT = 3


def read_reference(ref_fn):
    with open(ref_fn, 'r') as ref_file:
        genome_name = ref_file.readline().strip()[1:]
        chrom_name = ref_file.readline().strip()[4:]
        ref = ''.join([line.strip() for line in ref_file])
    return genome_name, chrom_name, ref


def process_line(line):
    """
    :param line:
    :return:
    """
    start_point, read, cigar = line.strip().split('\t')
    start_point = literal_eval(start_point)[1]
    instructions = cigar_to_instructions(cigar)

    ## We drop insertions at the beginning of a read
    first_instruction = instructions[0]
    first_instruction_type = first_instruction[1]
    if first_instruction_type == 'I':
        first_instruction_run = first_instruction[0]
        read = read[first_instruction_run:]
        instructions = instructions[1:]
    return start_point, read, instructions

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
            number_string = ''
            instruction = letter
            instruction_list.append([number, instruction])
        else:
            number_string += letter
    return instruction_list


def get_current_base(read_start, read_str, read_instructions):
    current_instruction = read_instructions[0]
    current_instruction_type = current_instruction[1]
    current_instruction_run = current_instruction[0]
    if current_instruction_type in 'MX':
        current_base, new_read_str = read_str[0], read_str[1:]
    elif current_instruction_type == 'D':
        current_base, new_read_str = '-', read_str
    else:  # We should NEVER run into an insertion by itself, by construction.
        raise ValueError

    if current_instruction_run == 1:  # If the current instruction run is length 1, we may need to append an insertion.
        if len(read_instructions) == 1:
            new_instructions = read_instructions[1:]
            new_read_str = read_str[1:]
            assert new_read_str == ''
        else:
            next_instruction = read_instructions[1]
            next_instruction_type = next_instruction[1]
            if next_instruction_type == 'I':
                next_instruction_run = next_instruction[0]
                current_base += new_read_str[:next_instruction_run]
                new_read_str = new_read_str[next_instruction_run:]
                new_instructions = read_instructions[2:]
            else:
                new_instructions = read_instructions[1:]
    else:
        new_instructions = read_instructions
        new_instructions[0][0] -= 1

    new_active_line = [read_start + 1, new_read_str, new_instructions]

    return current_base, new_active_line


def generate_consensus_base(bases, reference_base, ref_pos):
    """
    Returns the
    :param bases:
    :param reference_base:
    :return:
    """
    if not bases:
        return reference_base, None, None

    base_dicts = []
    for base in bases:
        for i in range(len(base)):
            if len(base_dicts) == i:
                base_dicts.append(defaultdict(float))
            this_base = base[i]
            base_dicts[i][this_base] += 1

    n_bases = sum(base_dicts[0].values())

    consensus_base = ''
    for d in base_dicts:
        total_reads = 0
        max_read = 0
        for k in d:
            value = d[k]
            total_reads += value
            if value > max_read:
                output_base = k
        if float(total_reads)/n_bases < .5:
            break
        else:
            consensus_base += output_base

    first_base = consensus_base[0]
    if first_base == reference_base:
        match_record = None
    elif first_base == '-':
        match_record = ('D', reference_base, ref_pos)
    else:
        match_record = ('X', reference_base, first_base, ref_pos)
    if len(consensus_base) == 1:
        insertion_record = None
    else:
        insertion_record = ('I', consensus_base[1:], ref_pos)

    return consensus_base, match_record, insertion_record


def generate_consensus_sequence(reads_fn, reference):
    consensus_sequence = []
    consensus_coverage = []
    records = defaultdict(list)
    active_lines = []  # An active line is of the form (start point, read, instructions
    with open(reads_fn, 'r') as pileup_file:
        current_line = pileup_file.readline()
        current_read_start, current_read_str, current_read_instruction = process_line(current_line)
        ref_pos = 0
        while True:
            # print ref_pos
            ref_base = reference[ref_pos]
            if current_read_start == ref_pos:
                active_lines.append([current_read_start, current_read_str, current_read_instruction])
                current_line = pileup_file.readline()
                if current_line:
                    current_read_start, current_read_str, current_read_instruction = process_line(current_line)
                    continue
                else:
                    pass
            current_bases_and_new_active_lines = [get_current_base(*active_line) for active_line in active_lines]
            if current_bases_and_new_active_lines:
                current_bases, active_lines = zip(*current_bases_and_new_active_lines)
                active_lines = [line for line in active_lines if line[1]]
            else:
                current_bases = []
                assert not active_lines

            base_coverage = len(current_bases)
            consensus_current_base, match_record, insert_record = generate_consensus_base(current_bases, ref_base, ref_pos)

            consensus_sequence.append(consensus_current_base)
            consensus_coverage.append(base_coverage)

            if match_record:
                record_type, record_info = match_record[0], match_record[1:]
                records[record_type].append(record_info)
            if insert_record:
                record_type, record_info = insert_record[0], insert_record[1:]
                records[record_type].append(record_info)
            ref_pos += 1
            if ref_pos == len(reference):
                break
    # for k in records:
    #     record_list = records[k]
    #     for record in record_list:
    #         print k, record
    #     print '\n\n-----------\n\n'
    return records


def write_pileup_output(output_dict, genome_name, chrom_name):
    output_fn = 'my_ans_' + genome_name + '.txt'
    print output_fn
    with open(output_fn, 'w') as output_file:
        first_line = '>' + genome_name + '\n'
        output_file.write(first_line)
        second_line = '>chr' + chrom_name + '\n'
        output_file.write(second_line)
        for k in output_dict:
            if k == 'I':
                output_file.write('>INSERT\n')
            elif k == 'D':
                output_file.write('>DELETE\n')
            elif k == 'X':
                output_file.write('>SNP\n')
            records = output_dict[k]
            for record in records:
                output_str = ','.join([chrom_name] + [str(x) for x in record]) + '\n'
                output_file.write(output_str)
    return output_fn


def zip_output(fn):
    zip_fn = fn + '.zip'
    with zipfile.ZipFile(zip_fn, 'w') as output_zip:
        output_zip.write(fn)
    return

if __name__ == "__main__":
    input_folder = './E1_genome'
    pileup_fn = os.path.join(input_folder, 'pileup.txt')
    ref_fn = os.path.join(input_folder, 'ref', 'ref_genomeE1.txt')
    genome_name, chrom_name, ref = read_reference(ref_fn)
    records_dict = generate_consensus_sequence(pileup_fn, ref)
    output_fn = write_pileup_output(records_dict, genome_name, chrom_name)
    zip_output(output_fn)