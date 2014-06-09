
"""
SEE THIS TIME WE'RE GONNA BE SMART SEE

WE MAP THE PERFECTLY MAPPED READS FIRST, THEN, SEE,
WE MAP ALL THE SHIT BETWEEN THE PERFECTLY MAPPED READS

"""
import math
from collections import defaultdict
from scipy.stats import binom
import numpy as np
import os
import zipfile
from create_hash_table import read_chromosome, hash_chromosome


KEY_LENGTH = 10
MAX_INDEL_DIST = 5


def get_min_number_of_matches(n_tests, ref_length, error_cutoff=.01,
                              indel_dist=MAX_INDEL_DIST, key_length=KEY_LENGTH):
    for i in range(1, n_tests):
        probability_of_improper_match = 1.0 - np.exp(ref_length *
                                                     binom.logcdf(i,
                                                                  n_tests,
                                                                  indel_dist*(4.0**(-key_length))))
        if probability_of_improper_match < error_cutoff:
            return i+1
    raise ValueError


def map_paired_end_read(paired_end_read, reference, genome_hash):
    end1, end2 = paired_end_read
    locs1 = map_read(end1, reference, genome_hash)
    locs2 = map_read(end2, reference, genome_hash)
    return locs1, locs2


def map_read(read, ref, genome_hash):
    all_locations = []
    for i in range(len(read) / KEY_LENGTH):
        offset = i * KEY_LENGTH
        key = read[offset: offset + KEY_LENGTH]
        locations = [(start_loc[0],  # The chromosome that the read is on
                      start_loc[1] - offset,  # The "implied start location" of the read
                      start_loc[1],  # The actual start location of the read
                      start_loc[1] + KEY_LENGTH)  # The actual end of the read.
                     for start_loc in genome_hash[key]]
        all_locations.append(locations)
    min_matches = get_min_number_of_matches(len(all_locations), len(ref))
    output = process_locations(all_locations, min_matches, read, ref)
    return output

def process_locations(location_list_list, min_matches, read, reference):
    """
    Returns the set of locations that are
    :param locations_list:
    :return:
    """
    if not any(location_list_list):
        return
    local_alignment_dict = {}
    for location_list in location_list_list:
        used_key_set = set([])
        for location in location_list:
            location_key = location[:2]
            if location_key in used_key_set:
                return None  # DON'T WANNA DEAL WITH THIS CASE
            elif location_key in local_alignment_dict:
                local_alignment_dict[location_key].append(location)
            else:
                chrom_dists = [(k, chromosomal_distance(location_key, k)) for k in local_alignment_dict]
                good_chrom_dists = [u for u in chrom_dists if u[1] < MAX_INDEL_DIST and u[0] not in used_key_set]
                if len(good_chrom_dists) == 0:  # No acceptable nearby location
                    local_alignment_dict[location_key] = [location]
                else:
                    good_key = sorted(good_chrom_dists, key=lambda p: p[0])[0][0]
                    local_alignment_dict[location_key] = local_alignment_dict[good_key] + [location]
                    del local_alignment_dict[good_key]
            used_key_set.add(location_key)
        for k in local_alignment_dict:
            if k in used_key_set:
                pass
            else:
                local_alignment_dict[k].append(None)
    local_alignment_dict_items = local_alignment_dict.items()
    sorted_longest_local_alignment_items = sorted(local_alignment_dict_items, key=lambda q: -len(q[1]))
    longest_local_alignment_item = sorted_longest_local_alignment_items[0]
    longest_local_alignment = longest_local_alignment_item[1]
    if len(longest_local_alignment) < min_matches:
        return None
    report = generate_report(longest_local_alignment, read, reference)
    return report


def generate_report(alignment,  read, reference):
    global PERFECT_MATCHES, MISMATCHES, DELETIONS, INSERTIONS, ENDS
    unmatched_read = ''
    unmatched_ref_start = None
    unmatched_ref_end = None
    for read_index in range(len(alignment)):
        alignment_chunk = alignment[read_index]
        hashed_sequence = read[read_index*KEY_LENGTH: (read_index+1)*KEY_LENGTH]
        if not alignment_chunk:
            unmatched_read += hashed_sequence
            if read_index == len(alignment) - 1:
                ENDS.append((unmatched_ref_start, unmatched_ref_end, unmatched_read))  # FUCK THE ENDS
                unmatched_read = ''

        else:  ## We have a perfect match!
            ref_start_loc, ref_end_loc = alignment_chunk[2:]
            unmatched_ref_end = ref_start_loc
            for j in range(ref_start_loc, ref_end_loc):
                PERFECT_MATCHES[j] += 1

            if unmatched_read:
                if not unmatched_ref_start:
                    ENDS.append((unmatched_ref_start, unmatched_ref_end, unmatched_read))
                else:
                    unmatched_ref_seq = reference[unmatched_ref_start: unmatched_ref_end]
                    ref_dist = len(unmatched_ref_seq)
                    if len(unmatched_read) == ref_dist:
                        # If the length of the unmatched read is equal to the distance in the reference,
                        # we should assume that there is a SNP or a basecall error (as opposed to
                        # offsetting insertions and deletions).
                        called_mismatches = call_mismatches(unmatched_ref_seq, unmatched_read, unmatched_ref_start)
                        MISMATCHES += called_mismatches
                    elif len(unmatched_read) > ref_dist:
                        # If the length of the unmatched read is less than the distance in the reference,
                        # there must be insertions (Or the read is really shittily mapped).
                        # print unmatched_ref_seq, unmatched_read, unmatched_ref_start
                        called_insertions = call_insertion(unmatched_ref_seq, unmatched_read, unmatched_ref_start)
                        INSERTIONS += called_insertions
                        # MISMATCHES += called_mismatches
                    else:
                        # If the length of the unmatched read is more than the distance in the reference,
                        # there must be deletions (Or the read is really shittily mapped).
                        assert len(unmatched_read) < ref_dist
                        called_deletions = call_deletion(unmatched_ref_seq, unmatched_read, unmatched_ref_start)
                        DELETIONS += called_deletions
                        # MISMATCHES += called_mismatches

            # If we hit a perfect match, the new start point of any future mismatched read.
            unmatched_ref_start = ref_end_loc
            unmatched_read = ''

    return alignment


def call_mismatches(ref_seq, read_seq, ref_start_index):
    mismatch_indices = get_mismatch_indices(ref_seq, read_seq)
    if len(mismatch_indices) > 2 + len(ref_seq)/10:
        called_mismatches = []
    else:
        called_mismatches = [(ref_seq[ix], read_seq[ix], ref_start_index + ix) for ix in mismatch_indices]
    return called_mismatches


def call_insertion(ref_seq, read_seq, ref_start_index):
    insertion_amount = len(read_seq) - len(ref_seq)
    if len(ref_seq) == 0:
        called_insertion = (read_seq,ref_start_index)
        return [called_insertion]
    mismatches_after_insertion = [get_mismatch_indices(read_seq[:ix] + read_seq[ix+insertion_amount:],
                                                       ref_seq) for ix in range(len(read_seq) - insertion_amount)]
    mismatch_counts = [len(mismatch_indices) for mismatch_indices in mismatches_after_insertion]

    min_mismatches = min(mismatch_counts)
    if min_mismatches > 1 + len(read_seq)/10:
        called_insertions = []
        # called_mismatches = []
    else:
        optimal_insertion_indices = [ix for ix in range(len(mismatch_counts)) if mismatch_counts[ix] == min_mismatches]
        called_insertions = []
        for optimal_insertion_index in optimal_insertion_indices:
            optimal_inserted_string = read_seq[optimal_insertion_index: optimal_insertion_index + insertion_amount]
            insertion_loc = ref_start_index + optimal_insertion_index
            called_insertion = (optimal_inserted_string, insertion_loc)
            called_insertions.append(called_insertion)  # We include all of the good ones, h
                                                        # oping that the best ones show up more often

            # optimal_remaining_string = ref_seq[:optimal_deletion_index] + \
            #                            ref_seq[optimal_deletion_index + deletion_amount]
            # called_mismatches = call_mismatches(optimal_remaining_string, read_seq, ref_start_index)
    return called_insertions


def call_deletion(ref_seq, read_seq, ref_start_index):
    """
    Attempts to call a single continuous deletion such that
    the ref and read sequence are most closely aligned.
    """
    deletion_amount = len(ref_seq) - len(read_seq)

    mismatches_after_deletion = [get_mismatch_indices(ref_seq[:ix] + ref_seq[ix+deletion_amount:],
                                                      read_seq) for ix in range(len(ref_seq)- deletion_amount)]
    mismatch_counts = [len(mismatch_indices) for mismatch_indices in mismatches_after_deletion]

    min_mismatches = min(mismatch_counts)
    if min_mismatches > 1 + len(ref_seq)/10:
        called_deletions = []
        # called_mismatches = []
    else:
        optimal_deletion_indices = [ix for ix in range(len(mismatch_counts)) if mismatch_counts[ix] == min_mismatches]
        called_deletions = []
        for optimal_deletion_index in optimal_deletion_indices:
            optimal_deleted_string = ref_seq[optimal_deletion_index: optimal_deletion_index + deletion_amount]
            deletion_loc = ref_start_index + optimal_deletion_index
            called_deletion = (optimal_deleted_string, deletion_loc)
            called_deletions.append(called_deletion)
            # optimal_remaining_string = ref_seq[:optimal_deletion_index] + \
            #                            ref_seq[optimal_deletion_index + deletion_amount]
            # called_mismatches = call_mismatches(optimal_remaining_string, read_seq, ref_start_index)
    return called_deletions


def get_mismatch_indices(seq1, seq2):
    assert len(seq1) == len(seq2)
    mismatch_indices = [ix for ix in range(len(seq1)) if seq1[ix] != seq2[ix]]
    return mismatch_indices


def chromosomal_distance(location1, location2):
    chromosome1 = location1[0]
    chromosome2 = location2[0]
    if chromosome1 == chromosome2:
        coord1 = location1[1]
        coord2 = location2[1]
        return abs(coord1 - coord2)
    else:
        return np.inf


def read_and_map_reads(reads_fn, genome_hash, ref_seq, pileup_fn):
    """
    Maps a whole bunch of paired-end reads to a reference genome.
    """
    pileup = []
    with open(reads_fn, 'r') as reads_file:
        reads_file.readline()
        reads_file.readline()
        count = 0
        for line in reads_file:
            paired_end_read = line.strip().split(',')
            # print count
            count += 1
            # if count <= 102:
            #     continue
            mapped1, mapped2 = map_paired_end_read(paired_end_read, ref_seq, genome_hash)
            if mapped1:
                pileup.append(mapped1)
            if mapped2:
                pileup.append(mapped2)
            # if count > 155:
            #     break
    sorted_pileup = sorted(pileup, key=lambda x: x[0]) ## Sort the pileup list by start position
    with open(pileup_fn, 'w') as pileup_file:
        for row in sorted_pileup:
            row_str = '\t'.join([str(x) for x in row]) + '\n'
            pileup_file.write(row_str)
    return


def process_indel_list(indel_list):
    indel_dict = defaultdict(int)
    for indel_record in indel_list:
        if len(indel_record[0]) > 5:
            continue
        indel_dict[indel_record] += 1

    good_indels = remove_conflicting_indels(indel_dict)
    return good_indels


def process_mismatch_list(mismatch_list, perfect_match_list):
    mismatch_dict = defaultdict(float)
    for mismatch_record in mismatch_list:
        mismatch_dict[mismatch_record] += 1

    SNPs = []
    for k in mismatch_dict:
        index = k[2]
        # print k, mismatch_dict[k], perfect_match_list[index]
        if perfect_match_list[index] >= mismatch_dict[k]:
            continue
        else:
            SNPs.append(k)
    SNPs = sorted(SNPs, key=lambda x: x[2])
    return SNPs

def remove_conflicting_indels(indels_dict):
    """
    Remove indels that overlap each other.
    :param indel_dict:
    :return:
    """
    sorted_indels = sorted(indels_dict.items(), key=lambda x: x[0][1])
    good_indels = []
    bad_indels = []
    for i in range(len(sorted_indels)):  # Fuck it, greedy algorithm.
        current_indel, current_indel_count = sorted_indels[i]
        if sorted_indels[i] in bad_indels:
            continue
        conflict_indels = [x for x in sorted_indels[i+1:] if x[0][1] < current_indel[1] + len(current_indel[0])]
        if len(conflict_indels) == 0:
            good_indels.append(current_indel)
        else:
            conflict_indel_counts = [x[1] for x in conflict_indels]
            if current_indel_count >= max(conflict_indel_counts):
                good_indels.append(current_indel)
            else:
                best_indel = conflict_indels[conflict_indel_counts.index(max(conflict_indel_counts))][0]
                good_indels.append(best_indel)
            bad_indels.extend(conflict_indels)
    return good_indels

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
    ref = 'ref/ref_genomeE1.txt'
    reads = 'reads/reads_genomeE1.txt'
    pileup_fn = os.path.join(input_folder, 'pileup.txt')
    ref_fn = os.path.join(input_folder, ref)
    reads_fn = os.path.join(input_folder, reads)
    index, seq = read_chromosome(ref_fn)
    seq_dict = {index: seq}
    my_hash = hash_chromosome(index, seq, KEY_LENGTH)
    PERFECT_MATCHES = [0 for i in range(len(seq))]
    DELETIONS, INSERTIONS, MISMATCHES, ENDS = [], [], [], []
    read_and_map_reads(reads_fn, my_hash, seq, pileup_fn)
    good_deletions = process_indel_list(DELETIONS)
    print '\n\n--------------------------------\n\n'
    for k in good_deletions:
        print k
    print '\n\n--------------------------------\n\n'
    good_insertions = process_indel_list(INSERTIONS)
    for k in good_insertions:
        print k
    print '\n\n--------------------------------\n\n'
    SNPs = process_mismatch_list(MISMATCHES, PERFECT_MATCHES)
    for k in SNPs:
        print k
    output_dict = {'I': good_insertions,
                   'D': good_deletions,
                   'X': SNPs}
    output_fn = write_pileup_output(output_dict, 'genomeE1', str(index))
    zip_output(output_fn)
    # for k in sorted(good_deletions, key=lambda x: x[1]):
    #     print k