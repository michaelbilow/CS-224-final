__author__ = 'michael'

from discodb import DiscoDB
from collections import defaultdict
import os


# def create_basic_hash(genome, key_length):
#     """
#     Builds a basic hash table (a dictionary) using
#     uncompressed strings for keys
#     """
#     genome_hash = defaultdict(list)
#     for i in range(len(genome) - key_length):
#         key = genome[i:i + key_length]
#         genome_hash[key].append(i)
#     return genome_hash
#
#
# def hash_reference(reference_fn, key_length):
#     with open(reference_fn, 'r') as reference_file:
#         first_line = reference_file.readline()
#         second_line = reference_file.readline()
#         genome = ''.join([line.strip() for line in reference_file])
#         return create_basic_hash(genome, key_length)

def read_chromosome(reference_fn):
    with open(reference_fn, 'r') as reference_file:
        print reference_file.readline()
        chr_index = int(reference_file.readline().strip()[4:])
        return chr_index, ''.join([line.strip() for line in reference_file])


def hash_chromosome(chromosome_id, chromosome_seq, key_length):
    genome_hash = defaultdict(list)
    for i in range(len(chromosome_seq) - key_length):
        key = chromosome_seq[i:i + key_length]
        genome_hash[key].append((chromosome_id, i))
    return genome_hash


if __name__ == "__main__":
    input_folder = './EE_genome'
    ref = 'ref/ref_genomeEExample.txt'
    ref_fn = os.path.join(input_folder, ref)
    key_size = 10
    index, seq = read_chromosome(ref_fn)
    my_hash = hash_chromosome(index, seq, key_size)
    k = 'AGATTTAAAC'
    print my_hash[k]
    k = 'CGGATATTAA'
    print my_hash[k]