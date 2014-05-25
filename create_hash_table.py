__author__ = 'michael'

from discodb import DiscoDB
from collections import defaultdict
import os

def create_basic_hash(genome,key_length):
    """
    Builds a basic hash table (a dictionary) using
    uncompressed strings for keys
    """
    genome_hash = defaultdict(list)
    for i in range(len(genome)-key_length):
        key = genome[i:i+key_length]
        genome_hash[key].append(i)
    return genome_hash

def hash_reference(reference_fn,key_length):
    with open(reference_fn,'r') as reference_file:
        first_line = reference_file.readline()
        second_line = reference_file.readline()
        genome = ''.join([line.strip() for line in reference_file])
        return create_basic_hash(genome,key_length)

if __name__ == "__main__":
    input_folder = './EE_genome'
    ref = 'ref_genomeEExample.txt'
    ref_fn = os.path.join(input_folder, ref)
    key_length = 10
    my_hash = hash_reference(ref_fn,key_length)
    k = 'AGATTTAAAC'
    print my_hash[k]
    k = 'CGGATATTAA'
    print my_hash[k]