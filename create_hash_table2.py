__author__ = 'michael'


"""
CREATES A HASH TABLE WITH
"""


import cPickle as pickle
from collections import defaultdict
import time
import os


MAX_SIZE = 10**4
START_STRINGS = [w + x + y + z for w in 'ACTG' for x in 'ACTG' for y in 'ACTG' for z in 'ACTG']
ALL_DICTS = {u: {'dict': defaultdict(list),
                 'size': 0,
                 'dict_count': 0} for u in START_STRINGS}
KEY_LENGTH = 15

def read_chromosome(reference_fn):
    with open(reference_fn, 'r') as reference_file:
        print reference_file.readline()
        chr_index = int(reference_file.readline().strip()[4:])
        return chr_index, ''.join([line.strip() for line in reference_file])


def hash_chromosome(reference_fn, output_folder):
    global ALL_DICTS
    chr_index, seq = read_chromosome(reference_fn)
    all_dicts_key_length = len(START_STRINGS[0])
    for i in range(len(seq) - KEY_LENGTH):
        if i % 10**6 == 0:
            print i
        key = seq[i: i+KEY_LENGTH]
        value = (chr_index, i)
        all_dicts_key = key[:all_dicts_key_length]
        ALL_DICTS[all_dicts_key]['dict'][key].append(value)
        ALL_DICTS[all_dicts_key]['size'] += 1
        if ALL_DICTS[all_dicts_key]['size'] >= MAX_SIZE:
            dict_count = ALL_DICTS[all_dicts_key]['dict_count']
            output_fn = '_'.join([all_dicts_key, str(dict_count).zfill(4)]) + '.pkl'
            output_fn = os.path.join(output_folder, output_fn)
            with open(output_fn, 'w') as output_file:
                pickle.dump(ALL_DICTS[all_dicts_key]['dict'], output_file)
            ALL_DICTS[all_dicts_key] = {'dict': defaultdict(list),
                                        'size': 0,
                                        'dict_count': dict_count + 1}
            print all_dicts_key, dict_count
    return


def hash_genome(genome_folder, output_folder):
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)
    chrom_fns = [os.path.join(genome_folder, x) for x in os.listdir(genome_folder)]
    for chrom_fn in chrom_fns[:1]:
        hash_chromosome(chrom_fn, output_folder)
    for k in ALL_DICTS:
        dict_count = ALL_DICTS[k]['dict_count']
        output_fn = '_'.join([k, str(dict_count).zfill(4)]) + '.pkl'
        output_fn = os.path.join(output_folder, output_fn)
        with open(output_fn, 'w') as output_file:
            pickle.dump(ALL_DICTS[k]['dict'], output_file)


def load_from_hash(hash_folder, start_seq):
    fns = [os.path.join(hash_folder,x) for x in os.listdir(hash_folder) if start_seq in x]
    output_dict = defaultdict(list)
    for fn in fns:
        new_dict = pickle.load(open(fn, 'r'))
        for k in new_dict:
            output_dict[k] += new_dict[k]
    return output_dict


if __name__ == "__main__":
    input_folder = './H_genome/ref'
    output_folder = './H_genome/hash'
    t = time.clock()
    # hash_genome(input_folder, output_folder)
    output_dict = load_from_hash(output_folder, 'AAA')
    print time.clock() - t
    pickle.dump(output_dict, open('test2.pkl','w'))