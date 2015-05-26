#!/usr/bin/python3
from itertools import product
import pickle
import os.path
import sys

"""
Dependencies:
python3 (calls to pickle don't seem to work in python2)

Set-up:
Before running, you'll need to create a directory in which to cache your
shift tables. Enter the full path to the directory in the global variable
'lookup_dir', below.

Running:
On the first execution, you'll need to generate all lookup_slices (tables
containing the shift value of all words for a particular word length. This is
RAM intensive. So, start with a modest value of, say, sslen=10 to make sure
everything is working OK.

To run from terminal, type
$ ./BioHash.py sslen num_words word_length
"""

# Magic names/numbers
lookup_dir = '/home/dmercer/var/bioinf3/'


class Experiment(object):
    """
    An Experiment holds a set of Kmers, as well as a look-up table for the
    shift rule on substrings of length sslen
    """
    def __init__(self, reads, sslen, alphabet):
        self.reads = reads
        self.sslen = sslen
        self.alphabet = alphabet
        self.lookup = make_lookup(sslen, alphabet, 'builtin'`)
        self.hashes = []
        self.naive_hashes = []

    def shift_hash(self, read):
        best = 'z'
        str_idx = 0
        iters = 0
        while str_idx <= (len(read) - self.sslen):
            sstring = ''
            for i in range(self.sslen):
                sstring = sstring + read[str_idx + i]
                iters += 1
                if best[i] > read[str_idx + i]:
                    best = read[str_idx:(str_idx + self.sslen)]
                    break
                elif best[i] < read[str_idx + i]:
                    break
            if self.lookup[best] == -1:  # iff best is 'a'*sslen
                return [best, iters]
            shift = self.lookup[sstring]
            # When str_idx > (len(read) - 2*self.sslen), it is possible for
            #   the shift rule to suggest a shift that is too large.
            # To prevent this, when close to the end, we shift conservatively
            # ex: for sslen = 8, 'acacgtctacacaaagcg' would fail otherwise
            if (str_idx + shift + self.sslen) > len(read):
                shift = 1
            str_idx += shift
        return [best, iters]

    def naive_hash(self, read):
        best = 'z'
        iters = 0
        for i in range(len(read) - self.sslen + 1):
            for j in range(0, self.sslen + 1):
                iters += 1
                if read[i+j:i+j+1] < best[j:j+1]:
                    best = read[i:i+self.sslen]
                    if best == 'a'*self.sslen:
                        return [best, iters]
                elif read[i+j:i+j+1] > best[j:j+1]:
                    break
        return [best, iters]

    def hash_all(self, hashfcn):
        hashes = []
        for read in self.reads:
            hashes.append(hashfcn(read))
        return hashes


class RollingDict(dict):
    def __init__(self, sslen, alphabet):
        self.current_hash = 0
        self.itemlist = []
        self.values = self.init_values(sslen, alphabet)

    def __getitem__(self, key):
        [length, index] = self.compute_hash(key)
        return self.values[length, index]

    def __setitem__(self, key, value):
        [length, index] = self.compute_hash(key)
        self.values[length, index] = value
        return self



    def init_values(self, sslen, alphabet):
        values = []
        for length in range(1, sslen + 1):
            values.append([0]*len(alphabet)**length)
        return values

    def update(self, aslice):
        length = len(aslice[0])
        self.values[length] = aslice
        return self


def make_lookup(sslen, alphabet, hashtype):
    if hashtype == 'builtin':
        lookup = {}
        filestr = 'table'
    else:
        lookup = RollingDict()
        filestr = 'rtable'
    for length in range(1, sslen + 1):
        lookup_path = lookup_dir + filestr + str(length) + '.p'
        if os.path.isfile(lookup_path):
            lookup_small = pickle.load(open(lookup_path, 'rb'))
        else:
            lookup_small = _make_lookup_slice(length, alphabet)
            pickle.dump(lookup_small, open(lookup_path, 'wb'))
        lookup.update(lookup_small)
        lookup['a'*length] = length
    lookup['a'*sslen] = -1  # for short circuit if perfect substring
    return lookup


def _make_lookup_slice(length, alphabet):
    lookup_slice = {}
    all_strs = [''.join(x) for x in product(alphabet, repeat=length)]
    for string in all_strs:
        lookup_slice[string] = _compute_lookup_value(string)
    return lookup_slice


def _compute_lookup_value(string):
    n = len(string)
    curstr = 'z'*n
    for i in range(n):
        if (string[i:n] + 'z'*i) < curstr:
            curstr = string[i:n] + 'z'*i
            value = i
    if value == 0:
        correction = _suffix_prefix_correction(string)
        if correction != -1:
            value = correction
        else:
            value = n
    return value


def _suffix_prefix_correction(string):
    n = len(string)
    best_suffix = 'z'
    correction = -1
    for i in range(1, n):
        if string[0:i] == string[n - i:n]:
            padded_suffix = string[0:i] + 'a'*i
            if padded_suffix <= best_suffix:
                best_suffix = padded_suffix
                correction = n - i
    return correction


if __name__ == "__main__":
    import BioTest
    # Read in command-line arguments
    sslen = int(sys.argv[1])
    num_words = int(sys.argv[2])
    word_length = int(sys.argv[3])
    # Default parameters
    alphabet = 'acgt'
    words = BioTest.gen_words(num_words, word_length, alphabet)
    exp = Experiment(words, sslen, alphabet)
    exp.naive_hashes = exp.hash_all(exp.naive_hash)
    exp.hashes = exp.hash_all(exp.shift_hash)
    [hashes, reads, wrong] = BioTest.hash_wrong(exp)
    BioTest.print_lines(hashes)
    BioTest.print_lines(reads)
    print('Number of incorrectly hashed words = ' + str(wrong))
    print('Shift Hash')
    print('Worst:')
    print(BioTest.get_worst_word(exp, hashtype='shift'))
    print('Best:')
    print(BioTest.get_best_word(exp, hashtype='shift'))
    print('Average:')
    print(BioTest.get_avg_comps(exp, hashtype='shift'))
    print('Naive Hash')
    print('Worst:')
    print(BioTest.get_worst_word(exp, hashtype='naive'))
    print('Best:')
    print(BioTest.get_best_word(exp, hashtype='naive'))
    print('Average:')
    print(BioTest.get_avg_comps(exp, hashtype='naive'))
