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
lookup_dir = '/home/dmercer/var/bioinf2/'


class Experiment(object):
    """
    An Experiment holds a set of Kmers, as well as a look-up table for the
    shift rule on substrings of length sslen
    """
    def __init__(self, reads, sslen, alphabet):
        self.reads = reads
        self.sslen = sslen
        self.alphabet = alphabet
        self.lookup = make_lookup(sslen, alphabet)

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


def make_lookup(sslen, alphabet):
    lookup = {}
    for length in range(1, sslen + 1):
        lookup_path = lookup_dir + 'table' + str(length) + '.p'
        if os.path.isfile(lookup_path):
            lookup_small = pickle.load(open(lookup_path, 'rb'))
        else:
            lookup_small = _make_lookup_slice(length, alphabet)
            pickle.dump(lookup_small, open(lookup_path, 'wb'))
        lookup.update(lookup_small)
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
        if (string[i:n] + 'z'*(n-i)) < curstr:
            curstr = string[i:n] + 'z'*(n-i)
            value = i
    if value == 0:
        correction = _suffix_prefix_correction(string)
        if correction != -1:
            value = correction
        else:
            value = n
    return value


def _suffix_prefix_correction(string):
    # This function currently induces O(nk) worst case time
    # Must figure out a way to simultaneously resolve the following two strs:
    # cccccccc - correct answer should be shift = 8
    # acagcaca - correct answer should be shift = 5
    n = len(string)
    best_suffix = 'z'
    correction = -1
    for i in range(1, n):
        if string[0:i] == string[n - i:n]:
            padded_suffix = string[0:i] + 'a'*(n-i)
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
    # words = BioTest.gen_words(num_words, word_length, alphabet)
    words = ['aacgcaacgcaacgcaacgcaacgcaacgcaacgcaacgcaacgcaacgcaacgcaacgcaacgcaa']
    exp = Experiment(words, sslen, alphabet)
    [hashes, wrong] = BioTest.hash_compare(exp)
    # [hashes, reads, wrong] = BioTest.hash_wrong(exp)
    # BioTest.print_lines(reads)
    BioTest.print_lines(hashes)
    print(wrong)
    # print(BioTest.get_worst_word(exp, exp.shift_hash))
    # Note: results of hashing should be stored in Experiment
    # print(BioTest.get_best_word(exp, exp.shift_hash))
    # Current configuration:
    # test 1,000,000 random strings of length 100
    # find 8-minimizer
    # result
    # 0 errors
    # worst_word was tgaccggaacaggaaaagtcctttccagcagtacgcgaaaaagtagcaggcggtcttggacgaaggtggcagaacagtgaaacacaaaaaaaggtattag
    # which required 126 comparisons