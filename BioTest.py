#!/usr/bin/python3
import random


def gen_words(num_words, word_length, alphabet):
    words = []
    for _ in range(num_words):
        words.append(rand_word(alphabet, word_length))
    return words


def rand_word(alphabet, word_length):
    return ''.join(random.SystemRandom().choice(alphabet) for _ in range(word_length))


def hash_wrong(exp):
    hashes = []
    reads = []
    wrong = 0
    for i in range(0, len(exp.reads)):
        hn = exp.naive_hashes[i]
        hs = exp.hashes[i]
        if hn[0] != hs[0]:
            hashes.append([hn, hs])
            reads.append(exp.reads[i])
            wrong += 1
    return hashes, reads, wrong


def get_worst_word(exp, hashtype='shift', **kwargs):
    worst_val = float('-inf')
    if hashtype == 'shift':
        hashes = exp.hashes
    else:
        hashes = exp.naive_hashes
    for idx, hash_val in enumerate(hashes):
        if hash_val[1] > worst_val:
            worst_word = exp.reads[idx]
            worst_hash = hash_val[0]
            worst_val = hash_val[1]
    return worst_word, worst_hash, worst_val


def get_avg_comps(exp, hashtype='shift', **kwargs):
    comps = 0
    if hashtype == 'shift':
        hashes = exp.hashes
    else:
        hashes = exp.naive_hashes
    for idx, hash_val in enumerate(hashes):
        comps += hash_val[1]
    avg_comps = comps/idx
    return avg_comps


def get_best_word(exp, hashtype='shift', **kwargs):
    best_word = 'z'
    best_val = float('inf')
    if hashtype == 'shift':
        hashes = exp.hashes
    else:
        hashes = exp.naive_hashes
    for idx, hash_val in enumerate(hashes):
        if hash_val[1] < best_val:
            best_word = exp.reads[idx]
            best_hash = hash_val[0]
            best_val = hash_val[1]
    return best_word, best_hash, best_val


def print_lines(lines):
    for line in lines:
        print(line)
