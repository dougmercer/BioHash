#!/usr/bin/python3
import random


def get_worst_word(exp, hash_fcn):
    worst_word = 'z'
    worst_val = float('-inf')
    for word in exp.reads:
        hword, val = hash_fcn(word)
        if val > worst_val:
            worst_val = val
            worst_word = word
    return worst_word, worst_val


def get_best_word(exp, hash_fcn):
    best_word = 'z'
    best_val = float('inf')
    for word in exp.reads:
        hword, val = hash_fcn(word)
        if val < best_val:
            best_val = val
            best_word = word
    return best_word, best_val


def gen_words(num_words, word_length, alphabet):
    words = []
    for _ in range(num_words):
        words.append(rand_word(alphabet, word_length))
    return words


def rand_word(alphabet, word_length):
    return ''.join(random.SystemRandom().choice(alphabet) for _ in range(word_length))


def hash_compare(exp):
    hashes = []
    wrongs = 0
    for read in exp.reads:
        hn = exp.naive_hash(read)
        hs = exp.shift_hash(read)
        hashes.append([hn, hs])
        if hn[0] != hs[0]:
            wrongs += 1
    return hashes, wrongs


def hash_wrong(exp):
    hashes = []
    reads = []
    wrong = 0
    for read in exp.reads:
        hn = exp.naive_hash(read)
        hs = exp.shift_hash(read)
        if hn[0] != hs[0]:
            hashes.append([hn, hs])
            reads.append(read)
            wrong += 1
    return hashes, reads, wrong


def print_lines(lines):
    for line in lines:
        print(line)
