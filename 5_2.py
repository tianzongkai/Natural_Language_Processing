#! /usr/bin/python
"""
Impement CKY algorithm
"""
__author__="Zongkai Tian <zt2218@columbia.edu>"
__date__ ="Mar 19, 2018"

import os
import numpy as np
import math
import copy

def create_rare_word_list_from_training_file(count_file):
    """
    read count file and return a list of rare words
    :param count_file: files created by run command line:
            python count cfg freqs.py parse train.dat > cfg.counts
    :return: a list of unique rare words
    """
    lines = []
    with open(count_file) as f:
        for line in f:
            thisline = line.split()
            if thisline[1] == 'UNARYRULE':
                lines.append(thisline[0])  #count
                lines.append(thisline[3])  #word
                lines.append(thisline[2])  # Part-of-speech tag

    word_counts = np.asarray(lines)
    word_counts = np.reshape(word_counts, (-1, 3)) # [count, word, tag]
    counts = word_counts[:,0].astype('int')
    #print counts
    words = word_counts[:, 1]
    rare_words = []
    total_rare = 0
    for w in words:
        word_sum = np.sum(counts[words[:] == w])
        #if w == 'Medical' :print 'medical counts', word_sum
        if word_sum < 5 and w not in rare_words:
            total_rare += word_sum
            rare_words.append(w)
    print 'rare words appear %r times' % total_rare
    return rare_words

def create_counts_iterator(count_file):
    """
    :param count_file: count file
    :return: an iterator
    """
    l = count_file.readline()
    while l:
        line = l.strip().split()
        if line[1] != 'WORDTAG':
            #print line
            yield line
        l = count_file.readline()

def build_rule_count_dict(counts_iterator):
    """
    build a dictionary with counts of each rule
    :param counts_iterator: an iterator of count file
    :return: {X: {w: 12, y1y2:13}}
    """
    rule_count_dict = {}
    for l in counts_iterator:
        if l[1] != 'NONTERMINAL':
            x = l[2]
            y = l[1] == 'UNARYRULE' and l[3] or l[3] + ' ' + l[4]
            # if l[1] == 'UNARYRULE':
            #     y = l[3]
            # else: # l[1] == 'BINARYRULE'
            #     y = l[3] + ' ' + l[4]
            if x not in rule_count_dict:
                rule_count_dict[x] = {}
            rule_count_dict[x][y] = int(l[0])
    return rule_count_dict

def build_para_dict(rule_count_dict):
    para_dict = copy.deepcopy(rule_count_dict)
    for key in para_dict:
        deno = sum(para_dict[key].values())
        for subkey, count in para_dict[key].iteritems():
            para_dict[key][subkey] = float(count)/deno
    return para_dict

def calculate_parameter():
    counts_iterator = create_counts_iterator(file('cfg_rare.counts'))
    rule_count_dic = build_rule_count_dict(counts_iterator)
    para_dict = build_para_dict(rule_count_dic)
    # for key in para_dict:
    #     # print para_dict[key].values()
    #     total = sum(para_dict[key].values())
    #     if total <0.99999999:
    #         print key, total
    return para_dict

def create_sentence_iterator(corpus_file):
    with open(corpus_file) as f:
        l = f.readline()
        while l:
            line = l.split() #turn a string of sentence into a list
            # add a place holder at index 0 to let the sentence starts
            # at index 1, so it can match subscript in the sudo-code
            line.insert(0, '')
            yield line
            l = f.readline()
def cky_algorithm(sentence, parameters):
    n = len(sentence) - 1
    base_case_pi = {}

    for l in range(1, n):
        for i in range(1, n - l + 1):
            j = i + l
            for x in parameters:
                if x != 'S':
                    pass
def init_memo_dict(n, parameters):
    memo_dict = {}
    for i in range(1, n):
        memo_dict[i] = {}
        for j in range(i + 1, n + 1):
            memo_dict[i][j] = {}
            for x in parameters:
                memo_dict[i][j][x] = -1.0
    return memo_dict

def pi(i, j, x, sentence, parameters, memo_dict):
    print 'i, j, x: ', i, j, x
    if i == j:
        w = sentence[i]
        if x in parameters and w in parameters[x]:
            print '11'
            return parameters[x][w]
        else:
            print '22'
            return 0.0
    else:
        if memo_dict[i][j][x] >= 0:
            print '33'
            return memo_dict[i][j][x]
        else:
            print '44'
            rules = parameters[x].keys()
            binary_rules = filter(lambda x: len(x.split()) == 2, rules)
            print binary_rules
            memo_dict[i][j][x] = \
                max([[parameters[x][r]
                      * pi(i, s, r.split()[0], sentence, parameters, memo_dict)
                      * pi(s + 1, j, r.split()[1], sentence, parameters, memo_dict)
                      for s in range(i, j)]
                    for r in binary_rules])
            return memo_dict[i][j][x]


if __name__ == "__main__":
    para_dict = calculate_parameter()
    rare_words_list = create_rare_word_list_from_training_file("cfg.counts")
    s = ['', 'limited', '_RARE_', '.']

    n = len(s) - 1
    print n
    memo_dict = init_memo_dict(n, para_dict)
    # print len(memo_dict[1][3].keys())
    # print memo_dict
    pi = pi(1, n, 'S', s, para_dict, memo_dict)
    print pi
    # sentense_iterator = create_sentence_iterator('parse_dev.dat')
    # for s in sentense_iterator:
    #     print s

