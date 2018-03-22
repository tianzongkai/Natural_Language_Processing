#! /usr/bin/python
"""
Impement CKY algorithm
"""
__author__="Zongkai Tian <zt2218@columbia.edu>"
__date__ ="Mar 19, 2018"

import sys, os
import numpy as np
import math
import copy
import itertools
import re
import time
import json, pprint

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
            if x not in rule_count_dict:
                rule_count_dict[x] = {}
            rule_count_dict[x][y] = int(l[0])
    return rule_count_dict

def build_para_dict(rule_count_dict):
    """
    build a dictionary with probability of each rule
    :param rule_count_dict: a dictionary with count of each rule
    :return: {X: {w: 0.12, y1y2: 0.13}}
    """
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
    return para_dict

def create_frequent_word_list_from_training_file(count_file):
    """
    read original count file and return a list of frequent words
    :param count_file: files created by run command line:
            python count_cfg_freqs.py parse_train.dat > cfg.counts
    :return: a list of unique frequent words
    """
    lines = []
    with open(count_file) as f:
        for line in f:
            thisline = line.split()
            if thisline[1] == 'UNARYRULE':
                lines.append(thisline[0])  # count
                lines.append(thisline[3])  # word
                lines.append(thisline[2])  # Part-of-speech tag

    word_counts = np.asarray(lines)
    word_counts = np.reshape(word_counts, (-1, 3)) # [count, word, tag]
    counts = word_counts[:,0].astype('int')
    #print counts
    words = word_counts[:, 1]
    frequent_words = []
    total_rare = 0
    for w in words:
        word_count = np.sum(counts[words[:] == w])
        if word_count >= 5 and w not in frequent_words:
            frequent_words.append(w)
    # print 'rare words appear %r times' % total_rare
    return frequent_words

def create_sentence_iterator(corpus_file):
    with open(corpus_file) as f:
        l = f.readline()
        while l:
            # turn a string of sentence into a list
            line = l.split()

            # add a place holder at index 0 to let the sentence starts
            # at index 1, so it can match subscripts in algorithm sudo-code
            line.insert(0, '')
            yield line
            l = f.readline()

def init_memo_dict(n, parameters):
    """
    initialize memoization dictionary for dynamic programming,
    all values are initialized to -1
    :param n: length of a sentence
    :param parameters: parameters of a PCFG, which is probability of each rule
    :return: a dictionary in format of {i:{j:{x: -1}}}
    """
    memo_dict = {}
    for i in range(1, n):
        memo_dict[i] = {}
        for j in range(i + 1, n + 1):
            memo_dict[i][j] = {}
            for x in parameters:
                memo_dict[i][j][x] = -1.0
    return memo_dict

def pi(i, j, x, sentence, parameters, memo_pi_dict, memo_bp_dict):
    """
    implementation of PCFG algorithm
    :param i: start index
    :param j: end index
    :param x: non-terminals in left-hand side of rule
    :param sentence: a sentence to parse, in format of list
    :param parameters: q(X -> Y1Y2) and q(x -> w), probability of each rule
    :param memo_pi_dict: lookup dictionary for calculated pi values; {i: {j: {X: pi}}}
    :param memo_bp_dict: backpointer dictionary of argmax s, Y1Y2;
            {i: {j: {X: [s, [Y1, Y2]]}}}
    :return: pi value for the given input
    """
    # print 'i, j, x: ', i, j, x
    # print memo_dict
    if i == j:
        w = sentence[i]
        if x in parameters and w in parameters[x]:
            # print '11: ', parameters[x][w]
            # print w, x
            return parameters[x][w]
        else:
            # print '22'
            return 0.0
    else:
        if memo_pi_dict[i][j][x] >= 0:
            # print '33: ', memo_dict[i][j][x]
            return memo_pi_dict[i][j][x]
        else:
            # print '44'
            # list of righ-hand side of binary & unary rules
            # binary rules are in the format of 'Y1 Y2'
            rules = parameters[x].keys()
            binary_rules = filter(lambda x: len(x.split()) == 2, rules) # Y1, Y2
            # print binary_rules
            if binary_rules == []: # there's no binary rule under input X
                memo_pi_dict[i][j][x] = 0.0
            else:
                sub_pi_matrix = np.asarray(
                    [[parameters[x][r]
                      * pi(i, s, r.split()[0], sentence, parameters, memo_pi_dict, memo_bp_dict)
                      * pi(s + 1, j, r.split()[1], sentence, parameters, memo_pi_dict, memo_bp_dict)
                      for s in range(i, j)]
                     for r in binary_rules]
                )
                memo_pi_dict[i][j][x] = np.amax(sub_pi_matrix)
                argmax_idx = np.argmax(sub_pi_matrix) # index of max value in flatten sub_pi matrix
                argmax_s = i + argmax_idx % (j - i)   # value of s to get max value of pi
                argmax_r = binary_rules[argmax_idx / (j - i)] # get the rule that gives max value of pi
                memo_bp_dict[i][j][x] = [argmax_s, argmax_r]
                # print i, argmax_s, j, argmax_r
            # print 'memo of %r, %r, %r: %r' % (i, j, x, memo_dict[i][j][x])
            return memo_pi_dict[i][j][x]

def build_parse_tree(sentence, i, j, x, bp_dict):
    if i == j:
        return sentence[i]
    bp = bp_dict[i][j][x]
    s = bp[0]
    rule_left = bp[1].split()[0]
    rule_right = bp[1].split()[1]
    return [rule_left,build_parse_tree(sentence, i, s, rule_left, bp_dict)],\
           [rule_right,build_parse_tree(sentence, s + 1, j, rule_right, bp_dict)]

def parse_corpus(dev_file, prediction_file):
    para_dict = calculate_parameter()
    frequent_words = create_frequent_word_list_from_training_file("cfg.counts")
    sentense_iterator = create_sentence_iterator(dev_file)
    newf = open(prediction_file, 'w+')
    total_s = 0
    total_not_s = 0
    for s in sentense_iterator:
        total_s += 1
        n = len(s) - 1
        for i in range(1, n + 1):
            if s[i] not in frequent_words:
                s[i] = '_RARE_'
        memo_pi_dict = init_memo_dict(n, para_dict) # memoization for pi value
        memo_bp_dict = copy.deepcopy(memo_pi_dict)  # memoization for backpointer
        x = 'S'
        prob = pi(1, n, x, s, para_dict, memo_pi_dict, memo_bp_dict)
        if prob == 0.0:
            max_pi = -1.0
            for key in memo_pi_dict[1][n]:
                pi_x = memo_pi_dict[1][n][key]
                if pi_x > max_pi:
                    max_pi = pi_x
                    x = key
        tree = str([x, build_parse_tree(s, 1, n, x, memo_bp_dict)])
        tree = re.sub(r"(?<!\')(\(|\))(?!\')", "", tree)
        tree = re.sub("'", "\"", tree)
        tree = re.sub("\"\"\"\"", "\"''\"", tree)
        tree = re.sub("\"\"", "\"'", tree)
        tree = re.sub("n\"t", "n't", tree)
        newf.write(tree + '\n')
        total_not_s += 1
        # print prob
    newf.close()

def testing():
    para_dict = calculate_parameter()
    # print para_dict
    frequent_words = create_frequent_word_list_from_training_file("cfg.counts")
    s = ['', 'The', 'complicated', 'language', 'in', 'the', 'huge',
         'new', 'law', 'has', 'muddied', 'the', 'fight', '.']
    for i in range(1, len(s)):
        if s[i] not in frequent_words:
            s[i] = '_RARE_'
    print s
    n = len(s) - 1
    prob = -2.0
    memo_pi_dict = init_memo_dict(n, para_dict)
    memo_bp_dict = copy.deepcopy(memo_pi_dict)
    prob = pi(1, n, 'S', s, para_dict, memo_pi_dict, memo_bp_dict)
    # print 'prob,', prob
    tree = ['']
    if memo_pi_dict[1][n]['S'] != 0:
        tree= str(['S', build_parse_tree(s, 1, n, 'S', memo_bp_dict)])
        tree = re.sub(r"(?<!\')(\(|\))(?!\')", "", tree)
        tree = re.sub("'", "\"", tree)
        s = json.loads(tree)
        pretty_print_tree(s)
    else:
        max_pi = -1.0
        x = 'X'
        for key in memo_pi_dict[1][n]:
            pi_value = memo_pi_dict[1][n][key]
            if pi_value > max_pi:
                max_pi = pi_value
                x = key
        tree = [x, build_parse_tree(s, 1, n, x, memo_bp_dict)]
    print 'tree:', tree

class Node:
  """
  Dummy class for python's pretty printer.
  """
  def __init__(self, name): self.name = name
  def __repr__(self): return self.name

def format_tree(tree):
  """
  Convert a tree with strings, to one with nodes.
  """
  tree[0] = Node(tree[0])
  if len(tree) == 2:
    tree[1] = Node(tree[1])
  elif len(tree) == 3:
    format_tree(tree[1])
    format_tree(tree[2])

def pretty_print_tree(tree):
  """
  Print out a tree with nice formatting.
  """
  format_tree(tree)
  print pprint.pformat(tree)

if __name__ == "__main__":
    # start =time.time()

    rare_file = sys.argv[1]
    dev_file = sys.argv[2]
    prediction_file = sys.argv[3]
    # rare_file = "parse_train.RARE.dat"
    # dev_file = "parse_dev.dat"
    # prediction_file ="q5_prediction_file"
    # print rare_file
    # print dev_file
    # print prediction_file

    # testing()
    parse_corpus(dev_file, prediction_file)
    # end = time.time()
    # print "running time %r s" % (end-start)

    # line_num = 0
    # with open(prediction_file) as f:
    #     l = f.readline()
    #     while l:
    #         line_num += 1
    #         print 'line num:', line_num
    #         pretty_print_tree(json.loads(l))
    #         l = f.readline()


