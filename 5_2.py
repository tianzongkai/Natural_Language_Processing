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

import copy

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
    for l in range(1, n):
        for i in range(1, n - l + 1):
            j = i + l

if __name__ == "__main__":
    sentense_iterator = create_sentence_iterator('parse_dev.dat')
    for s in sentense_iterator:
        print s

