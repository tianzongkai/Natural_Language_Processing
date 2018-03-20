#! /usr/bin/python
"""
replace low-frequency words in training file with '_RARE_'
and generate new training file - 'parse_train_rare.dat'

then execute
python count_cfg_freq.py parse_train_rare.dat > cfg_rare.counts
to generate new count file
"""
import sys, os
import numpy as np
import json
import types
import timeit

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

def edit_training_file(train_file, rare_file):
    """
    replace rare words with '_RARE_'
    :return:
    """
    rare_words_list = create_rare_word_list_from_training_file("cfg.counts")
    def modify_leaf(tree):
        """
        recursivly find the leaf level terminal words in a tree
        and replace it with "_RARE_" if the word is a rare word
        :param tree: a parse tree
        :return: the tree with low-frequency leaf words modified
        """
        for idx, item in enumerate(tree):
            if idx != 0:
                if isinstance(item, types.ListType):
                    modify_leaf(item)
                else:
                    if item in rare_words_list:
                        tree[idx] = '_RARE_'
        return tree

    newf = open(rare_file, 'w+')
    with open(train_file) as f:
        for line in f:
            tree = json.loads(line)
            modified_tree = modify_leaf(tree)
            newf.write(str(json.dumps(modified_tree)) + '\n')
    newf.close()
# edit_training_file()


if __name__ == "__main__":
    train_file = sys.argv[1]
    rare_file = sys.argv[2]
    os.system("python count_cfg_freq.py " + train_file + " > cfg.counts")
    edit_training_file(train_file, rare_file)