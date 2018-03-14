import os
import numpy as np
import timeit

def read_original_count_file(count_file):
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
    for w in words:
        word_sum = np.sum(counts[words[:] == w])
        #if w == 'Medical' :print 'medical counts', word_sum
        if word_sum < 5 and w not in rare_words:
            rare_words.append(w)
    #rare_words = np.char.add(np.char.add(words[counts[:,0]<5][:,0], ' ') ,words[counts[:, 0] < 5][:, 1])
    #print 'rare_words lenth', len(rare_words)
    return rare_words