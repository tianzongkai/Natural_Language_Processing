#! /usr/bin/python
"""
generate parameters q(X -> Y1Y2) and q(X -> w) for pcfg
"""

def create_counts_iterator(count_file):
    """
    Get an iterator object over the corpus file. The elements of the
    iterator contain (word, ne_tag) tuples. Blank lines, indicating
    sentence boundaries return (None, None).
    """
    l = count_file.readline()

    while l:
        line = l.strip().split()
        if line[1] != 'WORDTAG':
            #print line
            yield line
        l = count_file.readline()

def build_rule_dict(counts_iterator):
    rule_dict = {}
    for l in counts_iterator:
        if l[1] == 'UNARYRULE':
            if bigram_key not in bigram_dict:
                bigram_dict[bigram_key] = int(l[0])
        elif l[1] == 'BINARYRULE':
            trigram_key = l[2] + " " + l[3] + " " + l[4]
            if trigram_key not in trigram_dict:
                trigram_dict[trigram_key] = int(l[0])
    return bigram_dict, trigram_dict