'''
--------------------------------------------------------------------------------
<BE562 Project: Toehold Switch Motif Proposal>

Module contains functions to implement the expectation maximization algorithm

Written by Aiden Riley, Juan Montezo, Eric South
--------------------------------------------------------------------------------
'''
# Supporting Libaries
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

#Load in the scoring data
low_scores = pd.read_csv('low_score.csv')
high_scores = pd.read_csv('high_score.csv')


# Functions

# helper function to convert the dataframe into a matrix of useable sequences
# will take pandas low-score or high-score matrix as the input as well as the lenghth of the output sequence

# All functions after here will take a seq_array output as their input

def to_sequence_array(entry, length):
    output = np.empty(shape=(entry.shape[0], length), dtype='str')
    value = None
    for i in range(entry.shape[0]):
        value = entry.iloc[i, :]['pre_seq'] + entry.iloc[i, :]['promoter'] + entry.iloc[i, :]['loop1'] + \
                entry.iloc[i, :]['switch'] + entry.iloc[i, :]['loop2'] + entry.iloc[i, :]['stem1'] + entry.iloc[i, :][
                    'atg'] + entry.iloc[i, :]['stem2'] + entry.iloc[i, :]['linker'] + entry.iloc[i, :]['post_linker']
        for j in range(len(value)):
            output[i, j] = value[j]
    return output


# Function will take a set of sequences as input and return a Matrix of the maximum likelihood estimate for each motif at each position
# Each row is a nucleotide (A,G,C,T) and each column is the probablity of that nucleotide at that position
def maximum_likelihood_estimate(seq_array):
    output = np.zeros(shape=(4, seq_array.shape[1]))
    total = seq_array.shape[0]
    for i in range(seq_array.shape[1]):
        output[0, i] = sum(seq_array[:, i] == 'A') / total
        output[1, i] = sum(seq_array[:, i] == 'G') / total
        output[2, i] = sum(seq_array[:, i] == 'C') / total
        output[3, i] = sum(seq_array[:, i] == 'T') / total
    return output


# find motif of len k values along a sequence

def expectation_maximization(sequence, k):
    pass

if __name__ == '__main__':
    main()