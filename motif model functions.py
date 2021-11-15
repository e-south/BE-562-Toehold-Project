'''
--------------------------------------------------------------------------------
<BE562 Project: Toehold Switch Motif Proposal>

Module contains functions to implement the expectation maximization algorithm

Written by Aidan Riley, Juan Montezo, Eric South
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


high_scores = to_sequence_array(high_scores, 155)
low_scores = to_sequence_array(low_scores, 155)


# Function will take a set of sequences as input and return a Matrix of the maximum likelihood estimate or random guess for each motif at each position
# Each row is a nucleotide (A,G,C,T) and each column is the probablity of that nucleotide at that position
def init_guess_nucleotides(seq_array, random=False):
    if random == False:
        output = np.zeros(shape=(4, seq_array.shape[1]))
        total = seq_array.shape[0]
        for i in range(seq_array.shape[1]):
            output[0, i] = sum(seq_array[:, i] == 'A') / total
            output[1, i] = sum(seq_array[:, i] == 'G') / total
            output[2, i] = sum(seq_array[:, i] == 'C') / total
            output[3, i] = sum(seq_array[:, i] == 'T') / total
        return output

    if random == True:
        output = np.random.rand(4, seq_array.shape[1])
        for i in range(output.shape[1]):
            total = np.sum(output[:, i])
            for j in range(output.shape[0]):
                output[j, i] = output[j, i] / total
        return output


# find motif of len k values along a sequence
# Takes high_score or low_score as input along with the length of the motif and the stopping tolerance

def OOPS(seq_letters, k, tol=None):
    # Save the starting position probability array
    z_matrix = np.zeros(shape=(seq_letters.shape[0], seq_letters.shape[1]))

    # Save the nucleotide probability array separate from the sequence array
    # Can be either maximum likelihood or randomized initial guess
    seq_probs = init_guess_nucleotides(seq_letters, random=True)

    #     for k in range(seq)
    # Pop off one sequence to evaluate the probability of
    seq_to_train = seq_letters[0, :]

    # save the output probabilities of the motif
    output_probs = []

    # Assume nucleotide prob of 0.25 everywhere outside motif, probability of everything outside the motif is 0.25
    background = 0.25 ** (seq_probs.shape[1] - k)

    # Iterate through every position in the sequence, up to length of the sequence - k
    for i in range(seq_to_train.shape[0] - k):
        # Calculate the probability of this sequence being a motif, multiply by everything else outside the motif
        motif_probability = background * np.prod(seq_probs[0, i:i + k])
        output_probs.append(motif_probability)

    # Weighted sum for Z value at each index
    for i in range(seq_to_train.shape[0] - k):
        z_matrix[0, i] = output_probs[i] / (sum(output_probs))

    return z_matrix


x = OOPS(high_scores, 10)
sum(x[0, :])


if __name__ == '__main__':
    main()