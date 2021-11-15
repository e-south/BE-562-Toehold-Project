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
import random

# Load in the scoring data
low_scores = pd.read_csv('data/low_score.csv')
high_scores = pd.read_csv('data/high_score.csv')

nucleotide_dict = {'A': 0, 'G': 1, 'C': 2, 'T': 3}
nucleotide_list = ['A', 'G', 'C', 'T']


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
def init_guess_nucleotides(seq_array, length=None, random=False):
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

def OOPS(seq_letters, k, tol=None, random_init=True, iteration_cap=100):
    print('Training Model...')

    # Set the iteration tracker to 0
    iter_track = 0
    # Save the starting position probability array
    z_matrix = np.zeros(shape=(seq_letters.shape[0], seq_letters.shape[1] - k))

    # Save the nucleotide probability array separate from the sequence array
    # Can be either maximum likelihood or randomized initial guess
    if random_init:
        motif_probs = init_guess_nucleotides(seq_letters, random=True)
    else:
        motif_probs = init_guess_nucleotides(seq_letters, random=False)

    # Assume nucleotide prob of 0.25 everywhere outside motif, probability of everything outside the motif is 0.25
    background = 0.25 ** (motif_probs.shape[1] - k)
    while iter_track < iteration_cap:
        for l in range(seq_letters.shape[0]):

            # Pop off one sequence to evaluate the probability of
            seq_to_train = seq_letters[l, :]

            # save the output probabilities of the motif
            output_probs = []

            # Iterate through every position in the sequence, up to length of the sequence - k
            for i in range(seq_to_train.shape[0] - k):

                total_motif_probability = 1
                for j in range(k):
                    # Calculate the probability of this sequence from position i:k being a motif, multiply by everything else outside the motif

                    total_motif_probability *= motif_probs[nucleotide_dict[seq_to_train[i + k]], i + k]

                output_probs.append(total_motif_probability)

            # Weighted sum for Z value at each index, E-Step, using p matrix to estimate Z
            for i in range(seq_to_train.shape[0] - k):
                z_matrix[l, i] = output_probs[i] / (sum(output_probs))

        # M step, using z_matrix to re-estimate frequency array
        # Probability of character c in position k along the length of the sequence
        for i in range(seq_letters.shape[1] - k):
            for j in range(len(nucleotide_list)):
                character_indices = [seq_letters[:, i] == nucleotide_list[j]]
                numerator = sum(z_matrix[character_indices[0], i])

                denominator = sum(z_matrix[:, i])

                prob = numerator / denominator
                motif_probs[j, i] = prob
        iter_track += 1
    print('Done')
    return motif_probs

if __name__ == '__main__':
    main()