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

def to_sequence_array(entry, length, version='OOPS'):

    output = np.empty(shape=(entry.shape[0], length), dtype='str')
    ##NOTE: length needs to reflect the # of columns used; for OOPS: 45; for CUSTOM: 255
    print(entry.shape[0])
    value = None
    if version == 'OOPS':
        for i in range(entry.shape[0]):
            value = entry.iloc[i, :]['switch']+ entry.iloc[i, :]['stem1'] + entry.iloc[i, :]['stem2']
            for j in range(len(value)):
                output[i, j] = value[j]

    if version == 'CUSTOM':
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
def init_guess_nucleotides(seq_array,length):
    output = np.random.rand(4, length)
    for i in range(length):
        total = np.sum(output[:,i])
        output[0,i] = output[0,i] / total
        output[1,i] = output[1,i] / total
        output[2,i] = output[2,i] / total
        output[3,i] = output[3,i] / total
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
    motif_probs = init_guess_nucleotides(seq_letters, k)
    print(motif_probs.shape)

    # Assume nucleotide prob of 0.25 everywhere outside motif, probability of everything outside the motif is 0.25
    #     background = 0.25 ** (motif_probs.shape[1]-k)
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
                    total_motif_probability *= motif_probs[nucleotide_dict[seq_to_train[i + j]], j]

                output_probs.append(total_motif_probability)

            # Weighted sum for Z value at each index, E-Step, using p matrix to estimate Z
            for i in range(seq_to_train.shape[0] - k):
                z_matrix[l, i] = output_probs[i] / (sum(output_probs))

        # Hard assign the max value of Z matrix in each row to 1
        # May need to rethink this particular part of the code
        #         for i in range(z_matrix.shape[0]):
        #             max_idx = np.argmax(z_matrix[i,:])
        #             z_matrix[i,:] = 0
        #             z_matrix[i,max_idx] = 1

        # M step, using z_matrix to re-estimate frequency array
        # Probability of character c in position k along the length of the sequence
        new_freq = np.empty(shape=(z_matrix.shape[0], k), dtype='str')
        for i in range(z_matrix.shape[0]):
            letter_slice = seq_letters[i, np.argmax(z_matrix[i, :]):np.argmax(z_matrix[i, :]) + k]

            new_freq[i, :] = letter_slice
        for i in range(k):
            motif_probs[0, i] = (new_freq[:, i] == 'A').sum() / z_matrix.shape[0]
            motif_probs[1, i] = (new_freq[:, i] == 'G').sum() / z_matrix.shape[0]
            motif_probs[2, i] = (new_freq[:, i] == 'C').sum() / z_matrix.shape[0]
            motif_probs[3, i] = (new_freq[:, i] == 'T').sum() / z_matrix.shape[0]

        iter_track += 1
    print('Done')
    return motif_probs


#def preprocessing(min_high_treshold,max_low_treshold):
def preprocessing(UpperBound, LowerBound):
    dataset = pd.read_csv('data/Toehold_Dataset_Final_2019-10-23.csv')

    #Sort Values based on quality control score
    dataset.sort_values(['QC_ON_OFF'],inplace=True,ascending=True)

    #Filter out scores less than 2 as reported in the study
    dataset = dataset[dataset['QC_ON_OFF'] > 2]
    dataset.head()

    #Now Sort Based on ON_OFF ratio to create the poorly performing and strongly performing switch sets
    dataset.sort_values(['ON_OFF'],inplace=True,ascending=False)
    dataset.head()

    #Drop negative values as their meaning is unclear based on the study
    dataset = dataset[dataset['ON_OFF'] > 0]

    #Save each set to a new csv file
    #low_score = dataset[dataset['ON_OFF'] < max_low_treshold]
    #low_score.to_csv('data/low_score.csv')

    #high_score = dataset[dataset['ON_OFF']> min_high_treshold]
    #temp bracketed high_score

    #high_score = dataset[(LowerBound < dataset['ON_OFF'] < UpperBound)]
    high_score = dataset[(dataset['ON_OFF'] > LowerBound)&(dataset['ON_OFF']) < UpperBound ]
    high_score.to_csv('data/high_score.csv')

if __name__ == '__main__':
    main()