# ''
# --------------------------------------------------------------------------------
# <BE562 Project: Toehold Switch Motif Proposal>
#
# Module will provide functions for predicting secondary structure binding sites in a single stranded mRNA sequence
#
# Written by Aidan Riley, Juan Montezo, Eric South
# --------------------------------------------------------------------------------
# '''

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import random
import scipy.stats

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


# Function will take a set of sequences as input and return a Matrix of the maximum likelihood estimate or random guess for each motif at each position
# Each row is a nucleotide (A,G,C,T) and each column is the probablity of that nucleotide at that position
def init_guess_nucleotides(seq_array, length, rand=False):
    idx = 0
    if rand == False:
        output = np.zeros(shape=(4, length))
        total = seq_array.shape[0]
        for i in range(length):
            output[0, i] = sum(seq_array[:, idx + i] == 'A') / total
            output[1, i] = sum(seq_array[:, idx + i] == 'G') / total
            output[2, i] = sum(seq_array[:, idx + i] == 'C') / total
            output[3, i] = sum(seq_array[:, idx + i] == 'T') / total
        return output

    if rand == True:
        output = np.zeros(shape=(4, length))
        freqs = np.random.rand(4, length)
        for i in range(length):
            s = sum(freqs[:, i])
            freqs[:, i] = freqs[:, i] / s
            output[:, i] = freqs


high_scores = to_sequence_array(high_scores, 155)
low_scores = to_sequence_array(low_scores, 155)


# Helper functions to score the alignment of two sequences

# Returns compliment of a nucleotide
def compliment(nucleotide):
    if nucleotide == 'A':
        return 'T'
    if nucleotide == 'G':
        return 'C'
    if nucleotide == 'C':
        return 'G'
    if nucleotide == 'T':
        return 'A'


# Takes the reverse of sequence 2
def score(seq_1, seq_2):
    output = 0
    for i in range(len(seq_1)):
        if seq_1[i] == 'A' and seq_2[-(i + 1)] == 'T':
            output = output + 1
        elif seq_1[i] == 'G' and seq_2[-(i + 1)] == 'C':
            output += 1
        elif seq_1[i] == 'C' and seq_2[-(i + 1)] == 'G':
            output += 1
        elif seq_1[i] == 'T' and seq_2[-(i + 1)] == 'A':
            output += 1
        else:
            output = output - 2
    if output < 0:
        return 0
    else:
        return output


# Generate random set of sequence data with secondary strucutre binding determined by a poisson process

# I think the bug is coming from the last distance being the only one returned

def generate_SSS_data(n_samples, sequence_length, motif_length=4, mean_distance=15):
    output = np.empty(shape=(n_samples, sequence_length), dtype='str')
    # Populate the array with random letters
    z_positions = np.zeros(shape=(n_samples, 2))
    for i in range(n_samples):
        for j in range(sequence_length):
            value = np.random.uniform()
            if value <= 0.25:
                output[i, j] = 'A'
            if 0.25 < value <= 0.50:
                output[i, j] = 'G'
            if 0.5 < value <= 0.75:
                output[i, j] = 'C'
            if value > 0.75:
                output[i, j] = 'T'

        # Now assign the motif index positions randomly
        # Random initial index, poisson distribution for the distance between the 2

        distance = np.random.poisson(lam=mean_distance)

        z_1 = np.random.randint(0, (sequence_length - ((2 * motif_length) + distance)))

        z_2 = (distance + motif_length + z_1)

        # Make the sequence at position z_2 the reverse compliment of the sequence at z_1
        for k in range(motif_length):
            output[i, z_2 + ((motif_length) - k - 1)] = compliment(output[i, z_1 + k])

        z_positions[i, 0] = z_1
        z_positions[i, 1] = z_2
    return output, z_positions


# Secondary Structure Search using Expectation Maximization, and assuming poission distribution

def SSS(sequence_array, motif_length):
    z_matrix = np.zeros(shape=(sequence_array.shape[0], sequence_array.shape[1] - (2 * motif_length)))
    positions_matrix = np.zeros(shape=(sequence_array.shape[0], sequence_array.shape[1] - (2 * motif_length)))
    poisson_mean = 1

    # This is where the outer loop will go
    for x in range(sequence_array.shape[0]):
        seq_to_train = sequence_array[x, :]

        # E Step
        for i in range(sequence_array.shape[1] - (2 * motif_length)):
            z_1 = seq_to_train[i:i + motif_length]

            # First value is the index, second is the probability associated with that index
            max_prob = 0

            for j in range(i + motif_length, sequence_array.shape[1] - (2 * motif_length)):

                z_2 = seq_to_train[j:j + motif_length]

                probability = score(z_1, z_2)

                # This check made it run much faster by not appending the 0 entries
                if probability != 0:

                    weighted_probability = probability * scipy.stats.poisson.pmf((j - i), mu=poisson_mean)
                    if weighted_probability > max_prob:
                        z_matrix[x, i] = weighted_probability
                        positions_matrix[x, i] = j
                        max_prob = weighted_probability

    # M Step- Calculate a new poisson mean
    value_save = []
    for x in range(sequence_array.shape[0]):
        max_idx = np.argmax(z_matrix[x, :])

        distance = positions_matrix[x, max_idx] - max_idx

        value_save.append(distance)

    poisson_mean = np.mean(value_save)

    return poisson_mean + motif_length


a = SSS(test_data[0], 8)

print(a)