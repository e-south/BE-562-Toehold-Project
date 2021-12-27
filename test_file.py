import pandas as pd
import numpy as np
import csv

from motif_model_functions import OOPS, to_sequence_array, init_guess_nucleotides, preprocessing

low_scores = pd.read_csv('data/low_score.csv')
high_scores = pd.read_csv('data/high_score.csv')

high_scores = to_sequence_array(high_scores, 45)
low_scores = to_sequence_array(low_scores, 45)

## Changing variables for analysis
kmer_length = [44]
## initial high treshold: 0.85; initial low treshold: 0.05
#min_high_treshold = np.linspace(0.05,2,40)
LowerBound = np.linspace(0.10,0.90,9)
UpperBound = np.linspace(0.15,0.95,9)

##running preprocessing on various tresholds
for score_index in range(len(UpperBound)):

    input1 = UpperBound[score_index]
    input2 = LowerBound[score_index]
    preprocessing(input1,input2)

    ## reading the scores
#    low_scores = pd.read_csv('data/low_score.csv')
    high_scores = pd.read_csv('data/high_score.csv')
    high_scores = to_sequence_array(high_scores, 45)
#    low_scores = to_sequence_array(low_scores, 45)

    for k in kmer_length:
        # Evaluating OOPS for the k in k_mer_length
        x = OOPS(high_scores,k)
        # Writing to a csv file
        # field names
        fields = list(range(k))
        rows = x
        filename = "High_Scores_treshold_X_kmer_N.csv"
        filename = filename.replace('N', str(k))
        filename = filename.replace('X', str((input1+input2)/2))

        # writing to csv file
        with open(filename, 'w') as csvfile:
            # creating a csv writer object
            csvwriter = csv.writer(csvfile)
            # writing the fields
            csvwriter.writerow(fields)
            # writing the data rows
            csvwriter.writerows(rows)