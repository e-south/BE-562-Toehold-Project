'''
--------------------------------------------------------------------------------
<BE562 Project: Toehold Switch Motif Proposal>

Module will load the dataset from (Angenet-Mari et. al 2020) and save a file of high-scoring and low-scoring switch sequences

Written by Aiden Riley, Juan Montezo, Eric South
--------------------------------------------------------------------------------
'''
# Supporting Libaries
import pandas as pd
#Load in the dataset as a pandas dataframe
#Make sure the datafile is in the working directory
def main():
    dataset = pd.read_csv('Toehold_Dataset_Final_2019-10-23.csv')

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
    low_score = dataset[dataset['ON_OFF'] < 0.05]
    low_score.to_csv('low_score.csv')

    high_score = dataset[dataset['ON_OFF']> 0.85]
    high_score.to_csv('high_score.csv')








if __name__ == '__main__':
    main()