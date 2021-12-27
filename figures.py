"""
--------------------------------------------------------------------------------
Description:

Produce logos from an arbitrary numpy array position weight matrix of (4 x k-mer length)

Written by Aidan Riley, Juan Montezo, Eric South
--------------------------------------------------------------------------------
"""
from pathlib import Path
import glob
import os
import pandas as pd
import numpy as np
import logomaker


def generate_dummy_data(k_len=10):
    """
    Generate a dummy position weight matrix with 4 rows (ATCG) of length k.

    :param k_len: desired sequence length.
    :return wm: numpy array (4, k) consisting of probabilities summing to one.
    """
    wm = np.empty((4, k_len)).T  # Instantiate weight matrix
    for i in range(k_len):
        _ = np.random.random(size=4)  # Generate random probabilities
        foobar = _ / _.sum()
        wm[i] = foobar
    return wm.T


def data_import(file_name):
    """
    Load .csv file and reformat values as a numpy array.
    
    :param x: file path for csv file.
    :return y: numpy array (4, k).
    """ 
    data = pd.read_csv(file_name).values
    return data


def logomaker_reformat(weight_array):
    """
    Reformat numpy array to Dataframe with proper labels for logomaker calls.
    
    :param weight_array: numpy array (4, k) with position-specific probabilities.
    :return df: Dataframe with rows as position #s and columns as A, T, C, and G.
    """  
    df = pd.DataFrame(weight_array.T, columns = ['A','G','C', 'T'])
    # df = df.loc['A', 'T']
    # df = pd.DataFrame(weight_array.T, columns = ['G','C'])
    df.index.name = 'pos'

    return df




def generate_logo(data, title='Some File'):
    """
    Generate logos using the logomaker package.
    
    :param x: Dataframe of position-specific weighted probabilities
    :return y: figure
    """ 
    # Create logo object
    data_logo = logomaker.Logo(data,
                              shade_below=.5,
                              fade_below=.5,
                              font_name='Arial Rounded MT Bold')
    # Style using logo methods
    data_logo.style_spines(visible=False)
    data_logo.style_spines(spines=['left', 'bottom'], visible=True)
    data_logo.style_xticks(rotation=0, fmt='%d', anchor=0)
    # Style using axes methods
    data_logo.ax.set_ylabel('Saliency', labelpad=-1)
    data_logo.ax.set_title(title[:-4], size=18)
    data_logo.ax.xaxis.set_ticks_position('none')
    data_logo.ax.xaxis.set_tick_params(pad=-1)
    # Highlight functional regions of ARS1
    # data_logo.highlight_position_range(pmin=4, pmax=6, color='cyan')


#TODO: Generate histogram that compares threhold values with other features.
def histogram():
    pass


# LowResponse: switches picked with a score up to 0.10, 0.50, or 0.90
# HighResponse: switches picked with a score no less than 0.10, 0.50, 0.90
# BinnedResponse: switches picked with an average bin of 0.125, 0.525, 0.925 (range 0.05)



def main():
    """
    Main function which uploads data and generates logos.
    
    :param x: ...
    :return y: ...
    """
    path = Path('/Users/Shockwing/Dropbox/Boston/BU_curriculum/BE_562/project/BE-562-Toehold-Project/data/211207_2/')
    # p = Path('./data') / 'k_mers_211203'  # Specify relative file path
    # path = p.resolve()  # Generate absolute path
    os.chdir(str(path))
    
    for file_name in glob.glob('*.csv'):  # Extract all .xlsx files from cwd.
        data = data_import(file_name)
        data = logomaker_reformat(data)
        data_at = data[['A', 'T']]
        data_gc = data[['G', 'C']]
        generate_logo(data_at, title=f'at_{file_name}')
        generate_logo(data_gc, title=f'gc_{file_name}')
        
        
if __name__ == "__main__":
    app()