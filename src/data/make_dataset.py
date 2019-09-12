# -*- coding: utf-8 -*-
# Install and File Managment:
import click
import logging
from pathlib import Path
from dotenv import find_dotenv, load_dotenv
# Cleaning Datasets:
import pandas as pd
import numpy as np
from data_class import MutantDataset


def Clean_Skempi(path):
    '''
    Purpose:
        1. Loads SKEMPI CSV file.
        2. Calculates ddG
        3. For multiple measurements, keeps the median value
        4. Eliminates entries with mutations on both sides of the interface
    Input:
        path : Location of SKEMPI CSV file
    Output:
        SKEMPI_SingleSided : MutantDataset(pd.DataFrame)
    Note:
        Content and order subject to change with additional datasets.
        It is foreseeable that some steps may occur post combination.
    '''
  # Initialize class
    skempi = MutantDataset(path, sep=';')

  # Convert 'Temperature' comments/str's to numeric values. Default is 298
    skempi['Temperature'] = skempi['Temperature'].str.extract(r'(\d+)')
    skempi['Temperature'] = skempi.to_numeric('Temperature')
    skempi['Temperature'].fillna(value=298, inplace=True)  # 6665-6668 blank

  # Calculate free energies
    dropna_lst = ['Affinity_wt_parsed', 'Affinity_mut_parsed']
    skempi.dropna(subset=dropna_lst, inplace=True)
    skempi = skempi.solve_ddG('Affinity_wt_parsed', 'Affinity_mut_parsed')

  # Median and duplicate ddG/tmp values
    group_keys = ['#Pdb', 'Mutation(s)_PDB']
    skempi['ddgMedian'] = skempi.groupby(group_keys)['ddG'].transform('median')
    skempi = skempi.drop_duplicates(subset=[*group_keys, 'Temperature'],
                                    keep='first', inplace=False)

  # Flag multiple mutations in the same protein
    skempi['MutSplit'] = skempi['Mutation(s)_PDB'].str.split(',')
    skempi['NumMutations'] = skempi['MutSplit'].apply(len)

  # Extract Chains and remove cross chain mutations.
    skempi['CrossChain'] = skempi.find_cross_chains()
    SKEMPI_SingleSided = skempi[skempi.CrossChain == False]
    return SKEMPI_SingleSided


# @click commands removed, extra complexity and CLI unnecessary ATM
def main():
    """ Runs data processing scripts to turn raw data from (../raw) into
        cleaned data ready to be analyzed (saved in ../processed).

        1) clean each dataset to create consistant MutantDataset's
            1a) store indeviduals in ~/data/intermediate
        2) combine into uniform MutantDataset
            2a) store in ~/data/final
    """
  # Paths
    input_filepath = 'data/raw/'
    interim_filepath = 'data/interim/'
    output_filepath = 'data/processed/'

  # 1.0 – Clean skempi
    skempi_final = Clean_Skempi(input_filepath + 'skempi_2.0.csv')
    # Log Info
    NumProteins = skempi_final['#Pdb'].nunique()
    NumMutations = skempi_final['#Pdb'].count()
    print(f'There are {NumMutations} unique single sided'
          f'mutations in {NumProteins} proteins')
  # 1.0a save intermediate to interim
    skempi_final.to_csv(interim_filepath + 'skempi_final.csv')

  # 1.1 – Import Other
    # other_final = Clean_Other(input_filepath + 'other.csv')
  # 1.1a – save intermediate to intermediates
    # Other_final.to_csv('data/intermediate')

  # 2 – Combine datasets
    # code
  # 2a – save final to processed
    # combined.to_csv(data/processed')
    skempi_final.to_csv(output_filepath + 'skempi_final.csv')

    logger = logging.getLogger(__name__)
    logger.info('making final data set from raw data')


if __name__ == '__main__':
    log_fmt = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    logging.basicConfig(level=logging.INFO, format=log_fmt)

    # not used in this stub but often useful for finding various files
    project_dir = Path(__file__).resolve().parents[2]

    # find .env automagically by walking up directories until it's found, then
    # load up the .env entries as environment variables
    load_dotenv(find_dotenv())

    main()
