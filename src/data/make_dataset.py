# -*- coding: utf-8 -*-
import click
import logging
from pathlib import Path
from dotenv import find_dotenv, load_dotenv


def SKEMPItoPandas(SKEMPI_loc):
    '''
    Purpose:
        1. Loads SKEMPI CSV file.
        2. Calculates ddG
        3. For multiple measurements, keeps the median value
        4. Eliminates entries with mutations on both sides of the interface
    Input:
        SKEMPI_loc : Location of SKEMPI CSV file
    Output:
        SKEMPI_df : Pandas dataframe    
    '''
    import pandas as pd
    import numpy as np
    import re
    # fix this
    pd.options.mode.chained_assignment = None  # default='warn'

    # Constants
    R = 1.9872036e-3  # Ideal Gas Constant in kcal

    SKEMPI_df = pd.read_csv(SKEMPI_loc, sep=';')

    # Convert non numeric temperature comments to numeric values. Default is 298K 
    ConvertTemp = lambda x: int(re.search(r'\d+', x).group(0) or 298)
    BadTemps = SKEMPI_df.Temperature.str.isnumeric() == 0
    SKEMPI_df['Temperature'].loc[BadTemps] = SKEMPI_df['Temperature'].loc[BadTemps].map(ConvertTemp)
    SKEMPI_df['Temperature'] = pd.to_numeric(SKEMPI_df['Temperature'], errors='coerce')

    # Drop missing values
    SKEMPI_df.dropna(subset=['Affinity_wt_parsed'], inplace=True)
    SKEMPI_df.dropna(subset=['Affinity_mut_parsed'], inplace=True)

    # Calculate free energies
    SKEMPI_df['dgWT'] = -R*SKEMPI_df['Temperature']*np.log(SKEMPI_df['Affinity_wt_parsed'])
    SKEMPI_df['dgMut'] = -R*SKEMPI_df['Temperature']*np.log(SKEMPI_df['Affinity_mut_parsed'])
    SKEMPI_df['ddG'] = SKEMPI_df['dgWT']-SKEMPI_df['dgMut']

    # Create a key for unique mutations based on PDB and 
    SKEMPI_df['MutKey'] = SKEMPI_df['#Pdb']+'_'+SKEMPI_df['Mutation(s)_PDB']
    # Replace multiple measurements of the same mutation with the group mean
    # May consider grouping by experimental method as well
    SKEMPI_df['ddgMedian'] = SKEMPI_df.groupby('MutKey')['ddG'].transform('median')        
    SKEMPI_df = SKEMPI_df.drop_duplicates(subset=['MutKey', 'Temperature'], keep='first', inplace=False)

    # Flag multiple mutations in the same protein
    SKEMPI_df['NumMutations'] = SKEMPI_df['Mutation(s)_PDB'].str.count(',')+1 

    # Extract Chains and remove cross chain mutations. Chain is the second position in the mutation code
    SKEMPI_df['Prot1Chain'] = SKEMPI_df['#Pdb'].str.split('_').str[1]
    SKEMPI_df['Prot2Chain'] = SKEMPI_df['#Pdb'].str.split('_').str[2]
    SKEMPI_df['MutSplit'] = SKEMPI_df['Mutation(s)_PDB'].str.split(',')

    def ChainCheck(df):
        if df['NumMutations'] == 1:
            CrossChain = False
            return CrossChain
        else:
            Chain = df['MutSplit'][0][1]
            if Chain in df['Prot1Chain']:
                ChainSet = df['Prot1Chain']
            elif Chain in df['Prot2Chain']:
                ChainSet = df['Prot2Chain']
            for i in range(len(df['MutSplit'])):
                Chain = df['MutSplit'][i][1]
                if Chain in ChainSet:
                    CrossChain = False
                else:
                    CrossChain = True
                    break
        return CrossChain

    SKEMPI_df['CrossChain'] = SKEMPI_df.apply(ChainCheck, axis=1)
    SKEMPI_SingleSided = SKEMPI_df[SKEMPI_df.CrossChain == False]

    NumProteins = SKEMPI_SingleSided['#Pdb'].nunique()
    NumMutations = SKEMPI_SingleSided['#Pdb'].count()
    print("There are %s unique single sided mutations in %s proteins" % (NumMutations, NumProteins))             
    return SKEMPI_SingleSided


@click.command()
# fix this patchwork later
#@click.argument('input_filepath', type=click.Path(exists=True))
#@click.argument('output_filepath', type=click.Path())
def main(): #main(input_filepath, output_filepath):
    """ Runs data processing scripts to turn raw data from (../raw) into
        cleaned data ready to be analyzed (saved in ../processed).
    """
    input_filepath = 'data/raw/'
    SKEMPItoPandas(f'{input_filepath}skempi_2.0.csv')  # GENERALIZE!

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
