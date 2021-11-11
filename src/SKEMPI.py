# -*- coding: utf-8 -*-
"""
Created on Sat Jul 20 02:41:31 2019

@author: Jeeves
"""
# =============================================================================
# Purpose: 
#    1 Parses skempi csv file to return dataframe
#    2 Coverts Kd binding constants to ddG free energies
#    3 Summarizes multiple tests on same mutation by using the median of the results
#    4 Fills in missing data with standard values
#    5 Removes rows where mutations are on both side of the complex
# Input: The skempi_v2 csv file from the SKEMPI website in unaltered form
# Output: A dataframe containing all the relevant information from the SKEMPI database
# Dependencies: No dependencies on other scripts
# Assumptions: skempi_v2.csv is in the current directory    
# =============================================================================

import pandas as pd
import numpy as np
import re

class MutantDataSet(pd.DataFrame):
    def __init__(self, data, sep=',', index=None, columns=None, dtype=None, 
                 copy=True,):
        # Initialize subclass from DataFrame instance or csv path.
        # __init__ is always run when the class is first called
        if type(data)==str:
            data=pd.read_csv(data, sep=sep)
        super(MutantDataSet, self).__init__(data=data,
                                            index=index,
                                            columns=columns,
                                            dtype=dtype,
                                            copy=copy)

    def Mutations(self, row):
        # Returns dictionary of mutation identifiers.
        keys = ['initAA', 'chain', 'loc', 'mutAA']  # code key
        mut_codes = self.loc[row]['Mutation(s)_cleaned'].split(',')
        unzip_code = zip(*[re.findall('(\d+|.)', mut) for mut in mut_codes])
        mut_dct = dict(zip(keys, unzip_code))
        return mut_dct
  
    def to_numeric(self, keys):
        # converts column of single or list of keys to numeric values
        self[keys] = self[keys].apply(pd.to_numeric, errors='coerce')
        return self[keys]
          
    def gibbsEq(self, Kd_key, tmp_key='Temperature'): 
        R = 1.9872036e-3  # Ideal Gas Constant in kcal
        dG = -R * self[tmp_key] * np.log(self[Kd_key]) #log is ln in np
        return dG
    
    def solve_ddG(self, wild, mutant, tmp_key='Temperature'):
        self['dgWT'] = self.gibbsEq(wild, tmp_key)
        self['dgMut'] = self.gibbsEq(mutant, tmp_key)
        self['ddG'] = self['dgWT']-self['dgMut']
        return self
    
    def grouped_avg(self, group_keys, avg_key):   
        # should rename to grouped_med... 
        averaged = self.groupby(group_keys)[avg_key].transform('median')
        return averaged  # returns series
    
    def ChainCheck(self,df):
    # operates on each row (called by apply)
    # if only one mutation skip and go to next row
        if df['NumMutations'] == 1:
            crossChain = False
            return crossChain
        else:
           # find if the first partner is in the mutation
           Chain = df['MutSplit'][0][1]
           if Chain in df['Prot1Chain']:
               ChainSet = df['Prot1Chain']
           # if not use the second partner    
           elif Chain in df['Prot2Chain']:
               ChainSet = df['Prot2Chain']
           # if not use the second partner 
           for i in range(len(df['MutSplit'])):
               Chain = df['MutSplit'][i][1]
               if Chain in ChainSet:
                   crossChain = False
               else:
                   crossChain = True
                   break
        return crossChain
   
    def find_cross_chains(self):
    # Calls ChainCheck to find rows that have mutations on both side of the protein
    # Create 2 new columns with the chain IDs
    # note that the Prot1Chain can contain more than one chain
    # A->BCD A binds to the complex of BCD
        self['Prot1Chain'] = self['#Pdb'].str.split('_').str[1]
        self['Prot2Chain'] = self['#Pdb'].str.split('_').str[2]
        # Apply function Chaincheck to each row (axis=1) 
        crossChain = self.apply(self.ChainCheck, axis=1)
        return crossChain
  
    @property
    def _constructor(self):
        return MutantDataSet # Class Name


#Driver function
def clean_Skempi(path):
# Initialize class
    skempi = MutantDataSet(path, sep=';')

# Convert non-numeric temperature comments to numeric values. Default is 298K 
    skempi['Temperature'] = skempi['Temperature'].str.extract(r'(\d+)')
    skempi['Temperature'] = skempi.to_numeric('Temperature')
    skempi['Temperature'].fillna(value=298, inplace=True) #6665-6668 blank  ### TOGGLE ME ###
  
# Calculate free energies
    dropna_lst = ['Affinity_wt_parsed','Affinity_mut_parsed'] #, 'Temperature']
    skempi.dropna(subset=dropna_lst, inplace=True)
    skempi = skempi.solve_ddG('Affinity_wt_parsed', 'Affinity_mut_parsed')

# Find duplicate ddG/tmp values and replace with the median
    group_keys = ['#Pdb', 'Mutation(s)_PDB']
    skempi['ddgMedian'] = skempi.groupby(group_keys)['ddG'].transform('median')
    #skempi['ddgMedian'] = skempi.grouped_avg(group_keys, 'ddG')
    skempi = skempi.drop_duplicates(subset=[*group_keys,'Temperature'], 
                                    keep='first', inplace=False)

# Flag multiple mutations in the same protein
    skempi['MutSplit'] = skempi['Mutation(s)_PDB'].str.split(',')
    skempi['NumMutations'] = skempi['MutSplit'].apply(len)

# Extract Chains and remove cross chain mutations. 
    skempi['CrossChain'] = skempi.find_cross_chains()
    SKEMPI_SingleSided = skempi[skempi.CrossChain == False]
    return SKEMPI_SingleSided

# =============================================================================
# End of class definition
# driver program
# =============================================================================
skempi_final = clean_Skempi('skempi_v2.csv')
NumProteins = skempi_final['#Pdb'].nunique()
NumMutations = skempi_final['#Pdb'].count()
print("There are %s unique single sided mutations in %s proteins" % 
      (NumMutations, NumProteins))  

#1a – store skempi in ~/data/intermediate
##skempi_final.to_csv('~/data/intermediate')

#1 – clean Other 
  #other = MutantDataSet('other.csv')
  #1a – store Other in ~/data/intermediate


#2 – combine 
  #2a – store in ~/data/final  
