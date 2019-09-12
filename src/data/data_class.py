
import pandas as pd
import numpy as np


class MutantDataset(pd.DataFrame):
    '''<Subclassed Pandsas DataFrame>
        Given the potential of multiple sources for mutant datasets,
        this calss serves to improve clarity, debugging, and reusability

        Change Mutation(s)_PDB to Mutation(s)_cleaned?
        double check validity of fillna
    '''
    def __init__(self, data, sep=',', index=None, columns=None, dtype=None,
                 copy=True,):
        '''Initialize subclass from DataFrame instance or csv path.'''
        if type(data) == str:
            data = pd.read_csv(data, sep=sep)
        super(MutantDataset, self).__init__(data=data,
                                            index=index,
                                            columns=columns,
                                            dtype=dtype,
                                            copy=copy)

    def Mutations(self, row):
        ''' Returns dictionary of mutation identifiers.
            key: initAA , chain, loc, mutAA
        '''
        keys = ['initAA', 'chain', 'loc', 'mutAA']  # code key
        mut_codes = self.loc[row]['Mutation(s)_cleaned'].split(',')
        unzip_code = zip(*[re.findall('(\d+|.)', mut) for mut in mut_codes])
        mut_dct = dict(zip(keys, unzip_code))
        return mut_dct

    def to_numeric(self, keys):
        '''converts column of single or list of keys to numeric'''
        self[keys] = self[keys].apply(pd.to_numeric, errors='coerce')
        return self[keys]

    def gibbsEq(self, Kd_key, tmp_key='Temperature'):
        '''Gibbs Free Energy = R * Temp * ln(kd)'''
        R = 1.9872036e-3  # Ideal Gas Constant in kcal
        ΔG = R * self[tmp_key] * np.log(self[Kd_key])  # log is ln in np
        return ΔG

    def solve_ddG(self, wild, mutant, tmp_key='Temperature'):
        '''ddG is the changes in affinity upon mutation:
              ddG = dG_Mutant-dG_wild_Type
        '''
        self['dgWT'] = self.gibbsEq(wild, tmp_key)
        self['dgMut'] = self.gibbsEq(mutant, tmp_key)
        self['ddG'] = self['dgMut']-self['dgWT']
        return self

    def solve_ddG2(self, wild, mutant, tmp_key='Temperature'):
        '''ddG is the changes in affinity upon mutation:
              Gibbs Free Energy = R * Temp * ln(kd)  AKA: dG
              ddG = dG_Mutant-dG_wild_Type
        '''
        R = 1.9872036e-3  # Ideal Gas Constant in kcal
        gibbsEq = lambda tmp, kd: R * self[tmp] * np.log(self[kd])  # log is ln in np

        self['dgWT'] = gibbsEq(tmp_key, wild)
        self['dgMut'] = gibbsEq(tmp_key, mutant)
        self['ddG'] = self['dgMut']-self['dgWT']
        return self

    def _ChainCheck(self, df):
        '''Utalizes subtracted sets and xor to identify if
           mutated chains are unique to a single protein
        '''
        mutated_chains = set(i[1] for i in df['MutSplit'])
        prot1, prot2 = set(df['Prot1Chain']), set(df['Prot2Chain'])
        if df['NumMutations'] == 1 or len(mutated_chains) == 1:
            # Single mutant chain is unique to one protien, not cross chain
            return False
        elif bool(mutated_chains-prot1) != bool(mutated_chains-prot2):  # xor
            # Mutated chains are specific to one protein, not cross chain
            return False
        else:
            # A chain remained in both sets after subrtaction,
            # mutations not unique to single protein, cross chain
            return True

    def find_cross_chains(self):
        '''checks if mutation occur on more than one protein'''
        self['Prot1Chain'] = self['#Pdb'].str.split('_').str[1]
        self['Prot2Chain'] = self['#Pdb'].str.split('_').str[2]
        crossChain = self.apply(self._ChainCheck, axis=1)
        return crossChain

    @property
    def _constructor(self):
        return MutantDataset  # Class Name
