def Clean_Skempi(path):
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


class MutantDataset(pd.DataFrame):
    def __init__(self, data, sep=',', index=None, columns=None, dtype=None,
                 copy=True,):
        '''Initialize subclass from DataFrame instance or csv path.'''
        if type(data) == str:
            data = pd.read_csv(data, sep=sep)
        super(MutantDataset, self).__init__(data=data)

    def _ChainCheck(self, df):
        mutated_chains = set(i[1] for i in df['MutSplit'])
        prot1, prot2 = set(df['Prot1Chain']), set(df['Prot2Chain'])
        if df['NumMutations'] == 1 or len(mutated_chains) == 1:
            return False
        elif bool(mutated_chains-prot1) != bool(mutated_chains-prot2):  # xor
            return False
        else:
            return True
