# Defines the properties associated with one mutation
# Protein.py defines protein level properties


import pandas as pd
import numpy as np

ITSM_loc="C:\\Users\\jbren\\OneDrive\\Documents\\Optimus\\iPTM.csv"
ITSM=pd.read_csv (ITSM_loc,index_col=0,header=0)
Blosum62_loc="C:\\Users\\jbren\\OneDrive\\Documents\\Optimus\\blosumPTSM.csv"
Blosum62 = pd.read_csv (Blosum62_loc,index_col=0,header=0)
VolAA_loc="C:\\Users\\jbren\\OneDrive\\Documents\\Optimus\\volume.csv"
VolAA=pd.read_csv (VolAA_loc,index_col=0,header=0)
class MutationMethods:
    def __init__(self):
        self.AA_list={'A','G','P','K','R','T','Y','L','I','V','D','E','W','F','H','C','M','S','Q','N'}

        
    def getSubMatrix(self,mutcode,SM):          
        WT_AA=mutcode[0]
        mut_AA=mutcode[-1]
        if SM == 'ITSM':
            SMScore=np.log(ITSM.at[WT_AA,mut_AA])
        elif SM == 'Blosum':
            SMScore=np.log(Blosum62.at[WT_AA,mut_AA])
        return SMScore
        
    def getSMScore(self,mutset,SM):
        SMScore=0
        for mut in mutset:
            self.SMScore=self.getSubMatrix(mut,SM)+SMScore
        return SMScore
    
    def HasProline(self,mutset):
        HasProline=0
        for mut in mutset:
            WT_AA=mut[0]
            mut_AA=mut[-1]
            if (WT_AA == 'P') or (mut_AA =='P'):
                HasProline=1
            self.HasProline=HasProline
        return HasProline
        
    def HasGlycine(self,mutset):
        HasGlycine=0
        for mut in mutset:
            WT_AA=mut[0]
            mut_AA=mut[-1]
            if (WT_AA == 'G') or (mut_AA =='G'):
                HasGlycine=1
            self.HasGlycine=HasGlycine
        return HasGlycine
    
    def SmallToLarge(self,mutset):
        Small=999999
        for mut in mutset:
            WT_AA=mut[0]
            mut_AA=mut[-1]
            vol_change=VolAA.at[mut_AA,'vol']-VolAA.at[WT_AA,'vol']
            if vol_change<Small:
                Small=vol_change
            self.Small=Small
        return Small 

    def LargeToSmall(self,mutset):
        Large=0
        for mut in mutset:
            WT_AA=mut[0]
            mut_AA=mut[-1]
            vol_change=VolAA.at[mut_AA,'vol']-VolAA.at[WT_AA,'vol']
            if vol_change>Large:
                Large=vol_change
            self.Large=Large
        return Large
    
    # returns the secondary structure code from a DSSP object
    def calcSS(self,dssp,ch,mutset):
        self.SS={}
        for mut in mutset:
            AAindex=int(mut[1:-1])
            DSSPkey = (ch,(' ',AAindex,' '))
            SS_AA=dssp.property_dict[DSSPkey][2]
            self.SS[AAindex] = SS_AA
        return self.SS
    
    def relSASA(self,dssp,ch,mutset):
        self.relSASA={}
        for mut in mutset:
            AAindex=int(mut[1:-1])
            DSSPkey = (ch,(' ',AAindex,' '))
            relSASA_AA =dssp.property_dict[DSSPkey][2]
            self.relSASA[AAindex] = relSASA_AA
        return self.relSASA
           
    #Heinkoff weighting -not finished. Should convert to matrix and use numpy
    # def Heinkoff(self,MSA):
    #     l=len(MSA)
    #     i=0
    #     while i<len(MSA):
    #         j=0
    #         wt=[]
    #         while j<len(MSA[i]):
    #             s=set(MSA[:,j])
    #             r=len(s)

    #             wt_seq=0
    #             for aa in s:
    #                 n_aa=MSA[:,j].count(aa)
    #                 wt_seq=wt_seq+1/(r*n_aa*l)
    #        #     wt.append(wt_seq)
    #             j+=1
    #             score=0
    #             for aa in self.amino_acids:
    #                 score=score+wt_seq
    
    #         i+=1
    #     return wt
    
    