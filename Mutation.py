
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 13 13:26:33 2019

@author: JRB
"""

# =============================================================================
# PURPOSE: Calculates structural features from sequence
# Generates:
# Volumes for WT and Mutant amino acid and the Volume change
# Hydrophobicity for WT and Mutant for the mutation spot and the Hydrophobicity change
# Charge  for WT and Mutant for the mutation spotand the Charge change
# Flags the presence of proline in WT or Mutant for the mutation spot
# Flags the presence of glycine in WT or Mutant for the mutation spot
# =============================================================================
# DEPENDENCIES
# None
# =============================================================================
# INPUT PARAMETERS
# MutStr: A string in the form WT Amino Acid,Chain ID, Residue Number, Mutant Amino Acid
# =============================================================================


class MutationFeatures():
    # Defines the properties associated with one mutation
	# inherits the class ProteinMethods to defines protein level properties
# =============================================================================
    def __init__(self,MutStr):
	#Basic Information    
        self.WT_AA=MutStr[0]
        self.Mut_AA=MutStr[-1]
        self.Chain=MutStr[1]
        
        #Position in chain
        self.ResidueIndex=MutStr[2:-1]
        
    	#Scores
    	#Volume change associated with mutation. Could be set as either used as is as an integer or binned to form 
        self.VolWT=self.FindVolume(self.WT_AA)
        self.VolMut=self.FindVolume(self.Mut_AA)
        self.VolumeChange=self.VolMut-self.VolWT
        
        self.HPhobWT=self.FindHPhob(self.WT_AA)
        self.HPhobMut=self.FindHPhob(self.Mut_AA)
        self.HPhobChange=self.HPhobMut-self.HPhobWT
        
        #Prolines can introduce conformational changes not modeled by other methods
        self.ProlineWT=self.IsProline(self.WT_AA)
        self.ProlineMut=self.IsProline(self.Mut_AA)
        if (self.ProlineWT==1) or (self.ProlineMut==1):
            self.Proline=1
        else:
            self.Proline=0
            
        #Glycine can break secondary structure
        self.GlycineWT=self.IsGlycine(self.WT_AA)
        self.GlycineMut=self.IsGlycine(self.Mut_AA)
        if (self.GlycineWT==1) or (self.GlycineMut==1):
            self.Glycine=1
        else:
            self.Glycine=0
            
        self.ChargeWT=self.FindCharge(self.WT_AA)
        self.ChargeMut=self.FindCharge(self.Mut_AA)
        self.ChargeChange=self.ChargeMut-self.ChargeWT
        
# =============================================================================        
    def FindVolume(self,AA):
        switcher = {
            "G" : 60.1,
            "A" : 88.6,
            "S" : 89,
            "C" : 108.5,
            "D" : 111.1,
            "P" : 112.7,
            "N" : 114.1,
            "T" : 116.1,
            "E" : 138.4,
            "V" : 140,
            "Q" : 143.8,
            "H" : 153.2,
            "M" : 162.9,
            "I" : 166.7,
            "L" : 166.7,
            "K" : 168.6,
            "R" : 173.4,
            "F" : 189.9,
            "Y" : 193.6,
            "W" : 227.8,
        }
        Volume=switcher.get(AA,"Invalid Amino Acid")
        return Volume
    
    def FindCharge(self,AA):
        switcher = {
            "G" : 0,
            "A" : 0,
            "S" : 0,
            "C" : 0,
            "D" : -1,
            "P" : 0,
            "N" : 0,
            "T" : 0,
            "E" : -1,
            "V" : 0,
            "Q" : 0,
            #assume His is not ionized. Could include a pH dependent switch
            "H" : 0, 
            "M" : 0,
            "I" : 0,
            "L" : 0,
            "K" : 1,
            "R" : 1,
            "F" : 0,
            "Y" : 0,
            "W" : 0,
        }
        Charge=switcher.get(AA,"Invalid Amino Acid")
        return Charge    
    
    def FindHPhob(self,AA):
        #Kyte-Doolittle scale
        switcher = {
            "I" : 0.5,
            "V" : 4.2,
            "L" : 3.8,
            "F" : 2.8,
            "C" : 2.5,
            "M" : 1.9,
            "A" : 1.8,
            "W" : -0.9,
            "G" : -0.4,
            "T" : -0.7,
            "S" : -0.8,
            "Y" : -1.3,
            "P" : -1.6,
            "H" : -3.2,
            "N" : -3.5,
            "D" : -3.5,
            "Q" : -3.5,
            "E" : -3.5,
            "K" : -3.9,
            "R" : -4.5,
        }
        HPhob=switcher.get(AA,"Invalid Amino Acid")
        return HPhob
    
    def IsProline(self,AA):
        if AA == "P": 
            ProlineFlag=1
        else:
            ProlineFlag=0
        return ProlineFlag

    def IsGlycine(self,AA):
        if AA == "P": 
            GlycineFlag=1
        else:
            GlycineFlag=0
        return GlycineFlag    
# =============================================================================
# Testing
# =============================================================================
#MutStr="GI45P"
#Mut=MutationFeatures(MutStr)
#print(vars(Mut))
# =============================================================================
