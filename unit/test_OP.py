"""
-----------------------------
Optimus binding unit testing
-----------------------------

Version: 0.0.1

Authors: Sang Young Noh (Please add your name once you edit this file)


Using pytest to do the following

- Test main binaries (foldx, ialign) 

- Test mutation file 

--------

For this purpose, I have copied 3 pdb files as input in this folder as follows: ["1DAN.pdb","1DQJ.pdb","1DVF.pdb"].

For pytest to register each function to test, each function of relevance needs to start with test_something.py. 

I've added the sys module so that it allows easy importing of the mutation python file.



"""

# Standard modules

import os
import sys
import pytest
import subprocess
sys.path.append(os.path.abspath(os.path.join('..', 'src'))) # Reference to the src path

# Importing mutation module from src

import mutation
import blast
from mutation import foldx_ialign # python code containing the foldx mutation function
from blast import BlastScoring # TODO - need to rename this module

IalignPath = "/home/oohnohnoh1/Desktop/ACADEMIA/Papermaking/OPTIMUS_BIND/ialign/bin/ialign.pl" # Path to Ialign binary 
FoldxPath = "/home/oohnohnoh1/Desktop/ACADEMIA/Papermaking/OPTIMUS_BIND/FoldX/foldx"  # Path to Foldx binary

DataFrame = foldx_ialign.SKEMPItoPandas('skempi_v2.csv') # Imported from mutations 
PDBList = ["1DAN.pdb","1DQJ.pdb","1DVF.pdb"] # PDB list in the test folder

mutationFolder = ' ' # TODO
FixedWTFolder = ' ' # TODO

#pytest.fixture()
def test_callialign():
	"""
	Purpose:
	-------
	
	We would like to test if the pandas table parser  (SKEMPItoPandas) works accordingly
	
	"""
	WorkVal = False 

	DataFrame = foldx_ialign.SKEMPItoPandas('skempi_v2.csv')
	if DataFrame is not None:
		WorkVal = True
	assert(WorkVal == True)
	
def test_callfoldx():
	"""
	Purpose:
	-------
	Test function for callfoldx in the mutations python file
	
	"""
	val = ['F', 'F', 'F']
	for index, pdb in enumerate(PDBList):
		try:
			foldx_ialign.callfoldx(pdb)
			val[index] = 'T'
		except ValueError:
			print ("FAIL")
	# Test
	assert(val == ['T', 'T', 'T'])

def test_psiBlastScoring():
	pass

