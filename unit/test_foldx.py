"""
Optimus binding unit testing
-----------------------------

Version: 0.0.1

Authors:

etc. etc. 

Using pytest to do the following

- Test main binaries (foldx, ialign) 

- Test mutation file 

For this purpose, I have copied 5 pdb files as input in this folder as follows: [] 


"""

import os
import sys
import pytest
import subprocess

sys.path.append(os.path.abspath(os.path.join('..', 'src'))) # Reference to the src path

# importing mutation module from src

import mutation
from mutation import foldx_ialign # python code containing the foldx mutation function

IalignPath = "/home/oohnohnoh1/Desktop/ACADEMIA/Papermaking/OPTIMUS_BIND/ialign/bin/ialign.pl" # Path to Ialign binary 
FoldxPath = "/home/oohnohnoh1/Desktop/ACADEMIA/Papermaking/OPTIMUS_BIND/FoldX/foldx"  # Path to Foldx binary

DataFrame = foldx_ialign.SKEMPItoPandas('skempi_v2.csv') # Imported from mutations 
PDBList = ["1DAN.pdb","1DQJ.pdb","1DVF.pdb"] # PDB list in the test folder

mutationFolder = ' ' # TODO
FixedWTFolder = ' ' # TODO

#pytest.fixture()
#def test_callialign():
#	"""
#	Documentation
#	"""
#	DataFrame = SKEMPItoPandas('skempi_v2.csv')
#	for pdb in PDBList:
#		assert(foldx_align.GenerateMutations(DataFrame, pdb, '.')) # Run the generateMutations
#		# TODO
		

#def test_callialign():
#	"""
#	-------------
#	Documentation
#	-------------
#	"""
#	pass


def test_callfoldx():
	"""
	-------------
	Documentation
	-------------

	Test function for callfoldx in the mutations python file

	=> Passes all test cases - 11/08/2019 - SYN
	"""
	val = ['F', 'F', 'F']
	for index, pdb in enumerate(PDBList):
		try:
			foldx_ialign.callfoldx(pdb)
			val[index] = 'T'
		except ValueError:
			print ("FAIL")
	assert(val == ['T', 'T', 'T'])


	
