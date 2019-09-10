"""
----------------
| OPTIMUS BIND | 
----------------

------------------------------------------------------------------------------------------------

Version: 0.0.1

Last Updated: 04/09/2019

Description: Functions which generate a database to plug into the machine learning framework 
             (Needs to be expanded)

Contributors: Sang Young Noh. Jeffrey Brender, .. 

Contact: sangyoung123@googlemail.comes
-------------------------------------------------------------------------------------------------
"""
import sys
import os
from os import listdir
from os.path import isfile, join
import glob


import click
import logging
from pathlib import Path
from dotenv import find_dotenv, load_dotenv
import re

# Subprocess modules for calling on making the 

import subprocess
from subprocess import call

# Scipy Stack

import numpy as np
import pandas as pd
import scipy
import seaborn as sns
import matplotlib.pyplot as plt

# tqdm

import tqdm as tqdm 

# XML Parser

import xml.etree.ElementTree as ET

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
    SKEMPI_df['MutCleanSplit'] = SKEMPI_df['Mutation(s)_cleaned'].str.split(',')

	# SYN added - Added a pdb name column to make it easier to identiy pdb when it comes to implementing
	# mutations
	
    NAME = [] 
    for pdbname in SKEMPI_df['#Pdb']:
        name = pdbname.split('_')[0]
        NAME.append(name)
    SKEMPI_df['NAME'] = NAME

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

DataFrame = SKEMPItoPandas('skempi_v2.csv')


# -----------------------
# Calling perl for ialign
# -----------------------

# Ialign subprocess. Example input
# ialign.pl -w output 1lyl.pdb AC 12as.pdb AB

IalignPath = "/home/oohnohnoh1/Desktop/ACADEMIA/Papermaking/OPTIMUS_BIND/ialign/bin/ialign.pl"

# -------------
# Calling foldx
# -------------

# FoldX subprocess. Example input:
# foldx --command=Optimize --pdb=example.pdb

FoldxPath = "/home/oohnohnoh1/Desktop/ACADEMIA/Papermaking/OPTIMUS_BIND/FoldX/foldx"

#subprocess.call(IalignPath)

def align(PDB, mutfolder, output):
	"""
	Purpose: 
	
	
	#Run "run_align.pl" to generate interface multiple sequence alignment file

	#run_align.pl DIMER_PDB_FILE IS-SCORE_CUTOFF ALIGNMENT_OUTPUT_FILE

	#Example:

	#run_align.pl demo/1A22.pdb 0.5 demo/align.out
	
	Parameters
	----------
	

	"""
	print ("Running run_align...")
	alignPath = "/home/oohnohnoh1/Desktop/ACADEMIA/Papermaking/OPTIMUS_BIND/BindProfX/bin/XBindProf/run_align.pl"
	try: 
		p = subprocess.Popen("{} {}/{} {} {}/{}".format(alignPath, mutfolder, PDB, score, mutfolder, output), shell = True) # Ok this works
		stdout, stderr = p.communicate()
	except subprocess.CalledProcessError as e:
		print ("ERROR: Cannot read run_align properly. Please check the input file/chain/pdb files, or check the path to the ialign binary")
	else:
		print ("Running run_align sucessfully..")
	

def callialign(folder, pdb1, chain1, pdb2, chain2, mut, Ialignpath = None): # Default None for Ialignpath for now
	"""
	Purpose:
	
	Function to call ialign with subprocess

	Parameters
	----------
	folder: path/to/directory 
	    Path to directory where the WT and mutations are stored
	pdb1: pdb file name
	    The string of the pdb file
	chain1: chain name
	    Chain in pdb1 
	pdb2: pdb file name
	    The string of the pdb file
	chain2: chain name
	    Chain in pdb1 
	Ialignpath: placeholder
	    ---
	"""
	# TODO - add timeit 
	try:
		import os
		from os import path
	except ImportError:
		print ("ERROR: Cannot import Python modules")
	IalignPath = "/home/oohnohnoh1/Desktop/ACADEMIA/Papermaking/OPTIMUS_BIND/ialign/bin/ialign.pl"
	print ("Running ialign...")
	try: 
		p = subprocess.Popen("{} --w {} {} {} {} {} >> {}_align_output.dat".format(IalignPath,folder, pdb1, chain1, pdb2, chain2, mut), shell = True) # Ok this works
		stdout, stderr = p.communicate()
	except subprocess.CalledProcessError as e:
		print ("ERROR: Cannot read ialign properly. Please check the input file/chain/pdb files, or check the path to the ialign binary")
	else:
		print ("Running ialign sucessfully..")
	print ("Storing the ialign scores in an array..." )

	# TODO - move ialign output to folder, and write the dimer.lst file for run_align to read
	# align(PDB, mutfolder, output)
	try:
		ialign_data = function_to_read_ialign_profiles("{}_align_output.dat".format(mut), mut) # placeholder
	except ValueError: # It might not be ValueError - need to check
		print ("ERROR: cannot read ialign score data properly")
	
		
def callfoldx(pdb, foldxpath = None): # Default None for Ialignpath for now
	"""
	Purpose:
	
	Function to call folx with subprocess
    
	Parameters
	----------
	pdb: pdb file name
	    The string of the pdb file
    foldxpath: path/to/binary
	    The path to the foldx binary
	"""
	# We need to make sure the rotabase.txt file is in the pdb folder
	try:
		import errno
		import os
		from os import path
	except ImportError:
		print ("..")
	FoldxPath = "/home/oohnohnoh1/Desktop/ACADEMIA/Papermaking/OPTIMUS_BIND/FoldX/foldx"
	print ("---------------------")
	print ("... Running foldx ...")
	print ("---------------------")
	try:
		pdbstring = "--pdb={}".format(str(pdb))
		p = subprocess.Popen([FoldxPath, "--command=Optimize", pdbstring], stdout = subprocess.PIPE) # Need to check if this works -works
		stdout, stderr = p.communicate()
	except subprocess.CalledProcessError as e:
		print ("ERROR: Cannot run foldx properly. Please check the input file/chain/pdb files, or check the path to the foldx binary")             
	else:
		print ("------------------------------")
		print ("... Ran foldx successfully ...")
		print ("------------------------------")
		
	# -------------------------------------------------
	# Check if the optimized files have been produced |
	# -------------------------------------------------
	
	output_FX = "OP_{}.fxout".format(pdb.split('.')[0]) # FX name output
	output_PDB = "Optimized_{}.pdb".format(pdb.split('.')[0]) # Optimized PDB name output
	rotabase = "rotabase.txt" # rotabase.txt required for foldx

	# Exception conditionals for the rotabase file

	if os.path.isfile(rotabase) is True:
		print ("rotabase.txt is present.")
	else:
		raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), rotabase)
	
	# Exception conditionals for the FX file
	if os.path.isfile(output_FX) is True:
		print ("{} has been produced successfully".format(output_FX))
	else:
		raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), output_FX)
	
	# Exception conditionals for the PDB file 
	if os.path.isfile(output_PDB) is True:
		print ("{} has been produced successfully".format(output_PDB))
	else:
		raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), output_PDB)

def ReadEnergy(PATH):
	"""
	Purpose:
	
	placeholder

	Parameters
	----------
	PATH:
	   The path to where the fxout files are stored
	"""
	# Find fxout files and parse the results within them
	try:
		import itertools
		import pandas as pd
		import numpy as np
	except ImportError:
		print ("Error")
		
	COLNAMES = ["PDB", "Total", "BackHbond","SideHbond", "Energy_VdW", "Electro", "Energy_SolvP", "Energy_SolvH", "Energy_vdwclash", "Entropy_sidec", "Entropy_mainc", "A", "B", "cis_bond energy_torsion", "energy_torsion", "backbone_vdwclash", "helix dipole", "C", "disulfide", "kn_electrostatic", "partial covalent interactions", "Energy_Ionisation", "F", 'IS-score', 'P-value', 'Z-score', 'Number of aligned residues', 'Number of aligned Contacts', 'RMSD', 'Seq Identity']

	#iAlignEntries = ['IS-score', 'P-value', 'Z-score', 'Number of aligned residues', 'Number of aligned Contacts', 'RMSD', 'Seq Identity']

	OutputVec = [] # Output list to store the datafiles which we will use to convert into a pandas table
	FxoutFileArray = [] # pathway list for the fxout files, which we will store and read one by one 
	for file in os.listdir(PATH): # List the fxout files in the directory, and store them in the array 
		if file.endswith(".fxout"):
			FileLocation = os.path.join(PATH, file)
			FxoutFileArray.append(FileLocation)
			
	for file in FxoutFileArray:
		PDBName = (file.split('/'))[-1] # From the filepath, take the pdb fxout name 
		TextFile = open(str(file), "r")   # Open the file 
		data = []
		for num, line in enumerate(TextFile):
			if num == 4: # Optimized values
				data.append(PDBName)
				for col in line.split("\t"):
					data.append(col)
		OutputVec.append(data)
	return OutputVec

	
def function_to_read_ialign_profiles(outputdata, PDBname, foldPandasInput = None):
	"""
	Purpose:

	- IS-Score - brief explanation  
	- P-Value - ditto 
	- Z-Score - ditto 
	- Number of aligned residues - ditto 
	- Number of aligned contacts - ditto
	- RMSD - ditto
	- Seq Identity - ditto
	
	Parameters
	----------
	outputdata: pdb file name
	    Placeholder
	PDBname: path/to/binary
	    PLaceholder
	foldPandasInput: 
	    Placeholder
	"""	
	Entries = ['IS-score', 'P-value', 'Z-score', 'Number of aligned residues', 'Number of aligned Contacts', 'RMSD', 'Seq Identity']
	ialignOutput = open(str(outputdata), "r")
	ialignLines = ialignOutput.readlines()
	datBlock = None
	output = []
	#print (ialignLines)
	# Concatenate the data block where we have the data
	for index, entry in enumerate(ialignLines):
		if 'IS-score' in entry:
			datBlock = ialignLines[index:index+4]
			break
	#print (datBlock)
	datBlock  = ','.join(datBlock)
	datBlock = datBlock.split(',')
	print (datBlock)
	
	for entry in datBlock:
		output.append([entry.split(' = ')[0].strip(), float(entry.split(' = ')[1])])
	print ("The WT vs mutation ialign results for {} is as follows:".format(PDBname))
	return output
		   
def GenerateMutations(DataFrame, PDB, PATH):
	"""
	Purpose:
	
	This function returns the mutated pdb protein files 
	from skempi_v2 database (https://life.bsc.es/pid/skempi2/). 

	Both single mutations and multiple comma separated mutations 
	are taken in to account. 

	If there are multiple mutation indices for the same protein, 
	then this will generate multiple pdb files.

	Parameters
	----------
	DataFrame: pandas table 
	    The pandas table to read_csv
	PDB: str
	    The string of the pdb file
	"""
	try:
		from Bio.PDB.PDBIO import PDBIO
		from Bio.PDB.PDBParser import PDBParser
		from Bio.Data.IUPACData import protein_letters
		from Bio.SeqUtils.ProtParam import ProteinAnalysis
		from Bio.PDB.Polypeptide import PPBuilder  
		from Bio.PDB.Polypeptide import standard_aa_names # Standard amino acid names - https://biopython.org/DIST/docs/api/Bio.PDB.Polypeptide-module.html 
		from Bio.PDB.Polypeptide import aa1 #  aa1 = 'ACDEFGHIKLMNPQRSTVWY'
		from Bio.PDB.Polypeptide import aa3 #  aa3 = ['ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE',... ]
		import tqdm as tqdm # tqdm - useful for estimating computing times for long for loops 

	except ImportError:
		print ("ERROR: Need to check Biopython imports!")

	# Before running anything, call foldx on the WT to get the optimized structure to mutate
	title = PDB.split('.')
	name = title[0]
	callfoldx(PDB) # Call FoldX on the WT 

	# Path to where the WT PDBs are stored
	WTArray = []
	nameArray = []
	
	for file in os.listdir(PATH): # List the fxout files in the directory, and store them in the array 
		if file.endswith(".pdb"):
			FileLocation = os.path.join(PATH, file)
			WTArray.append(FileLocation)
			nameArray.append(file)

	# Subprocessing block for WT
	subprocess.Popen("mkdir {}".format(name), shell = True) # Make directory
	subprocess.Popen("mv OP_{}.fxout {}/.".format(name,name), shell = True)	 # Move optmiized fxout file to directory	   
	subprocess.Popen("mv Optimized_{}.pdb {}".format(name,PDB), shell = True) 	# Rename file from Optimized_PDB.pdb to the same name as the original file to make our lives easier
 		   

	MutationSpecies = [] # List to store the names of the mutated speices 
	AminoAcidListDict = {} # Dictionary to assign alpabetical letters to amino acids
	for index, code in enumerate(standard_aa_names):
		AminoAcidListDict[aa1[index]] = aa3[index] # Building the mutation dictionary for each code
	parser = PDBParser(PERMISSIVE=1) # Standard PDB parser
	PDBList = set()	

	for pdb in DataFrame['#Pdb']:
		pdbname = pdb.split('_')[0]
		string = "{}.pdb".format(pdbname)
		PDBList.add(string)
		
	if PDB not in PDBList:
		raise Exception("The PDB is not in the SKEMPI list") # Not in the PDB list we expect - i.e. from the SKEMPI list
	
	# Search for PDB mutations that contain the PDB string - e.g. the 1CSE mutations will have the format 1CSE_E_I
# where it indicates the mutations were made in the 1CSE E and I chains

	MutationList = DataFrame.loc[(DataFrame['NAME'] == PDB.split('.')[0])] # This should get the PDB mutations 
	MutationList = MutationList.reset_index()
	print (MutationList)
	
	# Make a dictionary (hash map) with the mutation name and the residue lists to change
	for index, entry in MutationList.iterrows():
		structure = parser.get_structure(str(title[0]),PDB) # reset structure each time 
		model = structure[0]  # Switch back to the unchanged one 
		for mut in entry['MutCleanSplit']:
			initAA, chain, loc, mutAA = re.findall('(\d+|.)', mut)
			# Check we are reading the right residue and index
			assert(model[chain][int(loc)].resname == AminoAcidListDict[str(initAA)]) # This will check that the model is the unmutated pdb
			print ("Mutating {} on index {} of chain {} to {}".format(AminoAcidListDict[str(initAA)], chain, loc , AminoAcidListDict[str(mutAA)]))
			model[chain][int(loc)].resname = AminoAcidListDict[str(mutAA)] # This command replaces the nonmutated species into the mutated one
			assert(model[chain][int(loc)].resname == AminoAcidListDict[str(mutAA)]) # This will check that the mutation was successful
		mutanttotalstring = '_'.join(entry['MutCleanSplit'])
		mutatedname = "{}_{}_{}.pdb".format(entry['#Pdb'], mutanttotalstring, index)
		MutationSpecies.append(mutatedname)
		io = PDBIO(structure)
		io.set_structure(model)
		io.save(mutatedname) # This should print out the name of protein, the mutaton list, and the index on the pandas file 
		print ("Produced new mutation PDB file {}".format(mutatedname)) # Printing out sign to say the pdb was produced

	# Call foldx on the mutatied species

	print (" -----------------------------------------------")
	print ("The following mutant species are to be optimized")
	print (" -----------------------------------------------")

	for mutant in MutationSpecies:
		print ("PDB file: {}".format(mutant))
	ANS = []

	for species in MutationSpecies:
		callfoldx(species) # Call FoldX on each mutated species
		subprocess.Popen("mv {} {}/.".format(species, name), shell = True)
		subprocess.Popen("mv OP_{}.fxout {}/.".format(species.split(".")[0], name), shell = True)
		subprocess.Popen("mv Optimized_{}.pdb {}/.".format(species.split(".")[0], name), shell = True)
		ANS.extend(ReadEnergy("{}/".format(name)))
		print (ANS)
		
	print ("Finished Optimization")
	print ("Running Ialign..")

	# This part needs to be fixed
	
#	OptimizedMutationList = ["Optimized_{}".format(mutant) for mutant in MutationSpecies] # Make new list for Optimized_PDB.pdb names

#	for species in OptimizedMutationList:
#		chain = species.split('_')[2]
#		output_name = species.split(".")[0]
#		print ("Running ialign for mutant PDB: {} with {}".format(species, name))
#		callialign('output_{}'.format(name), PDB, chain, species, chain, output_name)

		
def mapped_index(pdb, chain, index, basis='FASTA_index'):
	names = ["Residue", "Chain", "FASTA_index", "PDB_index"]
	map_data = pd.read_csv(f'PDBs/{pdb}.mapping', sep=r"\s*", 
                           header=None, names=names)
	map_data = map_data.loc[map_data['Chain'] == chain]
	map_data = map_data.set_index(basis)
	output = names[-2:][basis=='FASTA_index']
	return map_data[output][index]

	# Now that we have a dictionary, we need to make sure the mutations can be read whether a single
	# mutation or a comma separated mutation.
	
	#for index in range(0, len(MutationList)-1):
	#	MutationDict[MutationList['#Pdb'][index]] = MutationList['MutSplit'][index]
	# Parsing mutation output 
	# If we have multiple mutations that are comma separated,
	# then we will use this 
	#mutations = test_mutation.split(',')  # Thanks to the suggestion by Thomas. J. Card of the Optimus Prime project 
	#for mut in mutations:
	#	initAA, chain, loc, mutAA = re.findall('(\d+|.)', mut)
	


"""
Some notes on psi-blast for my own sake

Basic Local Alignment Search Tool (BLAST)  is a sequence similarity search program
used to compare a user's query to a database of sequences.  Given a DNA
or amino acid sequence, the BLAST heuristic algorithm finds short matches 
between two sequences and attempts to start alignments from these "hot spots". 
BLAST also provides statistical information about an alignment such as the "expectation"
value. Note that BLAST is not a single program, but a family of programs. 

All BLAST programs search for match between sequences, but there is a specialized 
BLAST program for each type of sequence search. 

---------------------------------------------------------------------------------
BLAST is one of the most widely used bioinformatics research tools,  since it has
several applications, here is a list of typical BLAST applications
---------------------------------------------------------------------------------

1. Following the discovery of a previously unknown gene in one species, search other genomes 
   to see if other species carry a similar gene.

2. Finding functinoal and evolutionary relationships between sequences

3 Search for consensus regulatory patterns such as promoter signals, splicing sites and transcription
  factor binding sites

4. Infer protein structure based on previously crystallized proteins

5. Help identify members of gene families

If you work in bioinformatics, chances are that you will need to run some BLAST querues or face the need to process BLAST queries
generated by you or by another person. Biopython proves tools for both taskss.

----

Blast and its variants searches protein and nucleic acid sequences using the BLAST or FASTA method. Both methods 
find similar protein or nucleic aicd chains inthe PDB. psi-blast is used to find more distantly related seuqences

Sequences can be search in two ways

- By PDB ID and Chain ID. Type in a PDB in the structure ID text box and select a chain ID from the pull-down menu. This is useful
  to find all sequences that are similar to the sequence from the specified chain.

- 

"""

def psiBlastScoring(PATH, PSIBLASTPATH = None):

	"""
	Biopython has a wrapper for each BLAST executable, so you can run a blast program from inside your 
	script. The wrapper for blastn 
	
	NcbiblastnCommandline(blast executable, program name, database, input file, ) .. 

	This function returns a tuple with two file objects. The first one is the actual result 
	the second one is the blast error message.

	The output is in XML format. This information can be parsed using the tools learned or with 
	the tools provided by Biopython. There is also a way to avoid dealing with the XML 
	output by forcing NCbiblastncommandline to use plain text as output. This is done by using 
		
	 ---------------------------------------------------------------
	| Links for resources I have been looking at to make this code: |
	 ---------------------------------------------------------------

	-> https://www.rcsb.org/pages/help/advancedsearch/sequence

	-> https://www.biostars.org/p/10419/

	-> https://biopython.org/DIST/docs/api/Bio.PDB.Polypeptide-module.html

	-> https://www.ncbi.nlm.nih.gov/books/NBK2590/ - Good Psiblast explanation
	--------------------------------------------------------------

	JB's instructions on psi-blast
	
	1. We find all similar interfaces.
	
	2. Make a MSA of structurally aligned sequences (multiple sequence alignment) 

	3. Form a score from the probability of a particular mutation showing up in the 
	   MSA 

	4. Evaluate the mutation with this score  - This is where we need the mutation


	-----------------------------------------------    
	|	What to extract from each PSIBLAST record |
	-----------------------------------------------

	Information required to build the blast tables according to  J. B. 
	
	1. The aligned sequences.

	2. the statistics for the alignment <Statistics set>.

	3. The species name for each.

	4. Hsp_evalue.

	Notes:
	------
	-> Should not require an HPC to run each blast computation
	-> The ialign work is run on the biowulf cluster in the NIH, from what I know

	Some biopython options for blast:
	---------------------------------------------------------
	blastn -> nucleotide vs nucleotide
	blastp -> protein vs protein 
	blastx -> translated nucleotide vs protein
	tblastn -> protein vs translcated nucleotide
	tblastx -> translated nucelotide vs translated nucleotide
	---------------------------------------------------------
	Parameters
	----------
	PATH: 
      Path to where the wild type PDBs are found 
	PSIBLASTPATH:
	  Path to where the psiblast binary is 
	"""
	try:
		from Bio.PDB.PDBIO import PDBIO
		from Bio.PDB.PDBParser import PDBParser # PDBparser
		from Bio.Data.IUPACData import protein_letters
		from Bio.SeqUtils.ProtParam import ProteinAnalysis
		from Bio.PDB.Polypeptide import PPBuilder  
		from Bio.PDB.Polypeptide import standard_aa_names # Standard amino acid names - https://biopython.org/DIST/docs/api/Bio.PDB.Polypeptide-module.html 
		from Bio.PDB.Polypeptide import aa1 #  aa1 = 'ACDEFGHIKLMNPQRSTVWY'
		from Bio.PDB.Polypeptide import aa3 #  aa3 = ['ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE',... ]
		from Bio import AlignIO

		# ----------------------------------
		# | Modules for implementing BLAST |
		# ----------------------------------

		from Bio.Blast.Applications import NcbipsiblastCommandline as psiblastn  # psiblast reader
		from Bio.Blast import NCBIXML # For reading the BLAST output 

		# ---------------------------------------------
		# | Boilerplate modules to read the sequences |
		# ---------------------------------------------

		from Bio.Seq import Seq
		from Bio.Seq import translate, transcribe, back_transcribe
		from Bio.PDB.Polypeptide import PPBuilder
		from Bio.Alphabet import IUPAC
		# ----------------
		# | Align module |
		# ----------------
		
		from Bio.Align import MultipleSeqAlignment
		from Bio.SeqRecord import SeqRecord
		from Bio.Alphabet import generic_protein
	except ImportError:
		print ("Error - cannot imoort BLAST python modules")

	# --------------
	# | PDB parser |
	# --------------
	
	# Read the WT pdbs and 
	BLAST_EXE = '/home/oohnohnoh1/Desktop/ACADEMIA/Papermaking/OPTIMUS_BIND/ncbi-blast-2.9.0+/bin/psiblast' # The example given is /home/sb/opt/ncbi-
	# First things - need 
	WTArray = []
	nameArray = []
	for file in os.listdir(PATH): # List the fxout files in the directory, and store them in the array 
		if file.endswith(".pdb"):
			if file[0] == '.':
				pass
			else:
				FileLocation = os.path.join(PATH, file)
				WTArray.append(FileLocation) # Array with the appended path and the pdb file
				parser = PDBParser(PERMISSIVE=1)
				strand_name = file.split('.')
				structure = parser.get_structure(str(strand_name[0]), FileLocation)
				model = structure[0]				
				ppb = PPBuilder()
				seq_rec = []
				subprocess.Popen("mkdir {}_fasta".format(strand_name[0]), shell = True)				
				for index, pp in enumerate(ppb.build_peptides(structure)):
					try:
						sequenceCreator = SeqRecord(Seq(str(pp.get_sequence()), generic_protein), id = str(model.get_list()[index].id))
						align = MultipleSeqAlignment([sequenceCreator])
						AlignIO.write(align,'{}_{}.fasta'.format(strand_name[0], str(model.get_list()[index].id)), 'fasta')
						blast_statistics = ['Statistics_db-num', 'Statistics_db-len', 'Statistics_hsp-len', 'Statistics_eff-space', 'Statistics_kappa', 'Statistics_lambda', 'Statistics_entropy', 'Hsp_evalue', 'Hsp_qseq', 'Hsp_hseq']
						subprocess.Popen("mv {}_{}.fasta {}_fasta/.".format(strand_name[0], str(model.get_list()[index].id), strand_name[0]), shell = True)
						
						psiblastn_cline = psiblastn(cmd = BLAST_EXE, query = '{}_fasta/{}_{}.fasta'.format(strand_name[0], strand_name[0], str(model.get_list()[index].id)), db = "/home/oohnohnoh1/Desktop/ACADEMIA/Papermaking/OPTIMUS_BIND/PANDAS_TABLE/db/cdd_delta", evalue = .0005, outfmt=5, out="{}_fasta/{}_{}.xml".format(strand_name[0], strand_name[0], str(model.get_list()[index].id) ))
						rh,eh = psiblastn_cline()
						tree = ET.parse('fastas/{}_fasta/{}_{}.xml'.format(strand_name[0], strand_name[0], str(model.get_list()[index].id)))
						root = tree.getroot()
						hspNumList = []
						
						for index in root.iter(): # Append all hsp_num
							if index.tag == 'Hit_num':
								hspNumList.append(index.text)
							else:
								pass
								#print (index.tag, index.text)
							#print (index_HIT, index.tag, index.text)
						print (hspNumList)
						hspNumList = []
						for m in root.iter():
							for i in hspNumList:
								if m.text == i:
									print (root.iter('Hsp_qseq').text)									
					except IndexError:
						print ("Error for {}!".format(file))
	return WTArray

"""

"""

#if __name__ == '__main__':
#    log_fmt = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
#    logging.basicConfig(level=logging.INFO, format=log_fmt)#
#   # not used in this stub but often useful for finding various files
#    project_dir = Path(__file__).resolve().parents[2]

   # find .env automagically by walking up directories until it's found, then
   # load up the .env entries as environment variables
#    load_dotenv(find_dotenv())
#    main()


#@click.command()
# fix this patchwork later
#@click.argument('input_filepath', type=click.Path(exists=True))
#@click.argument('output_filepath', type=click.Path())
#def main(): #main(input_filepath, output_filepath):
#    """ Runs data processing scripts to turn raw data from (../raw) into
#        cleaned data ready to be analyzed (saved in ../processed).
#    """
#    input_filepath = 'data/raw/'
#    SKEMPItoPandas('skempi_v2.csv')  # GENERALIZE!

#    logger = logging.getLogger(__name__)
#    logger.info('making final data set from raw data')
