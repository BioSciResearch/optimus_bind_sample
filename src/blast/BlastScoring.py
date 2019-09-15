"""
----------------
| OPTIMUS BIND | 
----------------

------------------------------------------------------------------------------------------------

Version: 0.0.1

Last Updated: 11/09/2019

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

	BLASTEXE = '/home/oohnohnoh1/Desktop/ACADEMIA/Papermaking/OPTIMUS_BIND/ncbi-blast-2.9.0+/bin/psiblast' # The example given is /home/sb/opt/ncbi-

	# First things - need 

	WTArray = []

	nameArray = []
	TotalDict = {}
	for file in os.listdir(PATH): # List the fxout files in the directory, and store them in the array 
		if file.endswith(".pdb"):
			if file[0] == '.':
				pass
			else:
				FileLocation = os.path.join(PATH, file)
				WTArray.append(FileLocation) # Array with the appended path and the pdb file
				Parser = PDBParser(PERMISSIVE=1)
				strandName = file.split('.')
				structure = Parser.get_structure(str(strandName[0]), FileLocation)
				model = structure[0] # PDB loader 
				ppb = PPBuilder()  # PDB builder 
				seqRec = []
				hspNumList = [] 
				subprocess.Popen("mkdir {}_fasta".format(strandName[0]), shell = True)				
				for index, pp in enumerate(ppb.build_peptides(structure)):
					try:
						sequenceCreator = SeqRecord(Seq(str(pp.get_sequence()), generic_protein), id = str(model.get_list()[index].id))
						align = MultipleSeqAlignment([sequenceCreator])
						NAME = "{}_{}".format(strandName[0], str(model.get_list()[index].id))
						AlignIO.write(align,'{}_{}.fasta'.format(strandName[0], str(model.get_list()[index].id)), 'fasta')
						blast_statistics = ['Hsp_evalue', 'Hsp_qseq', 'Hsp_hseq','Statistics_db-num', 'Statistics_db-len', 'Statistics_hsp-len', 'Statistics_eff-space', 'Statistics_kappa', 'Statistics_lambda', 'Statistics_entropy']
						subprocess.Popen("mv {}_{}.fasta {}_fasta/.".format(strandName[0], str(model.get_list()[index].id), strandName[0]), shell = True)
						
						psiblastnCline = psiblastn(cmd = BLASTEXE, query = '{}_fasta/{}_{}.fasta'.format(strandName[0], strandName[0], str(model.get_list()[index].id)), db = "/home/oohnohnoh1/Desktop/ACADEMIA/Papermaking/OPTIMUS_BIND/PANDAS_TABLE/db/cdd_delta", evalue = .0005, outfmt=5, out="{}_fasta/{}_{}.xml".format(strandName[0], strandName[0], str(model.get_list()[index].id) ))
						rh,eh = psiblastnCline()
						tree = ET.parse('fastas/{}_fasta/{}_{}.xml'.format(strandName[0], strandName[0], str(model.get_list()[index].id)))
						root = tree.getroot()
						A = 0  # Default index 
						TotalDict[NAME] = []
						for ind in root.iter(): # Append all hsp_num
							if ind.tag == 'Hit_num':
								hspNumList.append(ind.text)
								A = ind.text
							TotalDict[NAME].append([strandName[0], A, str(model.get_list()[index].id), ind.tag, ind.text])
						numList = [int(row[1]) for row in TotalDict[NAME]]
						numList = list(set(numList))
						numList.sort()
						if len(numList) == 0:
							print ("SINGLE")
							D_1 = [row[4] for row in TotalDict[NAME] if row[3] in blast_statistics]
							print (D_1)
						else:
							print ("MULTIPLE")
							for i in numList:
								if i == 0:
									pass
								else:
									D_1 = [row[4] for row in TotalDict[NAME] if row[3] in blast_statistics]
									print (D_1)
					except IndexError:
						print ("Error for {}!".format(file))
	subprocess.Popen("rm -r *_fasta", shell = True)
	return TotalDict

