# -*- coding: utf-8 -*-
"""
Created on Sun Oct 13 13:26:33 2019

@author: JRB
"""
# =============================================================================
# These should be set as global variables in an initialization script
# =============================================================================
#FOLDX_loc="/home/jbrender/optimus"
FOLDX_loc="C:/Users/Jeeves/Documents/optimus"
WT_loc="./SKEMPI2_PDBs/PDBs" #must be given as relative path from FOLDX executable
RepairedPDBs_loc="./Soft" #must be given as relative path from FOLDX executable
MutantPDBs_loc="./Soft" #must be given as relative path from FOLDX executable
# =============================================================================
 
import subprocess
import os
import shutil
#import linecache 
import numpy as np

class foldx:
    def __init__(self):
        self.FOLDX_loc=FOLDX_loc
        
# =============================================================================
# Fixes errors in WT structure prior to creating mutation        
# typical usage x.repair_pdb("1CSE.pdb")
# Outfile is then 1CSE_Repaired.pdb in self.out_dir       
# =============================================================================
    def repair_pdb(self,
                   pdb,
                   water=" --water=CRYSTAL",
                   ion=" --ionStrength=0.15",
                   pH=" --pH=7.3",
                   temp=" --temperature=298",
                   vdw=" --vdwDesign=1", #0=soft,1=medium, 2=hard repulsion
                   repair=" --repair_Interface=ONLY"):
        self.out_dir=RepairedPDBs_loc
        self.PDB_loc=WT_loc        
        cmd_str=(self.FOLDX_loc+"/foldx"
                 +" --command=RepairPDB"
                 +" --output-dir="+self.out_dir
                 +water
                 +temp
                 +pH
                 +ion
                 +vdw
                 +" --pdb-dir="+self.PDB_loc + " --pdb="+pdb)
        print(cmd_str)
        try:
            subprocess.run(cmd_str,shell=True)
        except:
            file = open(FOLDX_loc+"/repair_errors.txt","a+")  
            file.write(self.PDB_loc+pdb+"\n")
            
# =============================================================================
# Creates mutation. Optimize WT structure first with optimize_pdb
# Typical usage x.make_mut("Optimized_1CSE_Repair.pdb","LI38G","1CSE_LI38G")           
# =============================================================================
    def make_mut(self,
                   pdb,
                   mutations,
                   out,
                   water=" --water=CRYSTAL",
                   ion=" --ionStrength=0.15",
                   vdw=" --vdwDesign=1", #0=soft,1=medium, 2=hard repulsion                   
                   pH=" --pH=7.3",
                   temp=" --temperature=298"):
        self.out_dir=MutantPDBs_loc
        self.PDB_loc=RepairedPDBs_loc
        outfile=" --output-file="+out             
        if os.path.isfile(FOLDX_loc+"/individual_list.txt"): #This filename is hardcoded into Foldx
            os.remove(FOLDX_loc+"/individual_list.txt")
        file = open(FOLDX_loc+"/individual_list.txt","w+")  
        file.write(mutations+";")
        file.close()
        cmd_str=(self.FOLDX_loc+"/foldx"
                 +" --command=BuildModel"
                 +" --mutant-file="+FOLDX_loc+"/individual_list.txt"
                 +" --output-dir="+self.out_dir
                 +water
                 +temp
                 +pH
                 +ion
                 +vdw
                 +outfile #This doesn't work correctly in Foldx
                 +" --out-pdb=true"
                 +" --pdbHydrogens=true"
                 +" --pdb-dir="+self.PDB_loc
                 + " --pdb="+pdb)
        print(cmd_str)
        try:
            subprocess.run(cmd_str,shell=True)
# =============================================================================
# #Workaround for broken outfile piping            
            pdb_file_str=MutantPDBs_loc+"/"+pdb[0:-4]+"_1.pdb"
            new_pdb_str=MutantPDBs_loc+"/"+out+".pdb"
            shutil.copyfile( pdb_file_str, new_pdb_str)
            os.remove(pdb_file_str)
# =============================================================================
        except:
            file = open(FOLDX_loc+"/errors_mutation.txt","a+")  
            file.write(pdb+"\n")

# =============================================================================
# Small adjustments to structure to minimize energy
# Run for both WT and mutant structures
#typical usage x.optimize_pdb("1CSE_Repair.pdb",RepairedPDBs_loc,RepairedPDBs_loc)
#or x.optimize_pdb("1CSE_LI38G.pdb",MutantPDBs_loc,MutantPDBs_loc)
#output file is Optimized_1CSE_LI38G.pdb in outdir                              
# =============================================================================
    def optimize_pdb(self,
                   pdb,
                   pdbdir,
                   outdir,
                   water=" --water=CRYSTAL",
                   ion=" --ionStrength=0.15",
                   pH=" --pH=7.3",
                   vdw=" --vdwDesign=1", #0=soft,1=medium, 2=hard repulsion                   
                   temp=" --temperature=298"):               
        cmd_str=(self.FOLDX_loc+"/foldx"
                 +" --command=Optimize"
                 +" --output-dir="+outdir
                 +water
                 +temp
                 +pH
                 +ion
                 +vdw
                 +" --out-pdb=true"
                 +" --pdbHydrogens=true"
                 +" --pdb-dir="+pdbdir
                 + " --pdb="+pdb) 
        print(cmd_str)
        try:
            subprocess.run(cmd_str,shell=True)
        except:
            file = open(FOLDX_loc+"/optimize_errors.txt","a+")  
            file.write(pdbdir+pdb+"\n")
            
# =============================================================================
# Calculate energy terms and creates list to be added to dataframe
# Note: the energy terms in the diff file include
# intrachain terms arising from the mutation. Use this instead
# typical usage row=x.char_interface("Optimized_1CSE_LI38G.pdb",MutantPDBs_loc,MutantPDBs_loc,"E","I")
# typical usage row=x.char_interface("Optimized_1CSE_Repair.pdb",,RepairedPDBs_loc,RepairedPDBs_loc,"E","I")           
# =============================================================================
    def char_interface(self,
                   pdb,
                   pdbdir,
                   outdir,
                   chainA,
                   chainB,
                   water=" --water=CRYSTAL",
                   ion=" --ionStrength=0.15",
                   pH=" --pH=7.3",
                   vdw=" --vdwDesign=1", #0=soft,1=medium, 2=hard repulsion                   
                   temp=" --temperature=298",
                   key="999",
                   mutation="NA"):             
        cmd_str=(self.FOLDX_loc+"/foldx"
                 +"  --command=AnalyseComplex"
                 +" --output-dir="+outdir
                 +water
                 +temp
                 +pH
                 +ion
                 +vdw
                 +" --pdb-dir="+pdbdir
                 + " --pdb="+pdb) 
        print(cmd_str)
        try:
            subprocess.run(cmd_str,shell=True)
            name=outdir+"/Interaction_"+pdb[0:-4]+"_AC.fxout" #not sure of meaning of AC
            #nrg=linecache.getline(name, 10)
            with open(name) as file:
                fileAsList = file.readlines()
                energy=[0]*21
                for i in range(9, len(fileAsList)):
                    nrg=fileAsList[i]
                    nrgFields=nrg.split('\t')
                    if ((chainA.find(nrgFields[1])!=-1) and (chainB.find(nrgFields[2])!=-1)) or (
                        (chainB.find(nrgFields[1])!=-1) and (chainA.find(nrgFields[2])!=-1)):
                        energy=np.add(list(map(float, nrgFields[5:26])), energy)                       
                rowValue = {
                        'Key':key,
                        'Mutation':mutation,                        
                        'TotalFoldX':energy[0],
                        'Backbone_HBond':nrgFields[1],
                        'Sidechain_Hbond':nrgFields[2],
                        'Van_der_Waals':nrgFields[3],
                        'Electrostatics':nrgFields[4],
                        'Solvation_Polar':nrgFields[5],
                        'Solvation_Hydrophobic':nrgFields[6],
                        'clashes':nrgFields[7],
                        'entropy_sidechain':nrgFields[8],
                        'entropy_mainchain':nrgFields[9],
                        'sloop_entropy':nrgFields[10],
                        'mloop_entropy':nrgFields[11],
                        'cis_bond':nrgFields[12],
                        'torsional_clash':nrgFields[13],
                        'backbone_clash':nrgFields[14],
                        'helix_dipole':nrgFields[15],
                        'water_bridge':nrgFields[16],
                        'disulfide':nrgFields[17],
                        'kon':nrgFields[18],
                        'partial_covalent':nrgFields[19],
                        'Ion':nrgFields[20],
                        'Entropy Complex':nrgFields[21]}
            return rowValue
        except:
            file = open(FOLDX_loc+"/char_errors.txt","a+")  
            file.write(pdbdir+pdb+"\n")
            