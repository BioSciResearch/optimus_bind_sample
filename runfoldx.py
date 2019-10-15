# -*- coding: utf-8 -*-
"""
Created on Sun Oct 13 13:26:33 2019

@author: JRB
"""
# =============================================================================
# These should be set as global variables in an initialization script
# =============================================================================
FOLDX_loc="C:/Users/Jeeves/Documents/optimus"
WT_loc="./SKEMPI2_PDBs/PDBs" #must be given as relative path from FOLDX executable
RepairedPDBs_loc="./Repaired2" #must be given as relative path from FOLDX executable
MutantPDBs_loc="./Mutant_Structures2" #must be given as relative path from FOLDX executable
# =============================================================================
 
import subprocess
import os
import shutil
import linecache 

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
                 +"--output-dir="+self.out_dir
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
            file = open(FOLDX_loc+"/repair_errors.txt","w+")  
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
            file = open(FOLDX_loc+"/errors_mutation.txt","w+")  
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
            file = open(FOLDX_loc+"/optimize_errors.txt","w+")  
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
                   temp=" --temperature=298"):             
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
            nrg=linecache.getline(name, 10)
            nrgFields=nrg.split('\t')            
            rowValue = {#'ClashA:':nrgFields[3],
                        #'ClashB:':nrgFields[4],
                        'TotalFoldX':nrgFields[5],
                        'Backbone_HBond':nrgFields[6],
                        'Sidechain_Hbond':nrgFields[7],
                        'Van_der_Waals':nrgFields[8],
                        'Electrostatics':nrgFields[9],
                        'Solvation_Polar':nrgFields[10],
                        'Solvation_Hydrophobic':nrgFields[11],
                        'clashes':nrgFields[12],
                        'entropy_sidechain':nrgFields[13],
                        'entropy_mainchain':nrgFields[14],
                        'sloop_entropy':nrgFields[15],
                        'mloop_entropy':nrgFields[16],
                        'cis_bond':nrgFields[17],
                        'torsional_clash':nrgFields[18],
                        'backbone_clash':nrgFields[19],
                        'helix_dipole':nrgFields[20],
                        'water_bridge':nrgFields[21],
                        'disulfide':nrgFields[22],
                        'kon':nrgFields[23],
                        'partial_covalent':nrgFields[24],
                        'Ion':nrgFields[25],
                        'Entropy Complex':nrgFields[26]}
            return rowValue
        except:
            file = open(FOLDX_loc+"/char_errors.txt","w+")  
            file.write(pdbdir+pdb+"\n")
