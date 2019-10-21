# -*- coding: utf-8 -*-
"""
Created on Sun Oct 13 13:26:33 2019

@author: JRB
"""
# =============================================================================
# These should be set as global variables in an initialization script
# =============================================================================
FOLDX_loc="/home/jbrender/optimus"
WT_loc="./SKEMPI2_PDBs/PDBs" #must be given as relative path from FOLDX executable
RepairedPDBs_loc="./Medium" #must be given as relative path from FOLDX executable
MutantPDBs_loc="./Medium" #must be given as relative path from FOLDX executable
csv_loc="./Medium/featurelist" 

vdw="1"
# =============================================================================

import subprocess
import os
import shutil
import pandas as pd
from timeit import default_timer as timer
from runfoldx import foldx
from SKEMPI import skempi_final as db
from pathlib import Path

os.chdir(FOLDX_loc)
fx=foldx()

def writeRow(fileloc,headerFlag,row):
    with open(fileloc, 'a+') as outfile:
        if headerFlag==True:
            for key, value in row.items():
                outfile.write(key+",")
            outfile.write("\n")
            headerFlag=False
            for key, value in row.items():
                outfile.write(value+",")
            outfile.write("\n")
        else:
            for key, value in row.items():
                outfile.write(value+",")
            outfile.write("\n")
    return (headerFlag)            
# =============================================================================
# Creates optimized structures of the WT PFB and calculates the foldx energy
# features for the WT
# =============================================================================
def WT(db):
    pdbs= db['#Pdb'].str.slice(start=0, stop=4, step=1)
    chainA=db['#Pdb'].str.slice(start=5, stop=6, step=1)
    chainB=db['#Pdb'].str.slice(start=7, stop=8, step=1)
    uniq_pdbs=pdbs.unique()
    i=0
    WTheaderFlag=True    
    WTrows_list=[]
    dirname = Path(RepairedPDBs_loc)
    if not (dirname.is_dir()):
        os.mkdir(RepairedPDBs_loc)            
    while i<len(uniq_pdbs):
        WTpdb=uniq_pdbs[i]+".pdb"
        chA=chainA[i]
        chB=chainB[i]
        print(WTpdb+chA+chB)
        if not os.path.exists(WT_loc+"/"+WTpdb):
            with open(FOLDX_loc+"/missing_WT_pdb.txt","a+") as missingWTfile: 
                missingWTfile.write(WTpdb+"\n")
        else:
            start = timer()
            fx.repair_pdb(WTpdb)
            end = timer()
            repair_time=end-start;
            opt_loc=WTpdb[0:-4]+"_Repair.pdb"
            fx.optimize_pdb(opt_loc,RepairedPDBs_loc,RepairedPDBs_loc)
            WTrow=fx.char_interface(opt_loc,RepairedPDBs_loc,RepairedPDBs_loc,chA,chB)
            with open(FOLDX_loc+"/FoldX_WTTimes.txt","a+") as wt_times:
                try:
                    wt_times.write(str(uniq_pdbs[i])+"\t"+str(repair_time)+"\n")
                except:
                    pass
            WTrow.update({'pdb':WTpdb[0:-4]})
            WTrows_list.append(WTrow)
            #writing to csv file allows data to be kept even if program crashes            
            WTheaderFlag=writeRow(csv_loc+"/WTFoldXNrgs.csv",WTheaderFlag,WTrow) 
        i+=1
    WTdf = pd.DataFrame(WTrows_list)
    return(WTdf)
    
# ============================================================================= 
# Creates optimized structures of each mutant and calculates the foldx energy
# features for each mutant
# =============================================================================   
def mut(db):
    mutations=db['Mutation(s)_cleaned']
    pdb=db['#Pdb'].str.slice(start=0, stop=4, step=1)
    chainA=db['#Pdb'].str.slice(start=5, stop=6, step=1)
    chainB=db['#Pdb'].str.slice(start=7, stop=8, step=1)
    i=0
    MutHeaderFlag=True
    MutRows_list=[]
    dirname = Path(MutantPDBs_loc)
    if not (dirname.is_dir()):
        os.mkdir(MutantPDBs_loc) 
    while i<len(mutations):
        if os.path.isfile(FOLDX_loc+"/individual_list.txt"): #The name of this file is hardcoded in FoldX
            os.remove(FOLDX_loc+"/individual_list.txt")
        file = open(FOLDX_loc+"/individual_list.txt","a+")  
        file.write(mutations.iloc[i]+";")
        file.close() 
        Temp=" --temperature="+str(db.Temperature.iloc[i])
        vdwRep=" --vdwDesign="+vdw
        mutpdb="Optimized_"+pdb.iloc[i]+"_Repair.pdb"
		#Files are keyed to database index
        outfile=str(db.index[i])+"_"+pdb.iloc[i]+"_"+mutations.iloc[i] 
        outpdb=outfile+".pdb"
		
		#Make mutations from optimized repaired WT structure		
        start = timer()
        fx.make_mut(mutpdb,mutations.iloc[i],outfile,temp=Temp,vdw=vdwRep)
        end = timer()
        mut_time=end-start;
		
		#Optimize the repair		
        start = timer()
        fx.optimize_pdb(outpdb,MutantPDBs_loc,MutantPDBs_loc,temp=Temp,vdw=vdwRep)
        end = timer()
        opt_time=end-start;
		
		#Report the energy terms		
        start = timer()
        MutRow=fx.char_interface(outpdb,MutantPDBs_loc,MutantPDBs_loc,chainA,chainB,temp=Temp,vdw=vdwRep)
        end = timer()
        char_time=end-start;
		
        mut_file = open(FOLDX_loc+"/FoldX_MutTimes.txt","a+")
        try:
            mut_file.write(str(db.index[i])+"\t"+pdb.iloc[i]+"\t"+mutations.iloc[i]+"\t"+str(mut_time)+"/t"+str(opt_time)+"\t"+str(char_time)+"\n")
        except:
            pass
        MutRow.update({'pdb':mutpdb[0:-4]})
        MutHeaderFlag=writeRow(csv_loc+"/MutFoldXNrgs.csv",MutHeaderFlag,MutRow)         
        MutRows_list.append(MutRow)        
        i+=1  
        
    Mutdf = pd.DataFrame(MutRows_list)
    return(Mutdf)
 # =============================================================================
 
wtdb=WT(db)
mutdb=mut(db)
wtdb.to_csv(csv_loc+"WTFoldX.csv")     
mutdb.to_csv(csv_loc+"MutFoldX.csv")           