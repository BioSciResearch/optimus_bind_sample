# -*- coding: utf-8 -*-
"""
Created on Sun Oct 13 13:26:33 2019

@author: JRB
"""

# =============================================================================
# PURPOSE: Loops through the SKEMPI database and runs FOLDX for each mutation along
# with the corresponding wild type. Skips PDBs too large to it in memory (set by maxfile)
#
# Generates at ./(Soft Medium or Hard)/featurelist :
# 1 MutantFoldXNrgs(_date).csv Wild type energy scores
# 2 WTFoldXNrgs(_date).csv Wild type energy scores
# 3 FoldX_MutTimes(_date).txt Calculation time per mutant, all inclusive
# 4 FoldX_WTTimes(_date).txt Calculation time per WT protein, all inclusive
# 5 missing_WT_pdb(_date).txt Skipped or missing WT pdbs
# 6 missing_Mut_pdb(_date).txt Missing mutant pdbs
# =============================================================================
# DEPENDENCIES
# runfoldx.py
# SKEMPI.py
# skempi_v2.csv in run directory
#
# =============================================================================
#
# INPUT PARAMETERS
# These should be set as global variables in an initialization script
# =============================================================================
# Changeable Parameters
# =============================================================================
# Sets the van der Waal repulsion parameter 0=soft 1=medium 2=hard
vdw="0"
# Sets the maximum filesize of the pdb. Larger than this the program notes the pdb and skips
maxfile=750000
#
# =============================================================================
# Mostly fixed Parmeters
# =============================================================================
switcher = {
	"0": "Soft",
	"1": "Medium",
	"2": "Hard",
}
vdwloc=switcher.get(vdw,"Invalid vdw parameter")
FOLDX_loc="/home/jbrender/optimus"
WT_loc="./SKEMPI2_PDBs/PDBs" #must be given as relative path from FOLDX executable
RepairedPDBs_loc="./"+vdwloc+"/Repaired" #must be given as relative path from FOLDX executable
MutantPDBs_loc="./"+vdwloc+"/Mutant" #must be given as relative path from FOLDX executable
csv_loc="./"+vdwloc+"/featurelist" 
# =============================================================================



import os
import datetime
import pandas as pd
from timeit import default_timer as timer
from runfoldx import foldx
from SKEMPI import skempi_final as db
from pathlib import Path

os.chdir(FOLDX_loc)
fx=foldx(vdw)
date= datetime.datetime.now()
datestr='_'+str(date.month)+'_'+str(date.day)+'_'+str(date.year)

# =============================================================================
# Begin function definition
# =============================================================================
def writeRow(fileloc,headerFlag,row):
    with open(fileloc, 'a+') as outfile:
        if headerFlag==True:
            for rowkey, value in row.items():
                outfile.write(rowkey+";")
            outfile.write("\n")
            headerFlag=False
            for rowkey, value in row.items():
                outfile.write(str(value)+";")
            outfile.write("\n")
        else:
            for rowkey, value in row.items():
                outfile.write(str(value)+";")
            outfile.write("\n")
    return (headerFlag)            
# =============================================================================
# Creates optimized structures of the WT PFB and calculates the foldx energy
# features for the WT
# =============================================================================
def WT(db):
    pdbs=db[['PDB', 'chainA','chainB']]
    uniq_pdbs=pdbs.drop_duplicates()
    i=0
    WTheaderFlag=True    
    WTrows_list=[]
    dirname = Path(RepairedPDBs_loc)
    if not (dirname.is_dir()):
        os.mkdir(RepairedPDBs_loc)            
    while i<len(uniq_pdbs):
        WTpdb=uniq_pdbs['PDB'].iloc[i]+".pdb"
        #patch to skip overly large pdb files due to linode's memory restrictions
        if os.path.getsize(WT_loc+"/"+WTpdb)>maxfile:
            with open(FOLDX_loc+"/missing_WT_pdb"+datestr+".txt","a+") as missingWTfile: 
                missingWTfile.write(WTpdb+"\n")
            i+=1
            continue
        chA=uniq_pdbs['chainA'].iloc[i]
        chB=uniq_pdbs['chainB'].iloc[i]
        print(WTpdb+chA+chB)
        if not os.path.exists(WT_loc+"/"+WTpdb):
            with open(FOLDX_loc+"/missing_WT_pdb"+datestr+".txt","a+") as missingWTfile: 
                missingWTfile.write(WTpdb+"\n")
        else:
            if os.path.isfile(RepairedPDBs_loc+"/"+WTpdb[0:-4]+"_Repair.pdb") is False:	
                start = timer()
                fx.repair_pdb(WTpdb)
                end = timer()
                repair_time=end-start;
            opt_loc=WTpdb[0:-4]+"_Repair.pdb"
            fx.optimize_pdb(opt_loc,RepairedPDBs_loc,RepairedPDBs_loc)
            WTrow=fx.char_interface(opt_loc,RepairedPDBs_loc,RepairedPDBs_loc,chA,chB)
            with open(FOLDX_loc+"/FoldX_WTTimes"+datestr+".txt","a+") as wt_times:
                try:
                    wt_times.write(str(uniq_pdbs[i])+"\t"+str(repair_time)+"\n")
                except:
                    pass
            WTrow.update({'pdb':WTpdb[0:-4]})
            WTrows_list.append(WTrow)
            #writing to csv file allows data to be kept even if program crashes            
            WTheaderFlag=writeRow(csv_loc+"/WTFoldXNrgs"+datestr+".csv",WTheaderFlag,WTrow) 
        i+=1
    WTdf = pd.DataFrame(WTrows_list)
    return(WTdf)
    
# ============================================================================= 
# Creates optimized structures of each mutant and calculates the foldx energy
# features for each mutant
# =============================================================================   
def mut(db):
    mutations=db['Mutation(s)_cleaned']
    pdb=db['PDB']
    chainA=db['chainA']
    chainB=db['chainB']
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
        outfile=str(db.Key.iloc[i])+"_"+pdb.iloc[i]+"_"+mutations.iloc[i] 
        outpdb=outfile+".pdb"
		#Make mutations from optimized repaired WT structure if not done already
        if (os.path.isfile(MutantPDBs_loc+"/"+pdb.iloc[i]+"_"+mutations.iloc[i]+"_Repair.pdb") is False) and (
			os.path.isfile(RepairedPDBs_loc+"/"+pdb.iloc[i]+"_Repair.pdb") is True):	
            start = timer()
            fx.make_mut(mutpdb,mutations.iloc[i],outfile,temp=Temp,vdw=vdwRep)
            end = timer()
            mut_time=end-start;
            mut_time_file = open(FOLDX_loc+"/FoldX_MutTimes"+datestr+".txt","a+")
            try:
                mut_time_file.write(str(db.index[i])+"\t"+pdb.iloc[i]+"\t"+mutations.iloc[i]+"\t"+str(mut_time)+"\n")
            except:
                pass
		#If Repaired WT PDB file cannot be made skip mutations and record mutation as missing
        elif os.path.isfile(RepairedPDBs_loc+"/"+pdb.iloc[i]+"_Repair.pdb") is False:
            with open(FOLDX_loc+"/missing_Mut_pdb"+datestr+".txt","a+") as missingMutfile: 
                missingMutfile.write(db.Key.iloc[i]+"\t"+pdb.iloc[i]+"\t"+mutations.iloc[i]+"\n")
            i+=1
            continue

		#Optimize the repair if not already
        if os.path.isfile(MutantPDBs_loc+"/Optimized"+pdb.iloc[i]+"_"+mutations.iloc[i] +"_Repair.pdb") is False:			
            start = timer()
            fx.optimize_pdb(outpdb,MutantPDBs_loc,MutantPDBs_loc,temp=Temp,vdw=vdwRep)
            end = timer()
            opt_time=end-start;
            mut_opt_time_file = open(FOLDX_loc+"/FoldX_MutOptTimes"+datestr+".txt","a+")
            try:
                mut_opt_time_file.write(str(db.index[i])+"\t"+pdb.iloc[i]+"\t"+mutations.iloc[i]+"\t"+str(opt_time)+"\n")
            except:
                pass		
		
        #Report the energy terms			
        start = timer()
        MutRow=fx.char_interface(outpdb,MutantPDBs_loc,MutantPDBs_loc,chainA.iloc[i],chainB.iloc[i],temp=Temp,vdw=vdwRep,key=db.Key.iloc[i], mutation=mutations.iloc[i])
        end = timer()
        char_time=end-start;
        mut_char_time_file = open(FOLDX_loc+"/FoldX_MutCharTimes"+datestr+".txt","a+")
        try:
            mut_char_time_file.write(str(db.index[i])+"\t"+pdb.iloc[i]+"\t"+mutations.iloc[i]+"\t"+str(char_time)+"\n")
        except:
            pass
        MutRow.update({'pdb':mutpdb[0:-4]})
        writeRow(csv_loc+"/MutFoldXNrgs"+datestr+".csv",MutHeaderFlag,MutRow)
        MutHeaderFlag=True         
        MutRows_list.append(MutRow)        
        i+=1  
        
    Mutdf = pd.DataFrame(MutRows_list)
    return(Mutdf)
	
# =============================================================================
# Main program driver
# =============================================================================
#wtdb=WT(db)
mutdb=mut(db)
##wtdb.to_csv(csv_loc+"WTFoldX.csv",sep=';')     
mutdb.to_csv(csv_loc+"MutFoldX.csv",sep=';')           
=======
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
