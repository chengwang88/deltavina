"""
Vina Features with 6, 10, and 58 Features option
"""

__author__ = "Cheng Wang"
__copyright__ = "Copyright 2016, NYU"
__license__ = ""

#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------

import os, sys
import numpy as np

#-----------------------------------------------------------------------------
# Code
#-----------------------------------------------------------------------------

MGLPY = "/Users/chengwang/sw/mgltools_i86Darwin9_1.5.6/bin/python"
MGLUTIL = "/Users/chengwang/sw/mgltools_i86Darwin9_1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/"
VINADIR = "/Users/chengwang/Dropbox/ws/scoref/vina/v1/build/mac/release/"

def runVina(protpdbqt, ligpdbqt):
    """Run modified AutoDock Vina program with Vina score and 58 features
    
    Parameters
    ----------
    protpdbqt : str
        PDBQT file name of protein
    ligpdbqt : str
        PDBQT file name of ligand
        
    Returns
    ----------
    vinalist : list[float]
        Vina score and 58 features by Vina
        
    """
    cmd = VINADIR + "vina --receptor " + protpdbqt + " --ligand " + ligpdbqt + \
          " --score_only > score_v1.tmp"
    os.system(cmd)
    
    vinalist = []
    with open("score_v1.tmp", "r") as f:
        for lines in f:
            if lines[0:4] in ["Affi", "Term"]:
                vinalist.append(float(lines.split()[2]))
    
    return vinalist
    

def prepareProt(inprot, protpqdbt):
    """Prepare protein PDBQT file by MGLTools
   
    Parameters
    ----------
    inprot : str
        Input file name of protein
    protpqdbt : str
        Output PDBQT file name of protein
    
    """
    cmd = MGLPY + " "  + MGLUTIL + "prepare_receptor4.py -r "  + inprot + \
          " -o " + protpqdbt + " -U 'nphs' > out1.tmp"
    os.system(cmd)



def prepareLig(inlig, ligpdbqt):
    """Prepare ligand PDBQT file by MGLTools
    
    Parameters
    ----------
    inlig : str
        Input file name of ligand
    ligpqdbt : str
        Output PDBQT file name of ligand
        
    """
    cmd = MGLPY + " "  + MGLUTIL + "prepare_ligand4.py -l " + inlig  + \
          " -o " + ligpdbqt +  " -U 'nphs' > out2.tmp"
    os.system(cmd)


def featureVina(inprot, inlig, num = 10):
    """Get Vina score and Vina features
    
    Parameters
    ----------
    inprot : str
        Input file name of protein
    inlig : str
        Input file name of ligand
    num : int
        Number of features to output 
       
    """
    fprot, __ = os.path.splitext(inprot)
    flig, __ = os.path.splitext(inlig)

    protpdbqt = fprot + ".pdbqt"
    ligpdbqt = flig + ".pdbqt"
    
    prepareProt(inprot, protpdbqt)
    prepareLig(inlig, ligpdbqt)
    os.system('rm *.tmp')
    vinalist = runVina(protpdbqt, ligpdbqt)
    print len(vinalist)
    return vinalist
    
    
if __name__ == "__main__":
    
    testdir = 'tests/1a42/'
    inprot = testdir + '1a42_protein_proc_se.pdb'
    inlig = testdir + '1a42_ligand_fix.mol2'

    print featureVina(inprot, inlig)