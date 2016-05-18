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


def featureVina(inprot, inlig):
    """Get Vina score and Vina features
    
    The first returned value is vina score not features
    
    Parameters
    ----------
    inprot : str
        Input file name of protein
    inlig : str
        Input file name of ligand
    num : int (default 10)
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
    
    # convert vina score and interaction term to be pKd unit
    c = -0.73349
    vinalist[0:53] = [c * i for i in vinalist[0:53]]
            
    return vinalist
    

class vina:
    """Vina score and vina features
    
    """
    
    def __init__(self, prot, lig):
        """Vina Socre and Vina Features
        
        Parameters
        ----------
        prot : str
            protein structure
        lig : str
            ligand structure
        
        """
        self.prot = prot
        self.lig = lig
        
        vinalist = featureVina(self.prot, self.lig)     
        
        self.vinaScore = vinalist[0]
        self.vinaFeatures = vinalist[1:]


    def features(self, num = 10):
        """Get subset of features
        
        Parameters
        ----------
        num : int (default 10)
            number of features to retrieve
        
        """
        idx6 = [10, 28, 37, 44, 49, 52]
        idx10 = [0, 2, 52, 54, 53, 55, 3, 51, 57, 47]
        
        if num == 10:
            return [self.vinaFeatures[i] for i in idx10]
        elif num == 6:
            return [self.vinaFeatures[i] for i in idx6]
        elif num == 58:
            return self.vinaFeatures
            
    def info(self, num = 10):
        """
        """
        
        vinafeatname = ('ad4_solvation(charge=T)', 
                        'ad4_solvation(charge=F)',
                        'electrostatic(x=1)', 
                        'electrostatic(x=2)', 
                        'gauss(0,0.3)', 
                        'gauss(0.5,0.3)', 
                        'gauss(1,0.3)', 
                        'gauss(1.5,0.3)', 
                        'gauss(2,0.3)', 
                        'gauss(2.5,0.3)', 
                        'gauss(0,0.5)', 
                        'gauss(1,0.5)', 
                        'gauss(2,0.5)', 
                        'gauss(0,0.7)', 
                        'gauss(1,0.7)', 
                        'gauss(2,0.7)', 
                        'gauss(0,0.9)', 
                        'gauss(1,0.9)', 
                        'gauss(2,0.9)', 
                        'gauss(3,0.9)', 
                        'gauss(0,1.5)', 
                        'gauss(1,1.5)', 
                        'gauss(2,1.5)', 
                        'gauss(3,1.5)', 
                        'gauss(4,1.5)', 
                        'gauss(0,2)', 
                        'gauss(1,2)', 
                        'gauss(2,2)', 
                        'gauss(3,2)', 
                        'gauss(4,2)', 
                        'gauss(0,3)', 
                        'gauss(1,3)', 
                        'gauss(2,3)', 
                        'gauss(3,3)', 
                        'gauss(4,3)', 
                        'repulsion(0.4)', 
                        'repulsion(0.2)', 
                        'repulsion(0.0)', 
                        'repulsion(-0.2)', 
                        'repulsion(-0.4)',
                        'repulsion(-0.6)', 
                        'repulsion(-0.8)', 
                        'repulsion(-1.0)', 
                        'hydrophobic(0.5,1)', 
                        'hydrophobic(0.5,1.5)', 
                        'hydrophobic(0.5,2)', 
                        'hydrophobic(0.5,3)', 
                        'non_hydrophobic(0.5,1.5)', 
                        'vdw(4,8)', 
                        'non_dir_h_bond(-0.7,0)', 
                        'non_dir_h_bond(-0.7,0.2)', 
                        'non_dir_h_bond(-0.7,0.4)', 
                        'num_tors', 
                        'num_rotors', 
                        'num_heavy_atoms', 
                        'num_hydrophobic_atoms', 
                        'ligand_max_num_h_bonds', 
                        'ligand_length')
        
        
        idx6 = [10, 28, 37, 44, 49, 52]
        idx10 = [0, 2, 52, 54, 53, 55, 3, 51, 57, 47]
        
        if num == 10:
            for i,item in enumerate(idx10):
                print i+1, vinafeatname[item]
        elif num == 6:
            for i,item in enumerate(idx6):
                print i+1, vinafeatname[item]
        else:
            for i in range(58):
                print i+1, vinafeatname[i]  


if __name__ == "__main__":
    
    testdir = 'tests/1a42/'
    inprot = testdir + '1a42_protein_proc_se.pdb'
    inlig = testdir + '1a42_ligand_fix.mol2'

    v =  vina(inprot, inlig)
    print v.vinaScore
    print v.features(6)
    print v.features(10)
    print v.features(58)
    v.info(10)
    v.info(6)
    v.info(58)
    