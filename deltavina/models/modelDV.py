"""
DeltaVinaRF20 Scoring Function
"""

__author__ = "Cheng Wang"
__copyright__ = "Copyright 2016, NYU"
__license__ = ""

#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------

import os, sys

import deltavina
from deltavina.features import featureSASA
from deltavina.features import featureVina


#-----------------------------------------------------------------------------
# Code
#-----------------------------------------------------------------------------

def runDV(pdblist):
    """Give a list of pdblist
    Calculate features for each complex and write features to input.csv
    Run R script to calculate predicted pKd to output.csv

    Parameters
    ----------
    pdblist : list[list[str,str]]
        list of pdb or mol2 of protein and ligand
    """
    
    # feature calculation for each complex
    featlist = []
    for pdb in pdblist:
        prot = pdb[0]
        lig = pdb[1]
        vina = featureVina.vina(prot, lig)
        sasa = featureSASA.sasa(prot, lig)
        vinafeat = [str(i) for i in [vina.vinaScore] + vina.features(10)]
        sasafeat = [str(i) for i in sasa.sasalist]
        feat = vinafeat + sasafeat
        featlist.append(feat)
    
    # write feature to input.csv as input for R script
    header = "pdb,vina," + ",".join(['F' + str(i+1) for i in range(20)]) + "\n"
    f = open('input.csv', 'w')
    f.write(header)
    for idx, feat in enumerate(featlist):
        f.write('pdb' + str(idx) + ',' + ','.join(feat) + '\n')
    f.close()
    
    # path of deltavina module
    rfpath = os.path.dirname(deltavina.__file__)

    # write input R script
    f = open('DVRF20.R', 'w')
    f.write("""library(randomForest)

# load the model
load('%s/models/rffit.rda')

# input and output file name
infn = 'input.csv'
outfn = 'output.csv'

print(paste("Read input: ", infn))
# read in input as dataframe df
df = read.table(infn, header=T, stringsAsFactors = F, sep=',')

print(df)

# get features from df
feats = df[3:22]

# predict the binding affinity
pred = round(predict(rffit, newdata = feats) + df$vina,2)

# write output
output = data.frame(pdb = df$pdb, pred = pred)

print(paste("Write input: ", outfn))
write.table(output, outfn, sep=',', row.names = F, quote = F)

print("Done")
    
""" %(rfpath))
    f.close()
    
    # run R script
    cmd = "R CMD BATCH DVRF20.R"
    os.system(cmd)
    



if __name__ == "__main__":
    """
    Test of deltavina method
    """
    testdir = 'tests/1a42/'
    inprot = testdir + '1a42_protein_proc_se.pdb'
    inlig = testdir + '1a42_ligand_fix.mol2'
    pdblist = [[inprot, inlig]]
    runDV(pdblist)
