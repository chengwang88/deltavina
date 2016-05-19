"""
DeltaVinaRF20 Scoring Function
"""

__author__ = "Cheng Wang"
__copyright__ = "Copyright 2016, NYU"
__license__ = ""

#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------

import os, sys, time
sys.path.insert(0, '../../')

from deltavina.features import featureSASA
from deltavina.features import featureVina

#-----------------------------------------------------------------------------
# Code
#-----------------------------------------------------------------------------

RDIR = ""


def runDVRF20(prot, lig):
    """
    
    """
    
    vina = featureVina.vina(prot, lig)
    sasa = featureSASA.sasa(prot, lig)
    vinafeat = [str(i) for i in [vina.vinaScore] + vina.features(10)]
    sasafeat = [str(i) for i in sasa.sasalist]
    feat = vinafeat + sasafeat

    header = "pdb,vina," + ",".join(['F' + str(i+1) for i in range(20)]) + "\n"

    f = open('input.csv', 'w')
    f.write(header)
    f.write('idx,' + ','.join(feat) + '\n')
    f.close()
    
    f = open('DVRF20.R', 'w')
    f.write("""library(caret)
library(randomForest)

# load the model
load('rffit.rda')

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
    
""")
    f.close()
    
    cmd = "R CMD BATCH DVRF20.R"
    os.system(cmd)
    
if __name__ == "__main__":
    
    testdir = 'tests/1a42/'
    inprot = testdir + '1a42_protein_proc_se.pdb'
    inlig = testdir + '1a42_ligand_fix.mol2'

    runDVRF20(inprot, inlig)
