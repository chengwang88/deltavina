"""
Script to Run DeltaVinaRF20
"""

__author__ = "Cheng Wang"
__copyright__ = "Copyright 2016, NYU"
__license__ = ""

#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------

import os, sys
from optparse import OptionParser
from deltavina.models import modelDV

#-----------------------------------------------------------------------------
# Code
#-----------------------------------------------------------------------------

# options for read files
parser = OptionParser()
parser.add_option("-r", "--receptor",  help="receptor file name")
parser.add_option("-l", "--ligand", help="ligand file name")
parser.add_option("-g", "--group",  help="a list of receptor-ligand pairs  ")
(options, args) = parser.parse_args()
options =  options.__dict__

# assign args to pdblist
pdblist = []
if options['group'] is None:
    # one protein-ligand complex
    if options['receptor'] is None or options['ligand'] is None:
        print "Please specify protein and ligand"
    else:
        pdblist.append([options['receptor'], options['ligand']])
else:
    # one or more protein-ligand complex
    with open(options['group'], 'r') as f:
        for lines in f:
            pdblist.append(lines.split()[0:2])
print pdblist

modelDV.runDV(pdblist)


