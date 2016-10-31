DeltaVina
======

A scoring function for protein-ligand rescoring.

Project Organization
------------
    ├── LICENSE
    ├── README.md          <- README file.
    │
    ├── docs               <- A default Sphinx project; see sphinx-doc.org for details
    │
    ├── bin                <- Master script to run deltavina
    │   └── dvrf20.py
    │
    ├── src                <- Source code for use in this project.
    │   ├── __init__.py    <- Makes src a Python module
    │   │
    │   ├── features       <- Scripts to turn raw data into features for modeling
    │   │   ├── pharma.py  
    │   │   ├── featureSASA.py 
    │   │   ├── featureVina.py 
    │   │   └── tests
    │   │
    │   └── models         <- Trained model rffit.rda and scripts to run model
    │       ├── rffit.rda
    │       ├── modelDV.py
    │       └── tests
    │
    └── examples           <- Examples
    
