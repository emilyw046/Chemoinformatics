# -*- coding: utf-8 -*-
"""
Created on Mon Feb 20 12:34:01 2023

Project: Program that allows a user to input a query (target) molecule of interest
and search for similar structures in a specified set of molecules



@author: emily
"""

from rdkit import Chem

from rdkit.Chem import MACCSkeys
from rdkit import DataStructs
from rdkit.Chem import AllChem

import pandas as pd
def findNN(mol, mset, threshold = 0.6, nn = None, fp = 'maccs', radius = 2):
    '''
    Given the following parameters, find nearest neighbors of target molecule in query dataset


    Parameters
    ----------
    mol : mol
        1 target molecule.
    mset : mol
        set of molecules.
    threshold : float, optional
        minimum similarity threshold. The default is 0.6.
    nn : int, optional
        # of nearest neighbors to calculate. The default is None.
    fp : str, optional
        Fingerprint type. The default is 'maccs'.
    radius : int, optional
         radius for use in Morgan footprints . The default is 2.

    Returns
    -------
    df : dataframe
        Dataframe of name, SMILES, and similarity metrics of nn nearest neighbors to target molecule.
'''
    
    
    #first, determine fingerprints for mset and mol based on user input of fingerprint type
    if fp == 'maccs':
        fps = [MACCSkeys.GenMACCSKeys(x) for x in mset]
        target_fp = MACCSkeys.GenMACCSKeys(mol)
        
        
    elif fp == 'topol':
        fps = [Chem.RDKFingerprint(x) for x in mset]
        target_fp = Chem.RDKFingerprint(mol)
        
        
    elif fp == 'morgan':
        fps = [AllChem.GetMorganFingerprintAsBitVect(x,radius,nBits=1024) for x in mset]
        target_fp = AllChem.GetMorganFingerprintAsBitVect(mol,radius,nBits=1024)
    
    
    s = []  #similarity empty list
    names = [] #names empty list
    SMILES = [] #SMILES empty list
    
    for x in range(len(mset)): #loop through length of mset
        y = DataStructs.FingerprintSimilarity(fps[x],target_fp)
        print(y) #determine Tanimoto similarity bt target and mol in mset
        if y > threshold: #if similarity greater than user inputted threshold
            names.append(mset[x].GetProp('_Name')) #append name of molecule
            s.append(y) # append similarity of moelcule
            SMILES.append(Chem.MolToSmiles(mset[x]))  #append SMILES of molecule
   
    
    if nn is not None:#specify value of nearest neighbors
        neighb = nn
    else:
        neighb = len(s) #if nn is None, neighb is equal to length of s (all similarities above threshold)
        
    d = {'Name': names, 'SMILES': SMILES, 'Similarity':s} #create dictionary
    df = pd.DataFrame(d) #create dataframe from dictionary
    
    #sort values in dictionary from highest to lowest similarity, and only include neighb number of rows
    df = df.sort_values('Similarity', ascending=False)[:neighb]
    
    
    return df #tada!
    
if __name__ == '__main__': 
     
    target_smiles = 'CN1C=NC2=C1C(=O)N(C)C(=O)N2C' #caffeine 
    target = Chem.MolFromSmiles(target_smiles) 
     
    input_file = 'micropollutants_193.sdf' 
    mset = Chem.SDMolSupplier(input_file) 
     
    df = findNN(target, mset, threshold=0.5, fp='morgan') 
    print(df)