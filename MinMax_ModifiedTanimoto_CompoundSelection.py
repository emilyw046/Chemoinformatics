# -*- coding: utf-8 -*-
"""
Created on Sun Mar 12 10:54:55 2023

Project: Compare using the Modified Tanimoto Coefficient (M.A. Fligner et al, 2002) 
to the Tanimoto Coefficient for selecting diverse compounds via the MinMax Selection 
Algorithm (Snarey et al, 1997). Analysis of diversity of subset chosen is based on 
M.A. Fligner's methods of examining feature number and average pairwise dissimilarity
in the chosen subset. The use of the Modified Tanimoto Coefficient can be used to 
select an intial compound, select other compounds for the subset based on the intial 
compound, and evaluate the compound selection in the final subset. Comparison of Modified Tanimoto
vs Tanimoto performance, however, should focus on the differences in the selection process.


@author: emily
"""
from rdkit import Chem
from rdkit.Chem import MACCSkeys
from rdkit import DataStructs
from rdkit.Chem import AllChem
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from functools import partial




def find_MT(molA, molB):
    '''
    
       Find Modified Tanimoto Coefficient based on input of two molecule fingerprints
    
       Parameters
       ----------
       molA_fp : fingerprint (explicit bit vector)
           fingerprint of 1 molecule.
       molB_fp : fingerprint (explicit bit vector)
           fingerprint of another molecule.
    
       Returns
       -------
       MT : float
           Modified Tanimoto Coefficient (Similarity measure).
    
    '''

    
    a = sum(molA)
    b = sum(molB)
    c=d=0
    for m in range(len(molA)): 
     
        if molA[m] == 1 and molB[m] ==1:
            c +=1 
            # add 1 to c if 1 bit is found in molA and molB
        if molA[m] ==0 and molB[m] ==0:
            d+=1 
    
    p = (a+b)/(2*(a+b-c+d))
     #calculate proportion of "1" bits in both vectors combined
    T1 = c/(a+b-c) 
    #conventional Tanimoto coefficient
    T0 = d/(a+b-2*c+d) 
    #Similar to Tanimoto, but based on 0 bits rather than 1 bits
    MT = ((2-p)/3)*T1 + ((1+p)/3)*T0 #modified Tanimoto
    
    return MT

def fp_find(mset,  fp ='maccs', radius = 2):
    '''
    

    Parameters
    ----------
    mset : mol set
         set of RDKit molecules.

    fp : str, optional
        Desired fingerprint method. The default is 'maccs'.
    radius : int, optional
        only for morgan fingerprints. The default is 2.

    Returns
    -------
    fp_object: function
        Which fingerprinting method to use for compound selection


    '''
    
    if fp == 'maccs':
        fp_object = MACCSkeys.GenMACCSKeys
     
            
    elif fp == 'topol':
        fp_object = Chem.RDKFingerprint
        
        
    elif fp == 'morgan':
        fp_object = partial(AllChem.GetMorganFingerprintAsBitVect, radius = radius, nbits = 1024)
    
        
        
    return fp_object

def ChooseSubset(mset, fp_object, sim = 'T', method = 'MostDissim'):
    '''
    

    Parameters
    ----------
    mset : Mol set
        Set of molecules
        .
    fp_object : function
        Which fingerprinting method to use.
        Options are 'maccs', 'topol', or 'morgan'
    sim : str, optional
        Which method to use to determine similarity. The default is 'T'.
        Only options are 'T' or 'MT'
    method : str, optional
        Method of finding first molecule for subset. The default is 'MostDissim'.

    Returns
    -------
    subset_fp : BitVector list
        list of molecules' fingerprints in initial subset.
    subset : Mol Object
        Molecules in subset
    fps : BitVector List
        Fingerprints of all molecules in whole set.

    '''
    
    
    fps = [fp_object(m) for m in mset]
    dissim_set = []
    sum_diss = []
    for j in range(len(fps)):  #compare every molecule in subset to every other molecule in mset
        for l in range(len(fps)):
            if j < l: #don't compare the same molecule to itself and don't repeat
                if sim == 'MT' and method == 'MostDissim':
                    dissim_set.append(1 - (find_MT(fps[j],fps[l])) )#modified tanimoto dissimilarity
                elif sim == 'T' and method == 'MostDissim':
                    dissim_set.append( 1- (DataStructs.FingerprintSimilarity(fps[j],fps[l])))
                elif sim == 'MT' and method == 'MostSim':
                    dissim_set.append(find_MT(fps[j],fps[l]))
                elif sim == 'T' and method == 'MostDissim':
                    dissim_set.append(DataStructs.FingerprintSimilarity(fps[j],fps[l]))
        sum_diss.append(np.sum(dissim_set))
    index = sum_diss.index(max(sum_diss))
    subset_fp = [fps[index]]
    subset = [mset[index]]
    
    return subset_fp, subset, fps
    
        
        
    

   
def MinMax(mset, subset, fps, subset_fp, n=20, sim = 'T' ):
    '''
    

    Parameters
    ----------
    mset : mol set
        Set of molecules
    subset : mol set
        Initial subset of molcules
        .
    fps : BitVector
        fingerprints of whole set.
    subset_fp : BitVector
        fingerprints of initial subset.
    n : int, optional
        Desired number of molecules in subset. The default is 20.
    sim : str, optional
        Type of dissimilarity metric to calculate. The default is 'T'.
        Options are 'T' or 'MT'

    Returns
    -------
    names_subset : list
        list of all names in new subset.
    subset_fp : BitVector
        all fingerprints in new subset.

    '''
    #remove initial seed molecule fingerprints from set of all fingerprints
    fps = [e for e in fps if e not in subset_fp]
    
    #find names of molecules in mset
    names_mset_og = [x.GetProp('_Name') for x in mset]
    names_subset =[x.GetProp('_Name') for x in subset] #names in subset

    #remove initial seed molecule names from set of all names
    names_mset = [m for m in names_mset_og if m not in names_subset]
   

   
    while len(subset_fp) < n : #while subset_fp is less than desired length
        min_sim = [] #establish minimum dissimilarity list
     
        for m in range(len(fps)):
            sub_lst = [] #sublist of dissimilarities between each subset molecule and the mth molecule in mset
            for i in range(len(subset_fp)): 
                if sim == 'MT':
                    y = 1 - (find_MT(fps[m],subset_fp[i])) #modified tanimoto dissimilarity
                elif sim == 'T':
                    y =  1- (DataStructs.FingerprintSimilarity(fps[m],subset_fp[i])) #reg tanimoto dissimilarity
                
                sub_lst.append(y) #append dissimilarity metric to sub list
                
            minn = min(sub_lst) #find min of two dissimilarity in sublist
                
    
            min_sim.append(minn) #add minimum to minimum similarity list
           
        
        i_nm = min_sim.index(max(min_sim)) #find index of the greatest minimum similarity 
       
        subset_fp.append(fps[i_nm]) #add fp with index of greatest minimum similarity to subset fp
        
        #subset.append(mset[i_nm]) #add molecule with index of greatest min sim to subset
    
        fps.remove(fps[i_nm]) #remove fp with index of greatest minimum similarity from fps (whole mset)
        
        names_subset.append(names_mset[i_nm]) #add name of molecule to subset of names
        
        names_mset.remove(names_mset[i_nm]) #remove name of molecule from mset of names
        
    
    return names_subset,  subset_fp
    
 
def metrics(subset_fp, names_subset,sim = 'T'):
    '''
    

    Parameters
    ----------
    
    subset_fp : BitVector
        all fingerprints in new subset.
    names_subset : list
            list of all names in new subset.
    sim : str, optional
        Which dissimilarity metric to use. The default is 'T'.
        Other option is 'MT'

    Returns
    -------
    features : list
        list of number of features of all molecules in subset.
    dissim_subset : TYPE
        Dissimilarity measures of all molecules in subset.

    list
        List of various metrics calculated from the dissim_subset and features lists.

    '''
    dissim_subset = []   
    features = [] #also want to find proportion of features predicted by MT vs T
    sub = []
    for j in range(len(subset_fp)):  #compare every molecule in subset to every other molecule in subset
        for l in range(len(subset_fp)):
            if j < l: #don't compare the same molecule to itself and don't repeat
                if sim == 'MT':
                    dissim_subset.append(1 - (find_MT(subset_fp[j],subset_fp[l])) )#modified tanimoto dissimilarity
                elif sim == 'T':
                    dissim_subset.append( 1- (DataStructs.FingerprintSimilarity(subset_fp[j],subset_fp[l])))
                
                           
                features.append(sum(subset_fp[j])/(len(subset_fp[j]))) #num 1 bits/all possible bits
                
                sub.append([names_subset[j], names_subset[l]])
                    
                
    avg_dissim = np.mean(dissim_subset) #average dissimilarity between all molecules in subset
    
    std_dissim =np.std(dissim_subset)
    
    min_dissim = np.min(dissim_subset) #minimum dissimilarity between all molecules in subset
        
    propfeatures =np.mean(features) #average of all features present in subset (average proportion of 1s in bit vector)
    
    std_propfeatures =np.std(features)
    
    
    return features, dissim_subset, [avg_dissim, std_dissim, min_dissim, propfeatures, std_propfeatures]

if __name__ == '__main__': 
     
    
     
    input_file = 'micropollutants_193.sdf' 
    mset_M = Chem.SDMolSupplier(input_file)
    
    fp_object = fp_find(mset_M, fp = 'maccs')
    
    MT_subset_fp, subset_M, MT_fps = ChooseSubset(mset_M, fp_object, sim = 'T', method = 'MostDissim')
    #MT and T select the same initial molecule! for micropollutants, molecule was VERAPAMIL
    #final results (metrics) pretty close to random selection as subset
    
    MTnames, subset_fp_MT = MinMax(mset_M, subset_M, MT_fps, MT_subset_fp, sim = 'MT')
    MT_feat, MT_dissub, MT_metrics = metrics(MT_subset_fp, MTnames, sim = 'MT')
    TMT_feat, TMT_dissub,  TMT_metrics = metrics(MT_subset_fp, MTnames, sim = 'T')
    
    
    mset_T = Chem.SDMolSupplier(input_file)
    T_subset_fp, subset_T, T_fps = ChooseSubset(mset_T, fp_object, sim = 'T', method = 'MostDissim')
    
    Tnames,  subset_fp_T = MinMax(mset_T,subset_T, T_fps, T_subset_fp, sim ='T')
    T_feat, T_dissub,  T_metrics = metrics(T_subset_fp, Tnames, sim = 'T')
    MTT_feat, MTT_dissub,  MTT_metrics = metrics(T_subset_fp, Tnames, sim = 'MT')
    

   
    
    fig1, axes1 = plt.subplots(2,2)
    fig1.suptitle('Dissimilarity Distributions (maccs): MT vs T')
    axes1[0,0].hist(MT_dissub, bins = 10)
    axes1[0,0].set_xlabel('MT Selection, MT Dissimilarity')
    axes1[0,0].set_ylabel('Count')
    axes1[1,0].hist(TMT_dissub,bins=10)
    axes1[1,0].set_xlabel('MT Selection,T Dissimilarity')
    axes1[1,0].set_ylabel('Count')
    axes1[0,1].hist(MTT_dissub,bins=10)
    axes1[0,1].set_xlabel('T Selection, MT Dissimilarity')
    axes1[0,1].set_ylabel('Count')
    axes1[1,1].hist(T_dissub,bins=10)
    axes1[1,1].set_xlabel('T Selection, T Dissimilarity')
    axes1[1,1].set_ylabel('Count')
    plt.tight_layout()
    #plt.savefig('Morgan_distribution.png')
   
   
      
    table = {'Average Pairwise Dissimilarity': [MT_metrics[0], MTT_metrics[0], TMT_metrics[0],T_metrics[0]],
             'StdDev Pairwise Dissimilarity': [MT_metrics[1], MTT_metrics[1],TMT_metrics[1],T_metrics[1]],
             'Minimum Pairwise Dissimilarity': [MT_metrics[2],  MTT_metrics[2],TMT_metrics[2],T_metrics[2]], 
             'Proportion of Features Present': [MT_metrics[3], MTT_metrics[3],TMT_metrics[3],T_metrics[3]],
             'StdDev Features': [MT_metrics[4], MTT_metrics[4],TMT_metrics[4],T_metrics[4]]
             }
    table_df = pd.DataFrame(table, index = ["MT Selection, MT Eval", "T Selection, MT Eval", 
                                            "MT Selection, T Eval", "T Selction, T Eval"])
    print(table_df)
    
    
    #table_df.to_csv('Metrics_Table_Morgan.csv')
