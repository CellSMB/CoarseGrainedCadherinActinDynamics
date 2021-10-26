# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 16:49:38 2020

@author: qiliny
"""
import time
import numpy as np
import modules as mo
import newmodules as nw
from inputset import *
import ActinCadInter as ACI

# Initialize Cadherins
###############################################################################
# SYSTEM INITIALIZATION  
###############################################################################
# Col 0: Index,         Col 7: Cis donor        # -1 if off, index of its receptor if in cis           
# Col 1: x_pos,         Col 8: Cis receptor     # -1 if off, index of its donor if in cis               
# Col 2: y_pos,         Col 9: Trans Partner    # -1 if off, index of its partner if in trans  
# Col 3: Orientation,   Col 10: Stall State     # Stalled state due to cis. It can also be stalled when its trans partner is in cis
# Col 4: Mono,                                  # 0 if free, 1 if stalled 
# Col 5: X,             Col 11: Actin binding site  # -1 if off, index of the actin segment if connected         
# Col 6: S,             Col 12: Stall State     # 0 if free, if stalled or its trans or cis partner is in actin binding
#                       Col 13: ABD             # 1 if yes, if this E-cad is available for binding
# Actinseg = Actinseg_all_1
def InitCad():
    
    array = np.zeros((num_par,13)) 
    array[:,0] = np.arange(num_par)
    array[:,3] = np.random.rand(num_par)*2*np.pi 
    array[:,4] = 1 
    array[:,[7,8,9,11]] = -1

    for i in range(num_par):
        collision_1 = True
        
        while collision_1 == True:
            array[i,1:3] = np.random.rand(2)*side_len_domain
            
            array_all_hash = nw.CreateHash(mo.CreateChainingHashArray(array[:i+1,:],side_len_par_collision,num_part_side_collision))
            partition = (array[i,1]//side_len_par_collision+
                        (array[i,2]//side_len_par_collision)*num_part_side_collision)
            nearby_part_desti = nw.collision_hash[partition]
            
            nearby_Ecad_desti = ACI.CheckNearbyPartition(nearby_part_desti, array_all_hash) 
            collision_1 = ACI.CheckCollision(i, i, nearby_Ecad_desti, array, array)
    return array

def TransInteraction(M1, M2):
    ###########################################################################
    # Find Particles that meet conditions for trans and cis interaction
    ###########################################################################
    # Flow of the algorithm:
    # 1st part, partition the filed without considering periodic boundary
    # 2nd part, partition the filed with considering periodic boundary
    # 3rd part, combine results from both parts, and only keep 1 pairs. 
    trans_pairs = np.zeros((0,2))
    #--------------------------------------------------------------------------
    # Trans 
    #--------------------------------------------------------------------------
    
    M1_mono = M1[M1[:,4]==1,:]
    M2_mono = M2[M2[:,4]==1,:]
    
    # --------------1st part--------------
    # Without count Periodic Boundary Condition and inter partitions

    start_3 = time.time()
    M1_mono_hash = mo.CreateChainingHashArray(M1_mono,side_len_par_trans,num_part_side_trans)
    M2_mono_hash = mo.CreateChainingHashArray(M2_mono,side_len_par_trans,num_part_side_trans)
    time3 = time.time()-start_3
    
    start_4 = time.time()
    new_pairs = nw.FindBindingPairs(M1_mono_hash, M2_mono_hash, M1_mono, M2_mono, nw.trans_hash, trans_thre_sq, True)
    time4 = time.time()-start_4
#    print("new_pairs: {}".format(new_pairs))
    new_pairs = np.reshape(new_pairs,(-1,2))
    
    start_5 = time.time()
    if len(new_pairs!=0):
        trans_pairs = mo.AvoidDoubleInteraction(new_pairs)
    time5 = time.time()-start_5
    
    return trans_pairs

def CisInteraction(M1, M2):
    #--------------------------------------------------------------------------
    # Cis M1  
    #--------------------------------------------------------------------------
    # Flow of the Algorithm
    # 0th: Calculate position of all donor and receptor sites
    # 1st: find all cis pairs without counting PBS and inter partitions
    # 2nd: find all cis pairs counting PBS and inter partitions
    # 3rd: combine results from both parts, and only keep 1 pairs if see repeats. 

    M1_donor_bool = (M1[:,7]==-1)
    M1_recep_bool = (M1[:,8]==-1)
    M1_cis_pairs = np.zeros((0,2))
    
    if (M1_donor_bool.shape[0] != 0) & (M1_recep_bool.shape[0] != 0):

        M1_donor = M1[M1_donor_bool,:]
        M1_recep = M1[M1_recep_bool,:]
        
        M1_donor[:,1] += rad_ecad*np.cos(M1_donor[:,3])
        M1_donor[:,2] += rad_ecad*np.sin(M1_donor[:,3])
        M1_donor = mo.adjustPBC(M1_donor)
        
        M1_recep[:,1] -= rad_ecad*np.cos(M1_recep[:,3])
        M1_recep[:,2] -= rad_ecad*np.sin(M1_recep[:,3])
        M1_recep = mo.adjustPBC(M1_recep)
        # 1st part
        # Without counting Periodic Boundary Condition and inter partitions
        
        M1_donor_hash = mo.CreateChainingHashArray(M1_donor,side_len_par_cis,num_part_side_cis)
        
        M1_recep_hash = mo.CreateChainingHashArray(M1_recep,side_len_par_cis,num_part_side_cis)
        
        new_pairs = nw.FindBindingPairs(M1_donor_hash, M1_recep_hash, M1_donor, M1_recep, nw.cis_hash, cis_thre_sq)
        new_pairs = np.reshape(new_pairs,(-1,2))
        
        if len(new_pairs!=0):
            M1_cis_pairs = mo.AvoidDoubleInteraction(new_pairs)
 
    #--------------------------------------------------------------------------
    # Cis M2  
    #--------------------------------------------------------------------------
    # Flow of the Algorithm
    # 0th: Calculate position of all donor and receptor sites
    # 1st: find all cis pairs without counting PBS and inter partitions
    # 2nd: find all cis pairs counting PBS and inter partitions
    # 3rd: combine results from both parts, and only keep 1 pairs if see repeats. 
    M2_donor_bool = (M2[:,7]==-1)
    M2_recep_bool = (M2[:,8]==-1)
    M2_cis_pairs = np.zeros((0,2))

    if (len(M2_donor_bool) != 0) & (len(M2_recep_bool) != 0):
        
        M2_donor = M2[M2_donor_bool,:]
        M2_recep = M2[M2_recep_bool,:]
        
        M2_donor[:,1] += rad_ecad*np.cos(M2_donor[:,3])
        M2_donor[:,2] += rad_ecad*np.sin(M2_donor[:,3])
        M2_donor = mo.adjustPBC(M2_donor)
        
        M2_recep[:,1] -= rad_ecad*np.cos(M2_recep[:,3])
        M2_recep[:,2] -= rad_ecad*np.sin(M2_recep[:,3])
        M2_recep = mo.adjustPBC(M2_recep)

        # 1st part
        # Without counting Periodic Boundary Condition and inter partitions

        M2_donor_hash = mo.CreateChainingHashArray(M2_donor,side_len_par_cis,num_part_side_cis)
        
        M2_recep_hash = mo.CreateChainingHashArray(M2_recep,side_len_par_cis,num_part_side_cis)
        
        new_pairs = nw.FindBindingPairs(M2_donor_hash, M2_recep_hash, M2_donor, M2_recep, nw.cis_hash, cis_thre_sq)
        new_pairs = np.reshape(new_pairs,(-1,2))
        
        if len(new_pairs!=0):
            M2_cis_pairs = mo.AvoidDoubleInteraction(new_pairs)
    
    return M1_cis_pairs, M2_cis_pairs

    ###########################################################################
    # Update Association and Dissociation
    ###########################################################################

def UpdateStates(trans_pairs, M1_cis_pairs, M2_cis_pairs, M1, M2):
    #--------------------------------------------------------------------------
    # trans formation and dissociation
    #--------------------------------------------------------------------------
      
    # --------------MM to X or S--------------

    MM_prob = np.random.random(len(trans_pairs))
    new_MM_X = trans_pairs[MM_prob<p_k_mx,:]
    new_MM_S = trans_pairs[(MM_prob>p_k_mx)&(MM_prob<(p_k_mx+p_k_ms)),:]
    
    # --------------X to S or MM--------------
    All_X = M1[M1[:,5]==1,0].astype(int)
    X_prob = np.random.random(len(All_X))
    new_X_S = All_X[X_prob<p_k_xs]
    new_X_MM = All_X[(X_prob>p_k_xs)&(X_prob<(p_k_xs+p_k_xm))]

    # -------------- 1. S to X or MM --------------
    
    All_S = M1[M1[:,6]==1,0].astype(int)
    S_prob = np.random.random(len(All_S))
    new_S_X = All_S[(S_prob<p_k_sx)]
    new_S_M = All_S[(S_prob<(p_k_sx+p_k_sm))&(S_prob>p_k_sx)]
    
    # --------------Perform Changes--------------
    # Perform only dissociation
    
    if len(new_MM_X)!=0:
        M1,M2 = mo.StateChangeMM(new_MM_X,4,5,M1,M2)
    
    if len(new_MM_S)!=0:
        M1,M2 = mo.StateChangeMM(new_MM_S,4,6,M1,M2)

    if len(new_X_S)!=0:
        M1,M2 = mo.StateChangeTrans(new_X_S,5,6,M1,M2)

    if len(new_X_MM)!=0:
        M1,M2 = mo.StateChangeTrans(new_X_MM,5,4,M1,M2)

    if len(new_S_X)>0:
        M1,M2 = mo.StateChangeTrans(new_S_X, 6, 5, M1, M2)
    
    if len(new_S_M)>0:
        M1,M2 = mo.StateChangeTrans(new_S_M, 6, 4, M1, M2)

    #--------------------------------------------------------------------------
    # cis formation and dissociation
    #--------------------------------------------------------------------------
    # -------------- 1. M+M or S+S or S+M to cis --------------
    M1_new_cis = []
    
    if len(M1_cis_pairs) != 0:
        Cis_prob_M1  = np.random.random(len(M1_cis_pairs))
        M1_new_cis   = mo.NewCis(M1,Cis_prob_M1,M1_cis_pairs)
    
    M2_new_cis = []
    if len(M2_cis_pairs) != 0:
        Cis_prob_M2  = np.random.random(len(M2_cis_pairs))
        M2_new_cis   = mo.NewCis(M2,Cis_prob_M2,M2_cis_pairs)

    # -------------- 2. cis to S --------------
    All_Cis_M1 = M1[M1[:,7]!=-1,0]          # cis donor site of the Ecad
    All_Cis_M2 = M2[M2[:,7]!=-1,0]          # cis donor site of the Ecad
    
    Cis_prob_M1 = np.random.random(len(All_Cis_M1))
    Cis_prob_M2 = np.random.random(len(All_Cis_M2))
    
    new_dis_Cis_M1 = All_Cis_M1[Cis_prob_M1<p_k_ciss].astype(int)
    new_dis_Cis_M2 = All_Cis_M2[Cis_prob_M2<p_k_ciss].astype(int)
    
    # -------------- 3. Perform only dissociation --------------
    # M and S to cis
    if len(M1_new_cis)>0:
        M1 = mo.StateChangeCis(M1, M1_new_cis[:,0], M1_new_cis[:,1])
    
    if len(M2_new_cis)>0:
        M2 = mo.StateChangeCis(M2, M2_new_cis[:,0], M2_new_cis[:,1])
    
    # Cis to S
    if len(new_dis_Cis_M1) > 0:
        M1[M1[new_dis_Cis_M1,7].astype(int),8] = -1       # set recep site of its cis partner to be available
        M1[new_dis_Cis_M1,7] = -1                         # set its donor site to be available
    
    if len(new_dis_Cis_M2) > 0:    
        M2[M2[new_dis_Cis_M2,7].astype(int),8] = -1
        M2[new_dis_Cis_M2,7] = -1
        
    return M1, M2, new_MM_X, M1_new_cis, M2_new_cis

#--------------------------------------------------------------------------
# Orientation Update
#--------------------------------------------------------------------------
def OriUpdate(M1,M2):
    
    if Rotation_withActinBinding:
        # All free monomers
        M1_mono_moveable_index = M1[(M1[:,4]==1)&(M1[:,10]==0),0].astype(int)
        M2_mono_moveable_index = M2[(M2[:,4]==1)&(M2[:,10]==0),0].astype(int)
        
        # All free trans dimers
        M1_trans_moveable_index = M1[((M1[:,5]==1)|(M1[:,6]==1))&(M1[:,10]==0),0].astype(int)
        M2_trans_moveable_index = M1[M1_trans_moveable_index,9].astype(int)
    else:
        # All free monomers
        M1_mono_moveable_index = M1[(M1[:,4]==1)&(M1[:,10]==0)&(M1[:,11]==-1),0].astype(int)
        M2_mono_moveable_index = M2[(M2[:,4]==1)&(M2[:,10]==0)&(M2[:,11]==-1),0].astype(int)
        
        # All free trans dimers
        M1_trans_moveable_index = M1[((M1[:,5]==1)|(M1[:,6]==1))&(M1[:,10]==0)&(M1[:,11]==-1),0].astype(int)
        M2_trans_moveable_index = M1[M1_trans_moveable_index,9].astype(int)
        
    M1[M1_mono_moveable_index,3] = np.random.rand(len(M1_mono_moveable_index))*2*np.pi 
    M2[M2_mono_moveable_index,3] = np.random.rand(len(M2_mono_moveable_index))*2*np.pi 

    M1[M1_trans_moveable_index,3] = np.random.rand(len(M1_trans_moveable_index))*2*np.pi 
    M2[M2_trans_moveable_index,3] = M1[M1_trans_moveable_index,3]-np.pi/2
    M2[M2_trans_moveable_index[M2[M2_trans_moveable_index,3] < 0],3] += np.pi*2

    return M1,M2

#--------------------------------------------------------------------------
# Position Update
#--------------------------------------------------------------------------
def PosUpdate(M1,M2):
    # All free monomers
    M1_mono_moveable_index = M1[(M1[:,4]==1)&(M1[:,10]==0)&(M1[:,11]==-1),0].astype(int)
    M2_mono_moveable_index = M2[(M2[:,4]==1)&(M2[:,10]==0)&(M2[:,11]==-1),0].astype(int)
    
    # All free trans dimers
    M1_trans_moveable_index = M1[((M1[:,5]==1)|(M1[:,6]==1))&(M1[:,10]==0)&(M1[:,11]==-1),0].astype(int)
    M2_trans_moveable_index = M1[M1_trans_moveable_index,9].astype(int)

    mu = 0; sigma = 1

    M1[M1_mono_moveable_index,1] += steplength*np.random.normal(mu, sigma, len(M1_mono_moveable_index))
    M1[M1_mono_moveable_index,2] += steplength*np.random.normal(mu, sigma, len(M1_mono_moveable_index))
    
    M2[M2_mono_moveable_index,1] += steplength*np.random.normal(mu, sigma, len(M2_mono_moveable_index))
    M2[M2_mono_moveable_index,2] += steplength*np.random.normal(mu, sigma, len(M2_mono_moveable_index))
    
    # trans dimers' moving directions
    trans_delta_x = steplength*np.random.normal(mu, sigma, len(M1_trans_moveable_index))
    trans_delta_y = steplength*np.random.normal(mu, sigma, len(M1_trans_moveable_index))
    
    M1[M1_trans_moveable_index,1] += trans_delta_x
    M1[M1_trans_moveable_index,2] += trans_delta_y
    
    M2[M2_trans_moveable_index,1] += trans_delta_x
    M2[M2_trans_moveable_index,2] += trans_delta_y
    
    M1 = mo.adjustPBC(M1)
    M2 = mo.adjustPBC(M2)
    
    return M1, M2
    