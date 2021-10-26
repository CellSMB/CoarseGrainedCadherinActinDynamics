# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 14:38:57 2020

@author: qiliny
"""
import numpy as np
import itertools
import ctypes
import sys
from inputset import *
import copy
# array = M2_mono; side_len_par = side_len_par_trans;
# num_partition_sqrt = num_partition_sqrt_trans; len_hash = len_hash_trans_1

def CreateChainingHashArray(array,side_len_par,num_partition_sqrt):   
    partition = (array[:,1]//side_len_par+
               (array[:,2]//side_len_par)*num_partition_sqrt).astype(int)
    return partition

def InterPar(array, side_len_par):
    array_copy1 = np.copy(array)
    array_copy1[:,1] = array_copy1[:,1] + side_len_par/2
    array_copy1[:,2] = array_copy1[:,2] + side_len_par/2
    
    array_copy2 = np.copy(array_copy1)
    array_copy3 = np.copy(array_copy1)
    array_copy4 = np.copy(array_copy1)
    
    array_copy2[:,1] = array_copy2[:,1] + side_len_par/2
    
    array_copy3[:,2] = array_copy3[:,2] + side_len_par/2
    
    array_copy4[:,1] = array_copy4[:,1] + side_len_par/2
    array_copy4[:,2] = array_copy4[:,2] + side_len_par/2

    return array_copy1, array_copy2, array_copy3, array_copy4

# array = Actin_copy; side_len_par = side_len_par_EcadActin
# Count the Periodic Boundary Condition
def PBC_DoubleEdge(array, side_len_par):
    # Take periodic boundary into account

    most_left = array[array[:,1]<(side_len_par),:]
    most_left[:,1] += side_len_domain

    most_right = array[array[:,1]>(side_len_domain-side_len_par),:]
    most_right[:,1] -= side_len_domain

    most_bot = array[array[:,2]<(side_len_par),:]
    most_bot[:,2] += side_len_domain
    
    most_top = array[array[:,2]>(side_len_domain-side_len_par),:]
    most_top[:,2] -= side_len_domain
    
    # Shift particle to make sure corrdinates are positive
    # Array1 shift all particle to top right by a/2 based on original corrdinate
    # This makes all the corrdinate to be possitive for partition to work.

    array_copy = np.concatenate((array,most_left,most_right,most_bot,most_top), axis = 0)
    array_original_index = array_copy[:,0].astype(int)
    array_copy[:,0] = np.arange(len(array_copy))
    array_copy[:,1] = array_copy[:,1] + side_len_par
    array_copy[:,2] = array_copy[:,2] + side_len_par

    return array_copy,array_original_index

def PBC_Shift(array, side_len_par):
    # Shift particle to make sure corrdinates are positive
    # Array1 shift all particle to top right by a/2 based on original corrdinate

    array[:,1:3] += side_len_par

    return array

# np.where((array_copy1[:,1]>side_len_domain)|(array_copy1[:,2]>side_len_domain))

# Find binding pairs for either trans, cis or cadherin/actin
# array1,array2 are M1 and M2 for trans binding
# array1,array2 are M1/M1 or M2/M2 for cis binding
# array1,array2 are M1/actin or M2/actin for cad/actin binding
def FindBindingPairs(tmp_pairs,array0,array1, dist_thre, trans = False):
        true_pairs = []
        
        for i in range(len(tmp_pairs)):
            i0 = tmp_pairs[i,0]; i1 = tmp_pairs[i,1];
            
            x0 = array0[i0,1]; x1 = array1[i1,1]; 
            
            y0 = array0[i0,2]; y1 = array1[i1,2]; 
            
            angle0 = array0[i0,3]; angle1 = array1[i1,3];
            
            distance = np.sqrt((x0-x1)**2+(y0-y1)**2)
            
            judge_angle = False
            if trans:
                a0 = abs(angle0-(angle1+np.pi/2))
                a360 = abs(abs(angle0-(angle1+np.pi/2))-2*np.pi)
                if (a0 < angle_thre_s) or (a360 < angle_thre_s):
                    judge_angle = True
            else:
                a0 = abs(angle0-angle1)
                if (a0 < angle_thre_s) or (a0 > angle_thre_l):
                    judge_angle = True
                
            if (distance < dist_thre) & judge_angle:
                true_pairs.append([i0,i1])
        return true_pairs

def adjustPBC(array_ori):
    array = np.copy(array_ori)
    array[array[:,1]>side_len_domain,1] -= side_len_domain
    array[array[:,2]>side_len_domain,2] -= side_len_domain
    array[array[:,1]<0,1] += side_len_domain
    array[array[:,2]<0,2] += side_len_domain
    return array

# Find newly associated cis
# Find all stalled Ecad due to cis
    
# array,Cis_prob,Cis_pairs,side1_new_stall,side2_new_stall = M1,Cis_prob_M1,M1_cis_pairs,M1_new_stall,M2_new_stall
def NewCis(array,Cis_prob,Cis_pairs):
    
    successful_cis = Cis_prob<p_k_scis
    new_cis = Cis_pairs[successful_cis,:]

    # if len(new_S_cis_donor) != 0:
    #     side1_new_stall+[int(i) for i in new_S_cis_donor]
    #     side1_new_stall+[int(i) for i in new_S_cis_recep]
    #     side2_new_stall+[int(i) for i in array[new_S_cis_donor,9]]
    #     side2_new_stall+[int(i) for i in array[new_S_cis_recep,9]]
        
        
#        side1_new_stall.append(new_S_cis_donor.astype(int))
#        side1_new_stall.append(new_S_cis_recep.astype(int))
#        side2_new_stall.append(array[new_S_cis_donor,9].astype(int))
#        side2_new_stall.append(array[new_S_cis_recep,9].astype(int))
    
    return np.array(new_cis)
# ,side1_new_stall,side2_new_stall
    
def StateChangeMM(mono_i, i_off, i_on, array1, array2, pair_c=9):
    array1[mono_i[:,0],i_off] = 0; array2[mono_i[:,1],i_off] = 0
    array1[mono_i[:,0],i_on] = 1;  array2[mono_i[:,1],i_on] = 1
    
    array1[mono_i[:,0],pair_c] = mono_i[:,1]; array2[mono_i[:,1],pair_c] = mono_i[:,0]

    return array1, array2

def StateChangeTrans(mono_i, i_off, i_on, array1, array2, pair_col=9):
    array2_pair = array1[mono_i,9].astype(int)
    array1[mono_i,i_off] = 0; array2[array2_pair,i_off] = 0
    array1[mono_i,i_on] = 1;  array2[array2_pair,i_on] = 1
    
    # Set trans pair to -1, if the trans dimer break
    if i_on == 4:
        array1[mono_i,pair_col] = -1; array2[array2_pair,pair_col] = -1

    return array1, array2

def StateChangeCis(array, new_S_cis_donor, new_S_cis_recep):
    array[new_S_cis_donor, 7] = new_S_cis_recep
    array[new_S_cis_recep, 8] = new_S_cis_donor
    
    return array

def AvoidDoubleInteraction(array):
    _,index_1 = np.unique(array[:,1],return_index = True)
    array = array[index_1,:]
    _,index_0 = np.unique(array[:,0],return_index = True)
    array = array[index_0,:]
    
    return array


# array1,array2 = M1_copy_hash_1,Actin_copy_hash_1
def Possible_Pairs(array1,array2):
    # c concatenate Ecad and actin segments
    # Ecad index are turn be negative numbers
    # col 0: partition index, col 1: Particle index
    c = np.c_[np.r_[array1,array2],
              np.r_[-np.arange(len(array1))-100,np.arange(len(array2))]]
    
    # sorted array based on the partition index
    # 'stable' makes sure that the index of Ecad appear first
    e = c[np.argsort(c[:,0],kind = 'stable')] 
    i = 0
    all_index = []       
    # record the pairs in the same partition
    while i < len(e)-1:
        # if the current index is for Ecad
        if e[i,1] < 0:
            Ecad = i
            # when next one is not Ecad and Ecad is in the same partition as Actin 
            while (i < len(e)-1) and (e[i+1,0] == e[Ecad,0]) and (e[i+1,1]>0):
                all_index.append([-e[Ecad,1]-100,e[i+1,1]])
                i+=1
        i+=1
    # The index returned here are index where the PBC is counted in
    return np.reshape(np.array(all_index),(-1,2)).astype(int)

# pairs,array1,array2,thre = potential_pairs,M1_copy_4,Actin_copy_4,EcadActin_thre
def FindEcadActinPairs(pairs,array1,array2,thre):
    
    x1 = array1[pairs[:,0],1]; y1 = array1[pairs[:,0],2];
    
    x2 = array2[pairs[:,1],1]; y2 = array2[pairs[:,1],2];
    
    dist = np.sqrt((x1-x2)**2+(y1-y2)**2)
    
    return pairs[dist>thre,:]

class Point:
    def __init__(self, x, y): 
        self.x = x 
        self.y = y 
        
def ccw(A,B,C):
    return (C.y-A.y)*(B.x-A.x) > (B.y-A.y)*(C.x-A.x)

def intersect(A,B,C,D):
    return ccw(A,C,D) != ccw(B,C,D) and ccw(A,B,C) != ccw(A,B,D)

# pairs, array1, array2, array3 = potential_pairs, Ecad_movable_4, Ecad_desti_4, Actin_move_4

def FindIntersections(pairs, array1, array2, array3):
    prev_actin = array3[pairs[:,1],4].astype(int)
    next_actin = array3[pairs[:,1],5].astype(int)
    all_pairs = []
    
    # All pairs between E-cad and actin segments
    for i in range(len(pairs)):
        new_list_1 = [pairs[i,0],prev_actin[i],pairs[i,1]]
        new_list_2 = [pairs[i,0],pairs[i,1],next_actin[i]]
        all_pairs.append(new_list_1)
        all_pairs.append(new_list_2)
    all_pairs = np.unique(np.array(all_pairs),axis = 0)
    all_pairs = np.array(all_pairs)
    
    # Remove pairs that the actin is on the pointed end
    all_pairs = all_pairs[(all_pairs[:,1]!=-1)&(all_pairs[:,2]!=-1),:]
    
    if max(all_pairs[:,0]) > len(array1):
        print("error")
        print("Ecad_movable_4: ",len(array1))
        print("all_pairs: ",max(all_pairs[:,0]))
        sys.exit()

    # Determine if two line segments intersect
    intersection_index = []
    for i in range(len(all_pairs)):
        if all_pairs[i,0] == 630 or i ==630:
            print(len(all_pairs))
            print(i)
            print("Ecad_movable_4: ",len(array1))
            print("all_pairs: ",max(all_pairs[:,0]))

        cad1 = Point(array1[all_pairs[i,0],1],array1[all_pairs[i,0],2])
        cad2 = Point(array2[all_pairs[i,0],1],array2[all_pairs[i,0],2])
        Actin1 = Point(array3[all_pairs[i,1],1],array3[all_pairs[i,1],2])
        Actin2 = Point(array3[all_pairs[i,2],1],array3[all_pairs[i,2],2])
        intersection = intersect(cad1,cad2,Actin1,Actin2)
        if intersection:
            intersection_index.append(i)
    all_intersections = all_pairs[intersection_index,:]
    
    return all_intersections

# Determine the final destination of each ecadherin 
# Mode 1: Reject the movement if hop failed, adjust the destination to be the original position
# Mode 2: If hop failed, adjust the destination to be the intersection.

#all_inter = intersections; Ecad_ori = M1_movable; Ecad_desti = M1_desti; Ecad_desti_4 = M1_desti_4;

def FinalDesti(all_inter, Ecad_ori, Ecad_desti, Ecad_desti_4):
    # shift back to original index
    failed_hop_ori = Ecad_desti_4[failed_hop,-1].astype(int)
    # adjust destinations
    Ecad_desti[failed_hop_ori,1:3] = Ecad_ori[failed_hop_ori,1:3]
    
    return Ecad_desti

#ori_index = ori_index_M1

def appendindexhash(ori_index):
    index_hash = copy.deepcopy(par_index_hash)
    for i in range(num_par,len(ori_index),1):
        index_hash[ori_index[i]].append(i)
    return index_hash

def MSD_data(M1, M2, M1_prev, M2_prev, M1_msd, M2_msd):
    
    M1_msd[np.where((M1[:,1] - M1_prev[:,1])<-900)[0],0] += side_len_domain
    M1_msd[np.where((M1[:,1] - M1_prev[:,1])> 900)[0],0] -= side_len_domain
    M1_msd[np.where((M1[:,2] - M1_prev[:,2])<-900)[0],1] += side_len_domain
    M1_msd[np.where((M1[:,2] - M1_prev[:,2])> 900)[0],1] -= side_len_domain
    
    M2_msd[np.where((M2[:,1] - M2_prev[:,1])<-900)[0],0] += side_len_domain
    M2_msd[np.where((M2[:,1] - M2_prev[:,1])> 900)[0],0] -= side_len_domain
    M2_msd[np.where((M2[:,2] - M2_prev[:,2])<-900)[0],1] += side_len_domain
    M2_msd[np.where((M2[:,2] - M2_prev[:,2])> 900)[0],1] -= side_len_domain
    
    return M1_msd, M2_msd

def RecordEcad(t, M1, M2, file_M1, file_M2, file_record_M1, file_record_M2, MSD = False):
    # row1: total par, row2: total_mono, row3: total_trans, row4: total_cis
    # row5: total stalled
    if MSD == False:
        M1_len          = len(M1)
        M1_all_mono     = np.count_nonzero(M1[:,4])
        M1_all_trans    = np.count_nonzero(M1[:,5]+M1[:,6])
        M1_all_cis      = sum((M1[:,7]!=-1)|(M1[:,8]!=-1))
        M1_all_stalled  = np.count_nonzero(M1[:,10])
        
        M2_len          = len(M2)
        M2_all_mono     = np.count_nonzero(M2[:,4])
        M2_all_trans    = np.count_nonzero(M2[:,5]+M2[:,6])
        M2_all_cis      = sum((M2[:,7]!=-1)|(M2[:,8]!=-1))
        M2_all_stalled  = np.count_nonzero(M2[:,10])
        
        string_M1 =  "    ".join(np.array([t, M1_len,M1_all_mono,
                                        M1_all_trans,M1_all_cis,
                                        M1_all_stalled]).astype(str))+"\n"
        string_M2 =  "    ".join(np.array([t, M2_len,M2_all_mono,
                                        M2_all_trans,M2_all_cis,
                                        M2_all_stalled]).astype(str))+"\n"
        
        file_record_M1.write(string_M1)
        file_record_M2.write(string_M2)

    header = np.zeros(M1.shape[1]); header[0] = len(M1)
    string_header = "    ".join((header).astype(str))+"\n"
    
    file_M1.write(string_header)
    file_M2.write(string_header)
    for i in range(len(M1)):
        string_M1 =  "    ".join(np.around(M1[i,:],decimals=2).astype(str))+"\n"
        file_M1.write(string_M1)
    for i in range(len(M2)):
        string_M2 =  "    ".join(np.around(M2[i,:],decimals=2).astype(str))+"\n"
        file_M2.write(string_M2)
# Actinseg1, Actinseg2, Actinfila1, Actinfila2 = Actinseg_all_1, Actinseg_all_2, ActinFila_all_1, ActinFila_all_2
def RecordActin(t, Actinseg1, Actinseg2, Actinfila1, Actinfila2, file1, file2):
    header = np.zeros(Actinfila1.shape[1]); header[0] = len(Actinfila1)
    string_header = "    ".join((header).astype(str))+"\n"
    
    file1.write(string_header)
    file2.write(string_header)
    for i in range(len(Actinfila1)):
        actin_array = np.r_[Actinfila1[i,0:3],Actinseg1[(i+1)*num_seg_fila-1,1:4]]
        string1 =  "    ".join(np.around(actin_array,decimals=7).astype(str))+"\n"
        file1.write(string1)
    for i in range(len(Actinfila2)):
        actin_array = np.r_[Actinfila2[i,0:3],Actinseg2[(i+1)*num_seg_fila-1,1:4]]
        string2 =  "    ".join(np.around(actin_array,decimals=7).astype(str))+"\n"
        file2.write(string2)
    
def TransCisLT(M1,M2,M1_trans_LT, M2_trans_LT, M1_cis_LT, M2_cis_LT):
    M1_trans_state = (M1[:,5] != -1) | (M1[:,6] != -1)
    M2_trans_state = (M2[:,5] != -1) | (M2[:,6] != -1)
    M1_trans_LT += M1_trans_state
    M2_trans_LT += M2_trans_state
    
    M1_cis_state = (M1[:,7] != -1) | (M1[:,8] != -1)
    M2_cis_state = (M2[:,7] != -1) | (M2[:,8] != -1)
    M1_cis_LT += M1_cis_state
    M2_cis_LT += M2_cis_state
    
    return M1_trans_LT, M2_trans_LT, M1_cis_LT, M2_cis_LT
    
def RecordTransCisLT(t, 
                     M1_trans_LT, M2_trans_LT, M1_cis_LT, M2_cis_LT,
                     file_M1_trans_LT, file_M2_trans_LT, file_M1_cis_LT, file_M2_cis_LT):
    file_M1_trans_LT.write(' '.join(list(map(str,M1_trans_LT)))+'\n')
    file_M2_trans_LT.write(' '.join(list(map(str,M2_trans_LT)))+'\n')  

    file_M1_cis_LT.write(' '.join(list(map(str,M1_cis_LT)))+'\n')
    file_M2_cis_LT.write(' '.join(list(map(str,M2_cis_LT)))+'\n')    

