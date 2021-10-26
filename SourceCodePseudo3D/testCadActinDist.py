# -*- coding: utf-8 -*-
"""
Created on Sun Nov  8 12:17:50 2020

@author: Qilin Yu
"""
import newmodules as nw
import modules as mo
import ActinDynamics as AD
import random
import numpy as np
import math
import time
import sys
from inputset import *
import ctypes
import ActinCadInter as ACI
np.random.seed(1)
        
class Vector:
    def __init__(self, x, y): 
        self.x = x 
        self.y = y 
# i, nearby_part, array1_ori, array1_desti, array2, actin_binding = n, nearby_Actin_M1, M1_mono_ori, M1_mono_desti, Actinseg_all_1, actin_binding_M1
def CheckCrossActin(i, nearby_part, array1_ori, array1_desti, array2, actin_binding):

    if len(nearby_part) > 0:
        ori_i = int(array1_ori[i,0])
        # Point, t0 position of Ecad
        cad1 = mo.Point(array1_ori[i,1],array1_ori[i,2])
        # Point, t1 position of Ecad
        cad2 = mo.Point(array1_desti[i,1],array1_desti[i,2])
        
        # Sort actin segments into filaments, and do calculation only once 
        all_fila = set(array2[nearby_part,6].astype(int))
        all_seg_in_fila = dict()
        for fila in all_fila:
            all_seg_in_fila[int(fila)] = sorted([seg for seg in nearby_part if array2[seg,6] == fila])

        for fila in all_fila:
            if fila != array1_ori[i,11]:
                continue
            posi_end = (fila+1)*avg_len-1
            nega_end = fila*avg_len
            
            start = all_seg_in_fila[fila][0]
            end = all_seg_in_fila[fila][-1]

            angle = array2[start,3]
            
            Actin1x, Actin1y = ACI.ShiftAcrossBoundary(cad2.x,cad2.y,array2[start,1],array2[start,2])
            Actin1 = mo.Point(Actin1x, Actin1y);
            Actin2x, Actin2y = ACI.ShiftAcrossBoundary(cad2.x,cad2.y,array2[end,1],array2[end,2])
            Actin2 = mo.Point(Actin2x, Actin2y);
            
            Actin3 = mo.Point(Actin1.x-len_seg*math.cos(angle),Actin1.y-len_seg*math.sin(angle)) 
            Actin4 = mo.Point(Actin2.x+len_seg*math.cos(angle),Actin2.y+len_seg*math.sin(angle)) 
            
            slope = (Actin4.y-Actin3.y)/(Actin4.x-Actin3.x)
            intersect = (Actin4.y-Actin4.x*slope)
            dist = abs(-slope*cad1.x+cad1.y-intersect)/math.sqrt(slope**2+1)
            # First point is always the point found near the Ecad particle
            # First point is the reference for cross-boundary adjustment.
            if start == nega_end:
                AB = Vector(Actin4.x-Actin1.x,Actin4.y-Actin1.y)
                AC = Vector(cad2.x-Actin1.x,cad2.y-Actin1.y)
                radian = (AB.x*AC.x+AB.y*AC.y)
                
                if radian < 0:
                    return 0
                elif (radian > 0):
                   
                    if dist > EcadActin_bind_thre:
                        print("fila: {}, ori_i: {} ".format(fila, ori_i))
                        print(dist)
                        sys.exit()
                        
                    intersection = mo.intersect(cad1,cad2,Actin1,Actin4)
                    if intersection:
                        return 1    # Movement is stopped by Actin filament
    
                
            elif end == posi_end:
                BA = Vector(Actin1.x-Actin4.x,Actin1.y-Actin4.y)
                BC = Vector(cad2.x-Actin4.x,cad2.y-Actin4.y)
                radian = (BA.x*BC.x+BA.y*BC.y)                  
                
                if radian < 0:
                    return 0
                elif (radian > 0):
                    if dist > EcadActin_bind_thre:
                        print("fila: {}, ori_i: {} ".format(fila, ori_i))
                        print(dist)
                        sys.exit()
                        
                    intersection = mo.intersect(cad1,cad2,Actin3,Actin2)
                    
                    if intersection:
                        return 1    # Movement is stopped by Actin filament
                    

            else:
                
                if dist > EcadActin_bind_thre:
                    print("fila: {}, ori_i: {} ".format(fila, ori_i))
                    print(dist)
                    sys.exit()
                intersection = mo.intersect(cad1,cad2,Actin3,Actin4)

                if intersection:
                    return 1    # Movement is stopped by Actin filament
        # if intersections were found for all actin segment around the particle
        return 0    
    # if no actins are found in the surrounding partitions
    else:
        return 0        


def CheckError(M1, M2, Actinseg_all_1, Actinseg_all_2):     
    # --------------3st part--------------
    # update the hash after structural alignment
    M1_all_indexhash = mo.CreateChainingHashArray(M1,side_len_par_collision,num_part_side_collision)
    M2_all_indexhash = mo.CreateChainingHashArray(M2,side_len_par_collision,num_part_side_collision)
    
    M1_all_hash = nw.CreateHash(M1_all_indexhash)
    M2_all_hash = nw.CreateHash(M2_all_indexhash)
    
    Actin1_all_hash = nw.CreateHash(mo.CreateChainingHashArray(Actinseg_all_1,side_len_par_collision,num_part_side_collision))
    Actin2_all_hash = nw.CreateHash(mo.CreateChainingHashArray(Actinseg_all_2,side_len_par_collision,num_part_side_collision))
    
    # --------------4st part--------------
    # Get all the mobile E-cad and calculate their position in next step
    # st_9 = time.time()
    # Brownian dynamics algorithm 2 : Fixed length
    # Destination of movable partiples
    M1_mono_ori   = M1[(M1[:,4]==1)&(M1[:,10]==0)&(M1[:,11]!=-1),:];
    # M1_mono_ori   = M1[succ_binding,:]
    M1_mono_desti = np.copy(M1_mono_ori); len_M1_m = len(M1_mono_desti)
    M1_move_ori = np.random.random(len_M1_m)*2*np.pi
    M1_mono_desti[:,1] += steplength*np.cos(M1_move_ori)
    M1_mono_desti[:,2] += steplength*np.sin(M1_move_ori)
    
    M2_mono_ori   = M2[(M2[:,4]==1)&(M2[:,10]==0)&(M2[:,11]!=-1),:];
    M2_mono_desti = np.copy(M2_mono_ori); len_M2_m = len(M2_mono_desti)
    M2_move_ori = np.random.random(len_M2_m)*2*np.pi
    M2_mono_desti[:,1] += steplength*np.cos(M2_move_ori)
    M2_mono_desti[:,2] += steplength*np.sin(M2_move_ori)
    
    M1_trans_ori   = M1[((M1[:,5]==1)|(M1[:,6]==1))&(M1[:,10]==0)&(M1[:,11]!=-1),:]
    M1_trans_desti = np.copy(M1_trans_ori)
    len_trans = len(M1_trans_desti)
    trans_move_ori = np.random.random(len_trans)*2*np.pi
    M1_trans_desti[:,1] += steplength*np.cos(trans_move_ori)
    M1_trans_desti[:,2] += steplength*np.sin(trans_move_ori)
    
    M2_trans_ori   = M2[M1[M1_trans_desti[:,0].astype(int),9].astype(int),:]
    M2_trans_desti = np.copy(M2_trans_ori)
    M2_trans_desti[:,1] += steplength*np.cos(trans_move_ori)
    M2_trans_desti[:,2] += steplength*np.sin(trans_move_ori)
    
    # --------------5st part Determine Desti Partition--------------
    # Get hesh tables for all sets of mobile particles 
    # Hash of all particles, t0 position (KEY: Particle index, Value: Partition index)
    
    M1_mono_desti_adjustPBC = mo.adjustPBC(M1_mono_desti)
    M2_mono_desti_adjustPBC = mo.adjustPBC(M2_mono_desti)
    M1_trans_desti_adjustPBC = mo.adjustPBC(M1_trans_desti)
    M2_trans_desti_adjustPBC = mo.adjustPBC(M2_trans_desti)
    
    M1_mono_ori_hash = mo.CreateChainingHashArray(M1_mono_ori,side_len_par_collision,num_part_side_collision)
    M2_mono_ori_hash = mo.CreateChainingHashArray(M2_mono_ori,side_len_par_collision,num_part_side_collision)
    M1_trans_ori_hash = mo.CreateChainingHashArray(M1_trans_ori,side_len_par_collision,num_part_side_collision)
    M2_trans_ori_hash = mo.CreateChainingHashArray(M2_trans_ori,side_len_par_collision,num_part_side_collision)
    
    M1_mono_desti_adjustPBC_hash = mo.CreateChainingHashArray(M1_mono_desti_adjustPBC,side_len_par_collision,num_part_side_collision)
    M2_mono_desti_adjustPBC_hash = mo.CreateChainingHashArray(M2_mono_desti_adjustPBC,side_len_par_collision,num_part_side_collision)
    M1_trans_desti_adjustPBC_hash = mo.CreateChainingHashArray(M1_trans_desti_adjustPBC,side_len_par_collision,num_part_side_collision)
    M2_trans_desti_adjustPBC_hash = mo.CreateChainingHashArray(M2_trans_desti_adjustPBC,side_len_par_collision,num_part_side_collision)
    
    # Order of position update
    Move_order = np.random.permutation(np.arange(len_M1_m+len_M2_m+len_trans))
    
    # -------------- 6th part Position Update --------------
    # Iterate through all movable particles 
    # Volumn exclusion and Actin limitation lead to Movement Rejection 
    
    actin_binding_M1 = []
    actin_binding_M2 = []
    for i in Move_order:
    
        # M1 particles update
        if i < len_M1_m:
            
            n = i
            M1_i = int(M1_mono_desti[n,0])                  # monomer index
        
            # Find Ecad in the surrounding partitions 
            nearby_part_desti = nw.collision_hash[M1_mono_desti_adjustPBC_hash[n]]
            nearby_Ecad_desti = ACI.CheckNearbyPartition(nearby_part_desti, M1_all_hash) 
            
            nearby_Actin_ori = nw.EcadActin_hash[M1_mono_ori_hash[n]]
            nearby_Actin_M1 = ACI.CheckNearbyPartition(nearby_Actin_ori, Actin1_all_hash)        
            
            # if a collision is found, move to the next particle, without updating position
            if (CheckCrossActin(n, nearby_Actin_M1, M1_mono_ori, M1_mono_desti, Actinseg_all_1, actin_binding_M1) or
                ACI.CheckCollision(n, M1_i, nearby_Ecad_desti, M1_mono_desti_adjustPBC, M1)):
                continue
        
        # M2 particles Update
        elif (i > len_M1_m-1) & (i < (len_M1_m+len_M2_m)):
            n = i-len_M1_m
            M2_i = int(M2_mono_desti[n,0])                  # monomer index
                
            # Find Ecad in the surrounding partitions 
            nearby_part_desti = nw.collision_hash[M2_mono_desti_adjustPBC_hash[n]]
            nearby_Ecad_desti = ACI.CheckNearbyPartition(nearby_part_desti, M2_all_hash) 
            
            nearby_Actin_ori = nw.EcadActin_hash[M2_mono_ori_hash[n]]
            nearby_Actin_M2 = ACI.CheckNearbyPartition(nearby_Actin_ori, Actin2_all_hash)        
                  
            
            # if a collision is found, move to the next particle, without updating position
            if (CheckCrossActin(n, nearby_Actin_M2, M2_mono_ori, M2_mono_desti, Actinseg_all_2, actin_binding_M2) or
                ACI.CheckCollision(n, M2_i, nearby_Ecad_desti, M2_mono_desti_adjustPBC, M2)):
                continue
        
        # trans particles update
        else:
            n = i-len_M1_m-len_M2_m
            M1_i = int(M1_trans_desti[n,0])                     # monomer index
            M2_i = int(M2_trans_desti[n,0])                     # monomer index
            
            nearby_part_desti_1 = nw.collision_hash[M1_trans_desti_adjustPBC_hash[n]]
            nearby_Ecad_desti_1 = ACI.CheckNearbyPartition(nearby_part_desti_1, M1_all_hash) 
            
            nearby_part_desti_2 = nw.collision_hash[M2_trans_desti_adjustPBC_hash[n]]    
            nearby_Ecad_desti_2 = ACI.CheckNearbyPartition(nearby_part_desti_2, M2_all_hash)  
            
            nearby_Actin_ori_1 = nw.EcadActin_hash[M1_trans_ori_hash[n]]
            nearby_Actin_M1 = ACI.CheckNearbyPartition(nearby_Actin_ori_1, Actin1_all_hash) 
            
            nearby_Actin_ori_2 = nw.EcadActin_hash[M2_trans_ori_hash[n]]
            nearby_Actin_M2 = ACI.CheckNearbyPartition(nearby_Actin_ori_2, Actin2_all_hash)        
            
            actincrossing1 = CheckCrossActin(n, nearby_Actin_M1, M1_trans_ori, M1_trans_desti, Actinseg_all_1, actin_binding_M1)
            actincrossing2 = CheckCrossActin(n, nearby_Actin_M2, M2_trans_ori, M2_trans_desti, Actinseg_all_2, actin_binding_M2)
                
            if (actincrossing1 or actincrossing2 or 
                ACI.CheckCollision( n, M1_i, nearby_Ecad_desti_1, M1_trans_desti_adjustPBC, M1) or 
                ACI.CheckCollision( n, M2_i, nearby_Ecad_desti_2, M2_trans_desti_adjustPBC, M2)):
                continue

