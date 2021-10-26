# -*- coding: utf-8 -*-
"""
Created on Fri Mar  6 09:07:15 2020

@author: qiliny
"""
import newmodules as nw
import modules as mo
import ActinInitialization_psudo3d as AI
import random
import numpy as np
import math
import time
import sys
from inputset import *
import ctypes

#--------------------------------------------------------------------------
# Ecad Actin bond formation and dissociation
#--------------------------------------------------------------------------
# Ecad, Actinseg, EcadActin_pairs = M1, Actinseg_all_1, Ecad_Actin_binding1
def EcadActinAssDis(Ecad, Actinseg, EcadActin_pairs):
    # Find exiting Ecad Actin binding
    exist_EcadActin_pairs = Ecad[Ecad[:,11]!=-1,0].astype(int)
    
    #--------------- Association -------------
    if len(EcadActin_pairs)>0:
        
        # Association Probability (np array with random numbers)
        EcadActin_ass_prob = np.random.random((len(EcadActin_pairs),))
        
        # Successful binding
        ass_EcadActin = EcadActin_pairs[EcadActin_ass_prob<p_k_cadactin_a,:]
        
        # Perform Changes
        Ecad[ass_EcadActin[:,1],11] = ass_EcadActin[:,0]
        
    # -------------- Dissociation ------------- 
    
    if len(exist_EcadActin_pairs)>0:
        EcadActin_dis_prob = np.random.random((len(exist_EcadActin_pairs),))
        dis_EcadActin = exist_EcadActin_pairs[EcadActin_dis_prob<p_k_cadactin_d]
        # Perform Changes
        Ecad[dis_EcadActin,11] = -1
        
def ActinCrossBoundaryShift(Point1,Point2):
    if Point2.x-Point1.x > 3*collision_thre:
        Point2.x -= side_len_domain
    elif Point2.x-Point1.x < -3*collision_thre:
        Point2.x += side_len_domain
        
    if Point2.y-Point1.y > 3*collision_thre:
        Point2.y -= side_len_domain
    if Point2.y-Point1.y < -3*collision_thre:
        Point2.y += side_len_domain   

# Actin turnover
# Ecad, Actinseg, ActinFila = np.copy(M1), np.copy(Actinseg_all_1), np.copy(ActinFila_all_1)
def ActinTurnover(Ecad, Actinseg, ActinFila):
    actin_off = np.where(np.random.random(num_fila)<p_TO_r)[0]
    
    if len(actin_off) > 0:
        for actin_i in actin_off: 
            actin_start, actin_end = ActinFila[actin_i,4], ActinFila[actin_i,5]
            # unbind the Ecad/Actin
            Ecad[(Ecad[:,11] > actin_start-1) & (Ecad[:,11] < actin_end+1),11] = -1 
            
            # update the actin filament info
            ActinFila[actin_i,1:3] = np.random.rand(1,2)*side_len_domain
            ActinFila[actin_i,3] = np.random.rand()*2*np.pi
            
            Actinseg = AI.SegEachFila(Actinseg, ActinFila, actin_i)
            
    Actinseg = mo.adjustPBC(Actinseg)
    return Actinseg

def ShiftXY(x1,x2):
    
    if x2 > x1:
        x2 -= side_len_domain
    elif x2 < x1:
        x2 += side_len_domain 
    return x2

class Vector:
    def __init__(self, x, y): 
        self.x = x 
        self.y = y 

lib = ctypes.cdll.LoadLibrary('./libfun.so')
fun_coll = lib.collision
fun_coll.restype = ctypes.c_int
fun_coll.argtypes = [ctypes.c_double,ctypes.c_double,ctypes.c_double,
                      ctypes.c_double,ctypes.c_double]
#i, hash1, hash2, array1, array2 = n, nearby_part_desti, M1_all_hash, M1_mono_desti, M1_points
def CheckCollision(i, ori_i, nearby_part, array1, array2):
    # if any particle are in the nearby 9 partitions
    # Because length of step is about 1nm, original point must be in the nearby 9 partitions
    if len(nearby_part) > 1:
        nearby_part.remove(ori_i)
        x1 = array1[i,1]; y1 = array1[i,2];
        for j in nearby_part:
            
            x2 = array2[j,1]; y2 = array2[j,2] 
            collision = fun_coll(x1, y1, x2, y2, two_rad_ecad_sq)
            if collision:
                return 1
        return 0
    else:
        return 0   

# def CheckCollision(i, ori_i, nearby_part, array1, array2):
#     # if any particle are in the nearby 9 partitions
#     # Because length of step is about 1nm, original point must be in the nearby 9 partitions
#     if len(nearby_part) > 1:
#         nearby_part.remove(ori_i)
#         x1 = array1[i,1]; y1 = array1[i,2];
        
#         for j in nearby_part:

#             x2 = array2[j,1]; y2 = array2[j,2] 
#             distance = (x1-x2)**2+(y1-y2)**2
            
#             if distance < two_rad_ecad_sq:
#                 return 1    # Collision found  
#             # Check Cross boundary collision
#             if distance > collision_thre_sq*9:

#                 if abs(x1-x2)>3*collision_thre:
#                     x2 = ShiftXY(x1,x2)
#                 if abs(y1-y2)>3*collision_thre:    
#                     y2 = ShiftXY(y1,y2)

#                 distance = (x1-x2)**2+(y1-y2)**2
                
#                 if distance < two_rad_ecad_sq:
#                     return 1    
#         return 0            # no Collision
#     else:
#         return 0          


def TwoD_dist(x1,y1,x2,y2):
    x2,y2 = ShiftAcrossBoundary(x1,y1,x2,y2)
    return (x1-x2)**2+(y1-y2)**2

def actin_binding_condition(actin_binding, actin, ori_i, cad, intersection):
    # check 3D distance between cadherin to each actin node in the filament
    min_dist = 100000
    
    for i in range(len(actin)):
        dist_2D = TwoD_dist(cad.x, cad.y, actin[i,1],actin[i,2])
            
        if dist_2D <min_dist:
            min_dist = dist_2D
            actin_index = i
    if actin[actin_index,8] == 1:
        intersection = 0   # False intersect to the binding region
        actin_binding.append([int(actin[actin_index,0]), ori_i]) 
        
    return intersection
    
def ShiftAcrossBoundary(x1,y1,x2,y2):
    if abs(x1-x2)>3*collision_thre:
        x2 = ShiftXY(x1,x2)
    if abs(y1-y2)>3*collision_thre:    
        y2 = ShiftXY(y1,y2)

    return x2, y2

def CheckCrossActin(i, nearby_part, array1_ori, array1_desti, actinseg, actin_binding, collision):

    if len(nearby_part) > 0:
        ori_i = int(array1_ori[i,0])
        # Point, t0 position of Ecad
        cad1 = mo.Point(array1_ori[i,1],array1_ori[i,2])
        # Point, t1 position of Ecad
        cad2x, cad2y = ShiftAcrossBoundary(cad1.x,cad1.y,array1_desti[i,1],array1_desti[i,2])
        cad2 = mo.Point(cad2x, cad2y)
        
        # Sort actin segments into filaments, and do calculation only once 
        all_fila = set(actinseg[nearby_part,6].astype(int))
        all_seg_in_fila = dict()
        for fila in all_fila:
            all_seg_in_fila[int(fila)] = sorted([seg for seg in nearby_part if actinseg[seg,6] == fila])

        for fila in all_fila:

            posi_end = (fila+1)*num_seg_fila-1
            nega_end = fila*num_seg_fila
            
            start = all_seg_in_fila[fila][0]
            end = all_seg_in_fila[fila][-1]

            angle = actinseg[start,3]
            
            Actin1x, Actin1y = ShiftAcrossBoundary(cad1.x,cad1.y,actinseg[start,1],actinseg[start,2])
            Actin1 = mo.Point(Actin1x, Actin1y);
            Actin2x, Actin2y = ShiftAcrossBoundary(cad1.x,cad1.y,actinseg[end,1],actinseg[end,2])
            Actin2 = mo.Point(Actin2x, Actin2y);
            
            Actin3 = mo.Point(Actin1.x-len_seg*math.cos(angle),Actin1.y-len_seg*math.sin(angle)) 
            Actin4 = mo.Point(Actin2.x+len_seg*math.cos(angle),Actin2.y+len_seg*math.sin(angle)) 
            
            slope = (Actin4.y-Actin3.y)/(Actin4.x-Actin3.x)
            intersect = (Actin4.y-Actin4.x*slope)
            dist_ori   = abs(-slope*cad1.x+cad1.y-intersect)/math.sqrt(slope**2+1) # absolute distance, not distance square
            dist_desti = abs(-slope*cad2.x+cad2.y-intersect)/math.sqrt(slope**2+1)
            # First point is always the point found near the Ecad particle
            # First point is the reference for cross-boundary adjustment.
            if start == nega_end:
                
                intersection = mo.intersect(cad1,cad2,Actin1,Actin4)
                
                if not intersection and not collision and (dist_desti < EcadActin_bind_thre):
                    intersection = actin_binding_condition(actin_binding, actinseg[all_seg_in_fila[fila],:], ori_i, cad2, intersection)
                
                elif (intersection or collision) and (dist_ori < EcadActin_bind_thre):
                    intersection = actin_binding_condition(actin_binding, actinseg[all_seg_in_fila[fila],:], ori_i, cad2, intersection)
                    
                if intersection:
                    AB = Vector(Actin4.x-Actin1.x,Actin4.y-Actin1.y)
                    AC = Vector(cad2.x-Actin1.x,cad2.y-Actin1.y)
                    radian = (AB.x*AC.x+AB.y*AC.y)
                    if radian > 0:
                        return 1    # Movement is stopped by Actin filament
                    
            elif end == posi_end:
                             
                intersection = mo.intersect(cad1,cad2,Actin3,Actin2)
                
                if not intersection and not collision and (dist_desti < EcadActin_bind_thre):
                    intersection = actin_binding_condition(actin_binding, actinseg[all_seg_in_fila[fila],:], ori_i, cad2, intersection)
                
                elif (intersection or collision) and (dist_ori < EcadActin_bind_thre):
                    intersection = actin_binding_condition(actin_binding, actinseg[all_seg_in_fila[fila],:], ori_i, cad2, intersection)
                
                if intersection:
                    BA = Vector(Actin1.x-Actin4.x,Actin1.y-Actin4.y)
                    BC = Vector(cad2.x-Actin4.x,cad2.y-Actin4.y)
                    radian = (BA.x*BC.x+BA.y*BC.y)  
                    
                    if radian > 0:
                        return 1    # Movement is stopped by Actin filament
            
            # If either ends are not around
            else:
                    
                intersection = mo.intersect(cad1,cad2,Actin3,Actin4)

                if not intersection and not collision and (dist_desti < EcadActin_bind_thre):
                    intersection = actin_binding_condition(actin_binding, actinseg[all_seg_in_fila[fila],:], ori_i, cad2, intersection)
                elif (intersection or collision) and (dist_ori < EcadActin_bind_thre):
                    intersection = actin_binding_condition(actin_binding, actinseg[all_seg_in_fila[fila],:], ori_i, cad2, intersection)

                if intersection:
                    return 1    # Movement is stopped by Actin filament

        # if intersections were not found for all actin segment around the particle
        return 0    
    # if no actins are found in the surrounding partitions
    else:
        return 0        
    
def CheckNearbyPartition(hash1,hash2):         
    nearby_part = []
    for j in hash1:
        # if any actin are in the nearby 9 partitions
        if j in hash2:
            nearby_part += hash2[j]
    return nearby_part    

def PeriodicAngle(angle):
    if angle>np.pi*2:
        angle -= np.pi*2
    elif angle<0:
        angle += np.pi*2
        
    return angle

def CheckCollisionUpdatePos(i, ori_i,
                            cad_hash, cad_all_hash, Actin_hash,
                            cad_ori, cad_desti, cad_all, Actin):
    
    nearby_part_desti = nw.collision_hash[cad_hash[i]]      # nearby partition of the destination
    nearby_Ecad_desti = CheckNearbyPartition(nearby_part_desti, cad_all_hash)   # nearby cadherin
    
    nearby_Actin_ori = nw.EcadActin_hash[cad_hash[i]]
    nearby_Actin      = CheckNearbyPartition(nearby_Actin_ori, Actin_hash)     # nearby actin
    
    collision = CheckCollision(i, ori_i, nearby_Ecad_desti, cad_desti, cad_all)
    actincrossing = CheckCrossActin(i, nearby_Actin, cad_ori, cad_desti, Actin, [], collision)
    
    ecadcollision     = (collision or actincrossing)
    
    return ecadcollision

def CheckTransCollisionUpdatePos(i1, ori_i1,
                            cad_hash1, cad_all_hash1, Actin_hash1,
                            cad_ori1, cad_desti1, cad_all1, Actin1,
                            i2, ori_i2,
                            cad_hash2, cad_all_hash2, Actin_hash2,
                            cad_ori2, cad_desti2, cad_all2, Actin2):
    
    nearby_part_desti_1 = nw.collision_hash[cad_hash1[i1]]
    nearby_Ecad_desti_1 = CheckNearbyPartition(nearby_part_desti_1, cad_all_hash1) 
    
    nearby_part_desti_2 = nw.collision_hash[cad_hash2[i2]]    
    nearby_Ecad_desti_2 = CheckNearbyPartition(nearby_part_desti_2, cad_all_hash2)  
    
    nearby_Actin_ori_1 = nw.EcadActin_hash[cad_hash1[i1]]
    nearby_Actin_M1 = CheckNearbyPartition(nearby_Actin_ori_1, Actin_hash1) 
    
    nearby_Actin_ori_2 = nw.EcadActin_hash[cad_hash2[i2]]
    nearby_Actin_M2 = CheckNearbyPartition(nearby_Actin_ori_2, Actin_hash2)        
    
    collision1 = CheckCollision(i1, ori_i1, nearby_Ecad_desti_1, cad_desti1, cad_all1)
    collision2 = CheckCollision(i2, ori_i2, nearby_Ecad_desti_2, cad_desti2, cad_all2)
    collision = collision1 or collision2
    actincrossing1 = CheckCrossActin(i1, nearby_Actin_M1, cad_ori1, cad_desti1, Actin1, [], collision)
    actincrossing2 = CheckCrossActin(i2, nearby_Actin_M2, cad_ori2, cad_desti2, Actin2, [], collision)
    
    ecadcollision = (collision1 or collision2 or actincrossing1 or actincrossing2)
    
    return ecadcollision

def MobilityTransPartner(particle_array, cad_all):
    trans_partner   = np.copy(particle_array[:,9])
    trans_partner_mob = np.ones(len(trans_partner),dtype=bool)
    for i in range(len(trans_partner)):
        # this monomer doesn't have a trans partner
        if trans_partner[i] == -1:
            trans_partner_mob[i] = True
        # this monomer is not binding with actin
        elif cad_all[int(trans_partner[i]),11] ==-1:
            trans_partner_mob[i] = True
        
        else:
            trans_partner_mob[i] = False
    return trans_partner_mob

def StructuralAlignment(M1, M2, Actinseg_all_1, Actinseg_all_2, new_trans, new_cis_M1, new_cis_M2):
    # --------------1st part--------------
    # Get hesh tables for all particles in the system for collision check (due to volume exclusion) 
    # Hash of all particles, t0 position (KEY: Particle index, Value: Partition index)
    if (len(new_trans) > 0) or (len(new_cis_M1) > 0) or (len(new_cis_M2) > 0):

        M1_all_indexhash = mo.CreateChainingHashArray(M1,side_len_par_collision,num_part_side_collision)
        M2_all_indexhash = mo.CreateChainingHashArray(M2,side_len_par_collision,num_part_side_collision)
        
        M1_all_hash = nw.CreateHash(M1_all_indexhash)
        M2_all_hash = nw.CreateHash(M2_all_indexhash)
    
        Actin1_all_hash = nw.CreateHash(mo.CreateChainingHashArray(Actinseg_all_1,side_len_par_collision,num_part_side_collision))
        Actin2_all_hash = nw.CreateHash(mo.CreateChainingHashArray(Actinseg_all_2,side_len_par_collision,num_part_side_collision))
        
    # -------------- 2nd part Structural Alignment --------------
    
    # 2.1 trans alignment 

    if len(new_trans) > 0:
        M1_ori = M1[new_trans[:,0],:]
        M1_hash = mo.CreateChainingHashArray(M1_ori,side_len_par_collision,num_part_side_collision)
        M2_ori = M2[new_trans[:,1],:]
        M2_hash = mo.CreateChainingHashArray(M2_ori,side_len_par_collision,num_part_side_collision)
        
        mobile_M1 = ((M1[new_trans[:,0],10] == 0) & (M1[new_trans[:,0],11] == -1))
        mobile_M2 = ((M2[new_trans[:,1],10] == 0) & (M2[new_trans[:,1],11] == -1))

        for i in range(len(new_trans)):
            M1_i = new_trans[i,0]; M2_i = new_trans[i,1]
            M1_collision = CheckCollisionUpdatePos(i,M1_i,
                                                   M1_hash, M1_all_hash, Actin1_all_hash,
                                                   M1_ori, M2_ori, M1, Actinseg_all_1)
            M2_collision = CheckCollisionUpdatePos(i,M2_i,
                                                   M2_hash,M2_all_hash,Actin2_all_hash,
                                                   M2_ori, M1_ori,M2,Actinseg_all_2)
            if mobile_M1[i] and mobile_M2[i]:
                choice = random.choice([1,2])
                if choice == 1: 
                    if not M1_collision:
                        M1[M1_i,1:3] = M2_ori[i,1:3]; M1[M1_i,3] = PeriodicAngle(M2[M2_i,3]+np.pi/2)    
                    elif not M2_collision:
                        M2[M2_i,1:3] = M1_ori[i,1:3]; M2[M2_i,3] = PeriodicAngle(M1[M1_i,3]-np.pi/2) 
                else: 
                    if not M2_collision:
                        M2[M2_i,1:3] = M1_ori[i,1:3]; M2[M2_i,3] = PeriodicAngle(M1[M1_i,3]-np.pi/2)    
                    elif not M1_collision:
                        M1[M1_i,1:3] = M2_ori[i,1:3]; M1[M1_i,3] = PeriodicAngle(M2[M2_i,3]+np.pi/2) 
                        
            elif mobile_M1[i] and not M1_collision:
                M1[M1_i,1:3] = M2_ori[i,1:3]; M1[M1_i,3] = PeriodicAngle(M2[M2_i,3]+np.pi/2)
                
            elif mobile_M2[i] and not M2_collision:
                M2[M2_i,1:3] = M1_ori[i,1:3]; M2[M2_i,3] = PeriodicAngle(M1[M1_i,3]-np.pi/2)
        
    # 2.2 M1 Cis alignment 
    if len(new_cis_M1) > 0:
        donor_ori = M1[new_cis_M1[:,0],:]
        donor_hash = mo.CreateChainingHashArray(donor_ori,side_len_par_collision,num_part_side_collision)
        recep_ori = M1[new_cis_M1[:,1],:]
        recep_hash = mo.CreateChainingHashArray(recep_ori,side_len_par_collision,num_part_side_collision)
        
        donor_desti = np.copy(recep_ori)
        donor_desti[:,1] = recep_ori[:,1]-2*(rad_ecad+0.001)*np.cos(recep_ori[:,3])
        donor_desti[:,2] = recep_ori[:,2]-2*(rad_ecad+0.001)*np.sin(recep_ori[:,3])
        donor_desti_adjustPBC = mo.adjustPBC(donor_desti)
        
        recep_desti = np.copy(donor_ori)
        recep_desti[:,1] = donor_ori[:,1]+2*(rad_ecad+0.001)*np.cos(donor_ori[:,3])
        recep_desti[:,2] = donor_ori[:,2]+2*(rad_ecad+0.001)*np.sin(donor_ori[:,3])
        recep_desti_adjustPBC = mo.adjustPBC(recep_desti)
        
        donor_trans_mob = MobilityTransPartner(donor_ori,M2)
        recep_trans_mob = MobilityTransPartner(recep_ori,M2)
        
        # This monomer is not in another cis binding,
        # This monomer doesn't have a trans partner or 
        # The trans partner of this monomer is either not in a cis binding and a actin binding, 
        
        mobile_donor = ((donor_ori[:,10] == 0) & (donor_ori[:,11] == -1) & donor_trans_mob)
        mobile_recep = ((recep_ori[:,10] == 0) & (recep_ori[:,11] == -1) & recep_trans_mob)
        
        for i in range(len(new_cis_M1)):
            donor_i = new_cis_M1[i,0]; recep_i = new_cis_M1[i,1]
            donor_trans_collision = 0; recep_trans_collision = 0
            
            donor_trans_i   = int(donor_ori[i,9])
            if donor_trans_i == -1:
                donor_collision = CheckCollisionUpdatePos(i,donor_i,
                                                          donor_hash, M1_all_hash, Actin1_all_hash,
                                                          donor_ori,  donor_desti, M1, Actinseg_all_1)
            else:
                trans_ori = np.copy(donor_ori)
                trans_ori[i,:] = M2[donor_trans_i,:]
                donor_collision = CheckTransCollisionUpdatePos(i,donor_i,
                                                          donor_hash, M1_all_hash, Actin1_all_hash,
                                                          donor_ori,  donor_desti, M1, Actinseg_all_1,
                                                          i,donor_trans_i,
                                                          donor_hash, M2_all_hash, Actin2_all_hash,
                                                          trans_ori,  donor_desti, M2, Actinseg_all_2)
            
            recep_trans_i   = int(recep_ori[i,9])
            if recep_trans_i == -1:
                recep_collision = CheckCollisionUpdatePos(i,recep_i,
                                                          recep_hash, M1_all_hash, Actin1_all_hash,
                                                          recep_ori, recep_desti, M1, Actinseg_all_1) 
            else:
                trans_ori = np.copy(recep_ori)
                trans_ori[i,:] = M2[recep_trans_i,:]
                recep_collision = CheckTransCollisionUpdatePos(i,recep_i,
                                                          recep_hash, M1_all_hash, Actin1_all_hash,
                                                          recep_ori, recep_desti, M1, Actinseg_all_1,
                                                          i,recep_trans_i,
                                                          recep_hash, M2_all_hash, Actin2_all_hash,
                                                          trans_ori,  recep_desti, M2, Actinseg_all_2)

            # Level 1: check mobility
            if mobile_donor[i] and mobile_recep[i]:

                # Level 2: Check Collision
                if not donor_collision and not recep_collision:
                    choice = random.choice([1,2])
                    if choice == 1: 
                        
                        M1[donor_i,1:4] = donor_desti_adjustPBC[i,1:4]
                        # Level 4: Check if has trans partner
                        if donor_trans_i != -1:
                            M2[donor_trans_i,1:3] = donor_desti_adjustPBC[i,1:3]; M2[donor_trans_i,3] = PeriodicAngle(donor_desti[i,3]-np.pi/2) 
                    else: 

                        M1[recep_i,1:4] = recep_desti_adjustPBC[i,1:4]
                        # Level 4: Check if has trans partner
                        if recep_trans_i != -1:
                            M2[recep_trans_i,1:3] = recep_desti_adjustPBC[i,1:3]; M2[recep_trans_i,3] = PeriodicAngle(recep_desti[i,3]-np.pi/2) 
                
                elif not donor_collision:

                        M1[donor_i,1:4] = donor_desti_adjustPBC[i,1:4]
                        # Level 4: Check if has trans partner
                        if donor_trans_i != -1:
                            M2[donor_trans_i,1:3] = donor_desti_adjustPBC[i,1:3]; M2[donor_trans_i,3] = PeriodicAngle(donor_desti[i,3]-np.pi/2) 
                
                elif not recep_collision:  

                        M1[recep_i,1:4] = recep_desti_adjustPBC[i,1:4]
                        # Level 4: Check if has trans partner
                        if recep_trans_i != -1:
                            M2[recep_trans_i,1:3] = recep_desti_adjustPBC[i,1:3]; M2[recep_trans_i,3] = PeriodicAngle(recep_desti[i,3]-np.pi/2)     
                            
            elif mobile_donor[i]:
                if not donor_collision:
                    
                    M1[donor_i,1:4] = donor_desti_adjustPBC[i,1:4]
                    
                    # Level 4: Check if has trans partner
                    if donor_trans_i != -1:
                        M2[donor_trans_i,1:3] = donor_desti_adjustPBC[i,1:3]; M2[donor_trans_i,3] = PeriodicAngle(donor_desti[i,3]-np.pi/2) 
            
            elif mobile_recep[i]:
                if not recep_collision:   

                    M1[recep_i,1:4] = recep_desti_adjustPBC[i,1:4]
                    
                    # Level 4: Check if has trans partner
                    if recep_trans_i != -1:
                        M2[recep_trans_i,1:3] = recep_desti_adjustPBC[i,1:3]; M2[recep_trans_i,3] = PeriodicAngle(recep_desti[i,3]-np.pi/2)     
            
    # 2.3 M2 Cis alignment 
    if len(new_cis_M2) > 0:
        donor_ori = M2[new_cis_M2[:,0],:]
        donor_hash = mo.CreateChainingHashArray(donor_ori,side_len_par_collision,num_part_side_collision)
        recep_ori = M2[new_cis_M2[:,1],:]
        recep_hash = mo.CreateChainingHashArray(recep_ori,side_len_par_collision,num_part_side_collision)
        
        donor_desti = np.copy(recep_ori)
        donor_desti[:,1] = recep_ori[:,1]-2*(rad_ecad+0.001)*np.cos(recep_ori[:,3])
        donor_desti[:,2] = recep_ori[:,2]-2*(rad_ecad+0.001)*np.sin(recep_ori[:,3])
        donor_desti_adjustPBC = mo.adjustPBC(donor_desti)
        
        recep_desti = np.copy(donor_ori)
        recep_desti[:,1] = donor_ori[:,1]+2*(rad_ecad+0.001)*np.cos(donor_ori[:,3])
        recep_desti[:,2] = donor_ori[:,2]+2*(rad_ecad+0.001)*np.sin(donor_ori[:,3])
        recep_desti_adjustPBC = mo.adjustPBC(recep_desti)
        
        donor_trans_mob = MobilityTransPartner(donor_ori,M1)
        recep_trans_mob = MobilityTransPartner(recep_ori,M1)
        
        # This monomer is not in another cis binding,
        # This monomer doesn't have a trans partner or 
        # The trans partner of this monomer is either not in a cis binding and a actin binding, 
        
        mobile_donor = ((donor_ori[:,10] == 0) & (donor_ori[:,11] == -1) & donor_trans_mob)
        mobile_recep = ((recep_ori[:,10] == 0) & (recep_ori[:,11] == -1) & recep_trans_mob)
        
        for i in range(len(new_cis_M2)):
            donor_i = new_cis_M2[i,0]; recep_i = new_cis_M2[i,1]
            donor_trans_collision = 0; recep_trans_collision = 0
            
            donor_trans_i   = int(donor_ori[i,9])
            if donor_trans_i == -1:
                donor_collision = CheckCollisionUpdatePos(i,donor_i,
                                                          donor_hash, M2_all_hash, Actin2_all_hash,
                                                          donor_ori,  donor_desti, M2, Actinseg_all_2)
            else:
                trans_ori = np.copy(donor_ori)
                trans_ori[i,:] = M1[donor_trans_i,:]
                donor_collision = CheckTransCollisionUpdatePos(i,donor_i,
                                                          donor_hash, M2_all_hash, Actin2_all_hash,
                                                          donor_ori,  donor_desti, M2, Actinseg_all_2,
                                                          i,donor_trans_i,
                                                          donor_hash, M1_all_hash, Actin1_all_hash,
                                                          trans_ori,  donor_desti, M1, Actinseg_all_1)
            
            recep_trans_i   = int(recep_ori[i,9])
            if recep_trans_i == -1:
                recep_collision = CheckCollisionUpdatePos(i,recep_i,
                                                          recep_hash, M2_all_hash, Actin2_all_hash,
                                                          recep_ori, recep_desti, M2, Actinseg_all_2) 
            else:
                trans_ori = np.copy(recep_ori)
                trans_ori[i,:] = M1[recep_trans_i,:]
                recep_collision = CheckTransCollisionUpdatePos(i,recep_i,
                                                          recep_hash, M2_all_hash, Actin2_all_hash,
                                                          recep_ori, recep_desti, M2, Actinseg_all_2,
                                                          i,recep_trans_i,
                                                          recep_hash, M1_all_hash, Actin1_all_hash,
                                                          trans_ori,  recep_desti, M1, Actinseg_all_1)

            # Level 1: check mobility
            if mobile_donor[i] and mobile_recep[i]:
                
                # Level 2: Check Collision
                if (not donor_collision) and (not recep_collision):
                    choice = random.choice([1,2])
                    if choice == 1: 

                        M2[donor_i,1:4] = donor_desti_adjustPBC[i,1:4]
                        # Level 4: Check if has trans partner
                        if donor_trans_i != -1:
                            M1[donor_trans_i,1:3] = donor_desti_adjustPBC[i,1:3]; M1[donor_trans_i,3] = PeriodicAngle(donor_desti[i,3]-np.pi/2) 
                    else: 

                        M2[recep_i,1:4] = recep_desti_adjustPBC[i,1:4]
                        # Level 3: Check if has trans partner
                        if recep_trans_i != -1:
                            M1[recep_trans_i,1:3] = recep_desti_adjustPBC[i,1:3]; M1[recep_trans_i,3] = PeriodicAngle(recep_desti[i,3]-np.pi/2) 
                
                elif not donor_collision:
                    
                        M2[donor_i,1:4] = donor_desti_adjustPBC[i,1:4]
                        # Level 4: Check if has trans partner
                        if donor_trans_i != -1:
                            M1[donor_trans_i,1:3] = donor_desti_adjustPBC[i,1:3]; M1[donor_trans_i,3] = PeriodicAngle(donor_desti[i,3]-np.pi/2) 
                
                elif not recep_collision: 
                    
                        M2[recep_i,1:4] = recep_desti_adjustPBC[i,1:4]
                        # Level 4: Check if has trans partner
                        if recep_trans_i != -1:
                            M1[recep_trans_i,1:3] = recep_desti_adjustPBC[i,1:3]; M1[recep_trans_i,3] = PeriodicAngle(recep_desti[i,3]-np.pi/2)     
                            
            elif mobile_donor[i]:
                if not donor_collision:

                    M2[donor_i,1:4] = donor_desti_adjustPBC[i,1:4]
                    # Level 3: Check if has trans partner
                    if donor_trans_i != -1:
                        M1[donor_trans_i,1:3] = donor_desti_adjustPBC[i,1:3]; M1[donor_trans_i,3] = PeriodicAngle(donor_desti[i,3]-np.pi/2) 
            
            elif mobile_recep[i]:
                if not recep_collision:  

                    M2[recep_i,1:4] = recep_desti_adjustPBC[i,1:4]
                    # Level 3: Check if has trans partner
                    if recep_trans_i != -1:
                        M1[recep_trans_i,1:3] = recep_desti_adjustPBC[i,1:3]; M1[recep_trans_i,3] = PeriodicAngle(recep_desti[i,3]-np.pi/2)     
    
# Update immobile state 
def CisStallUpdate(M1,M2):   
    # All Ecad involved in Cis    
    all_cis_M1 = M1[(M1[:,7]!=-1)|(M1[:,8]!=-1),0].astype(int)
    all_cis_M2 = M2[(M2[:,7]!=-1)|(M2[:,8]!=-1),0].astype(int)
    
    # Stalled due to trans partner
    all_cis_M1_transpar = M1[all_cis_M1,9].astype(int)          # Stalled M2 due to its trans partner in M1
    all_cis_M2_transpar = M2[all_cis_M2,9].astype(int)          # Stalled M1 due to its trans partner in M2
    
    # All immobilized Ecad due to cis
    all_stall_M1 = np.append(all_cis_M1,all_cis_M2_transpar[all_cis_M2_transpar != -1])
    all_stall_M2 = np.append(all_cis_M2,all_cis_M1_transpar[all_cis_M1_transpar != -1])
    
    # Initial mobile state
    M1[:,10] = 0
    M2[:,10] = 0
    
    M1[all_stall_M1,10] = 1
    M2[all_stall_M2,10] = 1
    return M1, M2
    
def PositionUpdate(M1, M2, Actinseg_all_1, Actinseg_all_2):                            
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
    M1_mono_ori   = M1[(M1[:,4]==1)&(M1[:,10]==0)&(M1[:,11]==-1),0:4];
    M1_mono_desti = np.copy(M1_mono_ori); len_M1_m = len(M1_mono_desti)
    M1_move_ori = np.random.random(len_M1_m)*2*np.pi
    M1_mono_desti[:,1] += steplength*np.cos(M1_move_ori)
    M1_mono_desti[:,2] += steplength*np.sin(M1_move_ori)
    
    M2_mono_ori   = M2[(M2[:,4]==1)&(M2[:,10]==0)&(M2[:,11]==-1),0:4];
    M2_mono_desti = np.copy(M2_mono_ori); len_M2_m = len(M2_mono_desti)
    M2_move_ori = np.random.random(len_M2_m)*2*np.pi
    M2_mono_desti[:,1] += steplength*np.cos(M2_move_ori)
    M2_mono_desti[:,2] += steplength*np.sin(M2_move_ori)

    M1_trans_ori   = M1[((M1[:,5]==1)|(M1[:,6]==1))&(M1[:,10]==0)&(M1[:,11]==-1),:]
    M2_trans_index = M1_trans_ori[:,9].astype(int)
    M2_trans_movable = M2[M2_trans_index,11]==-1
    M1_trans_ori = M1_trans_ori[M2_trans_movable,:]
    
    M1_trans_desti = np.copy(M1_trans_ori)
    len_trans = len(M1_trans_desti)
    trans_move_ori = np.random.random(len_trans)*2*np.pi
    M1_trans_desti[:,1] += steplength*np.cos(trans_move_ori)
    M1_trans_desti[:,2] += steplength*np.sin(trans_move_ori)
    
    M2_trans_ori   = M2[M1[M1_trans_desti[:,0].astype(int),9].astype(int),0:4]
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
            nearby_Ecad_desti = CheckNearbyPartition(nearby_part_desti, M1_all_hash) 
            
            nearby_Actin_ori = nw.EcadActin_hash[M1_mono_ori_hash[n]]
            nearby_Actin_M1 = CheckNearbyPartition(nearby_Actin_ori, Actin1_all_hash)        
            
            collision = CheckCollision(n, M1_i, nearby_Ecad_desti, M1_mono_desti_adjustPBC, M1)
            actincrossing = CheckCrossActin(n, nearby_Actin_M1, M1_mono_ori, M1_mono_desti, Actinseg_all_1, actin_binding_M1, collision)
            
            # if a collision is found, move to the next particle, without updating position
            if (collision or actincrossing):
                continue

            # udpate position
            else:
                M1[M1_i,1:3] = M1_mono_desti_adjustPBC[n,1:3]   # Update position          

        # M2 particles Update
        elif (i > len_M1_m-1) & (i < (len_M1_m+len_M2_m)):
            n = i-len_M1_m
            M2_i = int(M2_mono_desti[n,0])                  # monomer index
                
            # Find Ecad in the surrounding partitions 
            nearby_part_desti = nw.collision_hash[M2_mono_desti_adjustPBC_hash[n]]
            nearby_Ecad_desti = CheckNearbyPartition(nearby_part_desti, M2_all_hash) 
            
            nearby_Actin_ori = nw.EcadActin_hash[M2_mono_ori_hash[n]]
            nearby_Actin_M2 = CheckNearbyPartition(nearby_Actin_ori, Actin2_all_hash)        
                  
            collision =  CheckCollision(n, M2_i, nearby_Ecad_desti, M2_mono_desti_adjustPBC, M2)
            actincrossing = CheckCrossActin(n, nearby_Actin_M2, M2_mono_ori, M2_mono_desti, Actinseg_all_2, actin_binding_M2, collision)
            # if a collision is found, move to the next particle, without updating position
            if (collision or actincrossing):
                continue

            # udpate position
            else:
                M2[M2_i,1:3] = M2_mono_desti_adjustPBC[n,1:3]   # Update position          

        # trans particles update
        else:
            n = i-len_M1_m-len_M2_m
            M1_i = int(M1_trans_desti[n,0])                     # monomer index
            M2_i = int(M2_trans_desti[n,0])                     # monomer index
            
            nearby_part_desti_1 = nw.collision_hash[M1_trans_desti_adjustPBC_hash[n]]
            nearby_Ecad_desti_1 = CheckNearbyPartition(nearby_part_desti_1, M1_all_hash) 
            
            nearby_part_desti_2 = nw.collision_hash[M2_trans_desti_adjustPBC_hash[n]]    
            nearby_Ecad_desti_2 = CheckNearbyPartition(nearby_part_desti_2, M2_all_hash)  
            
            nearby_Actin_ori_1 = nw.EcadActin_hash[M1_trans_ori_hash[n]]
            nearby_Actin_M1 = CheckNearbyPartition(nearby_Actin_ori_1, Actin1_all_hash) 
            
            nearby_Actin_ori_2 = nw.EcadActin_hash[M2_trans_ori_hash[n]]
            nearby_Actin_M2 = CheckNearbyPartition(nearby_Actin_ori_2, Actin2_all_hash)        
            
            collision1 = CheckCollision( n, M1_i, nearby_Ecad_desti_1, M1_trans_desti_adjustPBC, M1)
            collision2 = CheckCollision( n, M2_i, nearby_Ecad_desti_2, M2_trans_desti_adjustPBC, M2)
            collision = collision1 or collision2
            actincrossing1 = CheckCrossActin(n, nearby_Actin_M1, M1_trans_ori, M1_trans_desti, Actinseg_all_1, actin_binding_M1, collision)
            actincrossing2 = CheckCrossActin(n, nearby_Actin_M2, M2_trans_ori, M2_trans_desti, Actinseg_all_2, actin_binding_M2, collision)
                
            if (collision1 or collision2 or actincrossing1 or actincrossing2):
                continue

            else:
                #--------------------------------------------------------------
                # Update M1 trans position
                #--------------------------------------------------------------
                
                M1[M1_i,1:3] = M1_trans_desti_adjustPBC[n,1:3]      # Update position

                #--------------------------------------------------------------
                # Update M2 trans position
                #--------------------------------------------------------------

                M2[M2_i,1:3] = M2_trans_desti_adjustPBC[n,1:3]        

    return np.array(actin_binding_M1), np.array(actin_binding_M2)
