# -*- coding: utf-8 -*-
"""
Created on Thu Mar 19 14:29:54 2020

@author: Qilin Yu
"""
from inputset import *
import numpy as np
import time
import ctypes

def HashRetrivalTable(num_part_side_len):

    hash_table = dict()
    
    for i in range(num_part_side_len**2):
        
        # Most Bot
        if i < num_part_side_len:
            # Most left Bot
            if i == 0: 
                top = [2*num_part_side_len-1,i+num_part_side_len,i+num_part_side_len+1]
                mid = [num_part_side_len-1,i,i+1]
                bot = [num_part_side_len**2-1,
                       num_part_side_len*(num_part_side_len-1),
                       num_part_side_len*(num_part_side_len-1)+1]
            # Most right Bot
            elif i == num_part_side_len-1:
                top = [i+num_part_side_len-1,i+num_part_side_len,i+1]
                mid = [i-1,i,0]
                bot = [num_part_side_len**2-2,
                       num_part_side_len**2-1,
                       num_part_side_len*(num_part_side_len-1)]
                
            else:
                top = [i+num_part_side_len-1,
                       i+num_part_side_len,
                       i+num_part_side_len+1]
                mid = [i-1,i,i+1]
                bot = [i+num_part_side_len*(num_part_side_len-1)-1,
                       i+num_part_side_len*(num_part_side_len-1),
                       i+num_part_side_len*(num_part_side_len-1)+1]
            hash_table[i] = top+mid+bot
            
        # Most Top
        elif i > num_part_side_len**2-num_part_side_len-1:
            # Most left TOP
            if i == num_part_side_len**2-num_part_side_len: 
                top = [num_part_side_len-1,0,1]
                mid = [num_part_side_len**2-1,i,i+1]
                bot = [i-1,i-num_part_side_len,i-num_part_side_len+1]
            # Most right TOP
            elif i == num_part_side_len**2-1:
                top = [i-num_part_side_len*(num_part_side_len-1)-1,
                       i-num_part_side_len*(num_part_side_len-1),
                       0]
                mid = [i-1,i,i-num_part_side_len+1]
                bot = [i-num_part_side_len-1,i-num_part_side_len,i-2*num_part_side_len+1]
            else:
                top = [i-num_part_side_len*(num_part_side_len-1)-1,
                       i-num_part_side_len*(num_part_side_len-1),
                       i-num_part_side_len*(num_part_side_len-1)+1]
                mid = [i-1,i,i+1]
                bot = [i-num_part_side_len-1,i-num_part_side_len,i-num_part_side_len+1]
            hash_table[i] = top+mid+bot
        # Most left
        elif i%num_part_side_len == 0:
            top = [i+num_part_side_len*2-1,i+num_part_side_len,i+num_part_side_len+1]
            mid = [i+num_part_side_len-1,i,i+1]
            bot = [i-1,i-num_part_side_len,i-num_part_side_len+1]
            hash_table[i] = top+mid+bot
        # Most right
        elif i%num_part_side_len == (num_part_side_len-1):
            top = [i+num_part_side_len-1,i+num_part_side_len,i+1]
            mid = [i-1,i,i-num_part_side_len+1]
            bot = [i-num_part_side_len-1,i-num_part_side_len,i-2*num_part_side_len+1]
            hash_table[i] = top+mid+bot
        
        # In the middle
        else:
            top = [i+num_part_side_len-1,i+num_part_side_len,i+num_part_side_len+1]
            mid = [i-1,i,i+1]
            bot = [i-num_part_side_len-1,i-num_part_side_len,i-num_part_side_len+1]
            hash_table[i] = top+mid+bot 
            
    return hash_table

def HashRetrivalTable_CadActin(num_part_side_len):

    hash_table = dict()
    
    for i in range(num_part_side_len**2):
        
        # Most Bot
        if i < num_part_side_len:
            # Most left Bot
            if i == 0: 
                top = [i+num_part_side_len,i+num_part_side_len+1]
                mid = [i,i+1]
                bot = []
            # Most right Bot
            elif i == num_part_side_len-1:
                top = [i+num_part_side_len-1,i+num_part_side_len]
                mid = [i-1,i]
                bot = []
            else:
                top = [i+num_part_side_len-1,
                       i+num_part_side_len,
                       i+num_part_side_len+1]
                mid = [i-1,i,i+1]
                bot = []
            hash_table[i] = top+mid+bot
            
        # Most Top
        elif i > num_part_side_len**2-num_part_side_len-1:
            # Most left TOP
            if i == num_part_side_len**2-num_part_side_len: 
                top = []
                mid = [i,i+1]
                bot = [i-num_part_side_len,i-num_part_side_len+1]
            # Most right TOP
            elif i == num_part_side_len**2-1:
                top = []
                mid = [i-1,i]
                bot = [i-num_part_side_len-1,i-num_part_side_len]
            else:
                top = []
                mid = [i-1,i,i+1]
                bot = [i-num_part_side_len-1,i-num_part_side_len,i-num_part_side_len+1]
            hash_table[i] = top+mid+bot
        # Most left
        elif i%num_part_side_len == 0:
            top = [i+num_part_side_len,i+num_part_side_len+1]
            mid = [i,i+1]
            bot = [i-num_part_side_len,i-num_part_side_len+1]
            hash_table[i] = top+mid+bot
        # Most right
        elif i%num_part_side_len == (num_part_side_len-1):
            top = [i+num_part_side_len-1,i+num_part_side_len]
            mid = [i-1,i]
            bot = [i-num_part_side_len-1,i-num_part_side_len]
            hash_table[i] = top+mid+bot
        
        # In the middle
        else:
            top = [i+num_part_side_len-1,i+num_part_side_len,i+num_part_side_len+1]
            mid = [i-1,i,i+1]
            bot = [i-num_part_side_len-1,i-num_part_side_len,i-num_part_side_len+1]
            hash_table[i] = top+mid+bot 
            
    return hash_table

def CreateHash(partition):
    hash_part = dict()
    
    for i in range(len(partition)):
        
        if not partition[i] in hash_part:
            hash_part[partition[i]] = [i]
        else: 
            hash_part[partition[i]].append(i)
            
    return hash_part

def ShiftXY(x1,x2):
    
    if x2 > x1:
        x2 -= side_len_domain
    elif x2 < x1:
        x2 += side_len_domain 
    return x2

lib = ctypes.cdll.LoadLibrary('./libfun.so')
fun_other = lib.trueorfalse
fun_trans = lib.trueorfalsetrans
fun_other.restype = ctypes.c_int
fun_trans.restype = ctypes.c_int
fun_other.argtypes = [ctypes.c_double,ctypes.c_double,ctypes.c_double,
                      ctypes.c_double,ctypes.c_double,ctypes.c_double,
                      ctypes.c_double,ctypes.c_double,ctypes.c_double]
fun_trans.argtypes = [ctypes.c_double,ctypes.c_double,ctypes.c_double,
                      ctypes.c_double,ctypes.c_double,ctypes.c_double,
                      ctypes.c_double,ctypes.c_double,ctypes.c_double]    
def FindBindingPairs(hash_1, hash_2, array1, array2, hash_retrival, dist_thre, trans = False):

    if trans:
        fun = fun_trans
    else:
        fun = fun_other
    true_pairs = []
    hash_2 = CreateHash(hash_2)
    
    for i in range(len(hash_1)):
        hash_retr_temp = hash_retrival[hash_1[i]]
        nearby_part = []
        for j in hash_retr_temp:
            if j in hash_2:
                nearby_part += hash_2[j]

        if len(nearby_part) > 0:
            x1 = array1[i,1]; y1 = array1[i,2];

            for j in nearby_part:
                if trans == False and array2[j,0] == array1[i,0]:
                    continue

                x2 = array2[j,1]; y2 = array2[j,2]; 
                angle1 = array1[i,3]; angle2 = array2[j,3];

                
                trueorfalse = fun(x1, y1, x2, y2,
                                  angle1, angle2,
                                  dist_thre, angle_thre_s, angle_thre_l)
                if trueorfalse:
                    true_pairs.append([int(array1[i,0]),int(array2[j,0])])
    
    return true_pairs

# def FindBindingPairs(hash_1, hash_2, array1, array2, hash_retrival, dist_thre, trans = False):
    
#     true_pairs = []                                 # Initial the container
    
#     hash_2 = CreateHash(hash_2)                     # Create a dict
    
#     for i in range(len(hash_1)):
#         hash_retr_temp = hash_retrival[hash_1[i]]
#         nearby_part = []
#         for j in hash_retr_temp:
#             if j in hash_2:
#                 nearby_part += hash_2[j]
        
#         # if any particle are in the nearby 9 partitions
#         if len(nearby_part) > 0:
#             for j in nearby_part:
                
#                 angle1 = array1[i,3]; angle2 = array2[j,3];
#                 judge_angle = False
#                 # Judge if the angle meets the binding requirement
#                 if trans:
#                     a0 = abs(angle1-(angle2+np.pi/2))
#                     a360 = abs(abs(angle1-(angle2+np.pi/2))-2*np.pi)
#                     if (a0 < angle_thre_s) or (a360 < angle_thre_s):
#                         judge_angle = True
#                 else:
#                     a0 = abs(angle1-angle2)
#                     if (a0 < angle_thre_s) or (a0 > angle_thre_l):
#                         judge_angle = True
                        
#                 # Judge if the dist meets the binding requirement    
#                 if judge_angle:
                    
#                     x1 = array1[i,1]; y1 = array1[i,2];
                
#                     x2 = array2[j,1]; y2 = array2[j,2];
                    
#                     distance = (x1-x2)**2+(y1-y2)**2
                    
# #                    print("x1: {}, x2: {}, y1: {}, y2: {}, dist: {}".format(x1,x2,y1,y2,distance))
                    
#                     if distance < dist_thre:
#                         true_pairs.append([int(array1[i,0]),int(array2[j,0])])
#                     # Check Cross boundary collision
#                     if distance > dist_thre*9:
#                         if abs(x1-x2)>900:
#                             x2 = ShiftXY(x1,x2)
#                         if abs(y1-y2)>900:
#                             y2 = ShiftXY(y1,y2)
#                         distance = (x1-x2)**2+(y1-y2)**2
                        
#                         if distance < dist_thre:
#                             true_pairs.append([int(array1[i,0]),int(array2[j,0])])       
            
#     return true_pairs

  
def CreateAllHash():
    trans_hash = HashRetrivalTable(num_part_side_trans)
    cis_hash = HashRetrivalTable(num_part_side_cis)
    collision_hash = HashRetrivalTable(num_part_side_collision)
    EcadActin_hash = HashRetrivalTable(num_part_side_collision)
    return trans_hash, cis_hash, EcadActin_hash, collision_hash

trans_hash, cis_hash, EcadActin_hash, collision_hash = CreateAllHash()















