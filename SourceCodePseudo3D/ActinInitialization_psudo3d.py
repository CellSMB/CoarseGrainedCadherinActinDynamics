# -*- coding: utf-8 -*-
"""
Created on Sun Aug 30 19:27:40 2020

@author: Qilin Yu
"""

import math
import numpy as np
import modules as mo
from scipy.stats import truncnorm
from inputset import *

actin_binding_seg_index = np.arange(0,np.round(start_index),1,dtype=int)

def SegEachFila(Actinseg_all, ActinFila_all, i):
    
    for j in range(num_seg_fila):
        row = i*num_seg_fila+j
        Actinseg_all[row, 1] = ActinFila_all[i,1]+j*len_seg*np.cos(ActinFila_all[i,3])
        Actinseg_all[row, 2] = ActinFila_all[i,2]+j*len_seg*np.sin(ActinFila_all[i,3])
        Actinseg_all[row, 3] = ActinFila_all[i,3]
        
        if j == 0:
            Actinseg_all[row, 5] = row+1
        elif j == num_seg_fila-1:
            Actinseg_all[row, 4] = row-1
        else:
            Actinseg_all[row, 4] = row-1
            Actinseg_all[row, 5] = row+1
        Actinseg_all[row,6] = i
    
    Actinseg_all[i*num_seg_fila+actin_binding_seg_index,8] = 1
    
    return Actinseg_all

def Actinseg(ActinFila_all):    
    #------------------------------------------------------------------------------
    # Actinseg_all  Information of each actin segments (np.array)
    #------------------------------------------------------------------------------
    # Actin Network Initialization (Attributes of each actin segment)
    # Col 0: Index,             Col 4: Pointed end (index of next segment),         
    # Col 1: x_pos,             Col 5: Barbed end (index of previous segment),
    # Col 2: y_pos,             Col 6: Filament index
    # Col 3: Orientation,       Col 7: CCC binding or not
    #                           Col 8: binding site or not
   
    Actinseg_all = np.zeros((num_fila*num_seg_fila,9)) 
    Actinseg_all[:,4:] = -1
    Actinseg_all[:,0] = np.arange(num_fila*num_seg_fila)

    for i in range(num_fila):
        Actinseg_all = SegEachFila(Actinseg_all, ActinFila_all, i)
    
    Actinseg_all = mo.adjustPBC(Actinseg_all)
    
    return Actinseg_all

def InitActin():
    
    #------------------------------------------------------------------------------
    # ActinFila_all   Information of each actin filaments (np.array)
    #------------------------------------------------------------------------------
    # Col 0: Index,                                
    # Col 1: start point x_pos,         Col 4: Pointed end (index of small segments),               
    # Col 2: start point y_pos,         Col 5: Barbed end (index of small segments), 
    # Col 3: xy-plane Orientation, Phi
                  
    ActinFila_all = np.zeros((num_fila,6)) 
    ActinFila_all[:,0] = np.arange(num_fila)
    ActinFila_all[:,1:3] = np.random.rand(num_fila,2)*side_len_domain
    ActinFila_all[:,3] = np.random.rand(num_fila)*2*np.pi
    ActinFila_all[:,4] = np.arange(0,num_seg_fila*num_fila,num_seg_fila)
    ActinFila_all[:,5] = np.arange(0,num_seg_fila*num_fila,num_seg_fila)+num_seg_fila-1
    
    Actinseg_all = Actinseg(ActinFila_all)
    
    return ActinFila_all, Actinseg_all


ActinFila_all, Actinseg_all = InitActin()




