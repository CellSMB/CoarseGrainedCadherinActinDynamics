# -*- coding: utf-8 -*-
"""
Created on Fri Feb 28 12:50:06 2020

@author: qiliny
"""
import time
import sys
import numpy as np
import ActinInitialization_psudo3d as AI
import CadDynamics as CD
import ActinCadInter as ACI
from modules import RecordEcad, RecordActin, RecordTransCisLT, TransCisLT, MSD_data
from inputset import *

# Ecad and Actin network initialization
ActinFila_all_1, Actinseg_all_1 = AI.InitActin()
ActinFila_all_2, Actinseg_all_2 = AI.InitActin()

M1 = CD.InitCad()
M2 = CD.InitCad()
M1_MSD, M2_MSD = np.zeros((len(M1),2)),np.zeros((len(M2),2))
M1_prev, M2_prev = np.copy(M1),np.copy(M2)

M1_cis_LT = np.zeros((num_par))
M2_cis_LT = np.zeros((num_par))
M1_trans_LT = np.zeros((num_par))
M2_trans_LT = np.zeros((num_par))

# Recording
file_record_M1 = open("CadInfo_M1.txt","w")
file_record_M2 = open("CadInfo_M2.txt","w")
file_M1 = open("M1Info.txt","w")
file_M2 = open("M2Info.txt","w")
file_M1_MSD = open("M1Info_MSD.txt","w")
file_M2_MSD = open("M2Info_MSD.txt","w")
file_Actin1 = open("Actin1.txt","w")
file_Actin2 = open("Actin2.txt","w")

file_record_M1_s = open("CadInfo_M1_s.txt","w")
file_record_M2_s = open("CadInfo_M2_s.txt","w")
file_M1_s = open("M1Info_s.txt","w")
file_M2_s = open("M2Info_s.txt","w")
file_M1_MSD_s = open("M1Info_MSD_s.txt","w")
file_M2_MSD_s = open("M2Info_MSD_s.txt","w")

file_M1_cis_LT = open("M1Info_cis_LT.txt","w")
file_M2_cis_LT = open("M2Info_cis_LT.txt","w")
file_M1_trans_LT = open("M1Info_trans_LT.txt","w")
file_M2_trans_LT = open("M2Info_trans_LT.txt","w")

start_t = time.time()

for t in range(Total_steps):

    # -------------- Step 1 E-cad trans & cis interactions   ------------------
    trans_pairs                 = CD.TransInteraction(M1,M2)
    M1_cis_pairs, M2_cis_pairs  = CD.CisInteraction(M1,M2)
    
    # -------------- Step 2 States Update                    ------------------
    M1, M2, new_trans, new_cis_M1, new_cis_M2 = CD.UpdateStates(trans_pairs, M1_cis_pairs, M2_cis_pairs, M1, M2)

    # -------------- Step 3 StructuralAlignment & State & Position & orientation update
    ACI.StructuralAlignment(M1, M2, Actinseg_all_1, Actinseg_all_2, new_trans, new_cis_M1, new_cis_M2)
    
    M1, M2 = ACI.CisStallUpdate(M1,M2)
    M1, M2 = CD.OriUpdate(M1, M2)
    
    Ecad_Actin_binding1,Ecad_Actin_binding2  = ACI.PositionUpdate(M1, M2, Actinseg_all_1, Actinseg_all_2)

    # -------------- Step 4 Actin turnover                   ------------------
    if ActinTurnover_toggle:
        Actinseg_all_1 = ACI.ActinTurnover(M1, Actinseg_all_1, ActinFila_all_1)
        Actinseg_all_2 = ACI.ActinTurnover(M2, Actinseg_all_2, ActinFila_all_2)
    
    # -------------- Step 5 Actin association & Dissociation ------------------
    if ActinBinding_toggle:
        ACI.EcadActinAssDis(M1, Actinseg_all_1, Ecad_Actin_binding1)
        ACI.EcadActinAssDis(M2, Actinseg_all_2, Ecad_Actin_binding2)
   
    # -------------- Step 6 Recording & Cis lifetime & MSD   ------------------
    M1_trans_LT, M2_trans_LT, M1_cis_LT, M2_cis_LT = TransCisLT(M1,M2,M1_trans_LT, M2_trans_LT, M1_cis_LT, M2_cis_LT)
    M1_MSD, M2_MSD = MSD_data(M1, M2, M1_prev, M2_prev, M1_MSD, M2_MSD)

    if t%Record_freq == 0:
        # Original Data
        RecordEcad(t, M1, M2, file_M1, file_M2, file_record_M1, file_record_M2)
        # MSD data
        M1_MSD_temp, M2_MSD_temp = np.copy(M1),np.copy(M2)
        M1_MSD_temp[:,1:3] += M1_MSD; M2_MSD_temp[:,1:3] += M2_MSD
        RecordEcad(t, M1_MSD_temp, M2_MSD_temp, file_M1_MSD, file_M2_MSD, file_record_M1, file_record_M2, MSD = True)
        # Actin data
        RecordActin(t, Actinseg_all_1, Actinseg_all_2, ActinFila_all_1, ActinFila_all_2, file_Actin1, file_Actin2)
        # Trans Cis Lifetime
        RecordTransCisLT(t, M1_trans_LT, M2_trans_LT, M1_cis_LT, M2_cis_LT,
                         file_M1_trans_LT, file_M2_trans_LT, file_M1_cis_LT, file_M2_cis_LT)
        M1_trans_LT, M2_trans_LT, M1_cis_LT, M2_cis_LT = np.zeros((num_par)),np.zeros((num_par)),np.zeros((num_par)),np.zeros((num_par))
    
    if (t < Record_freq) and (t%100 ==0):
        RecordEcad(t, M1, M2, file_M1_s, file_M2_s, file_record_M1_s, file_record_M2_s)
        M1_MSD_temp, M2_MSD_temp = np.copy(M1),np.copy(M2)
        M1_MSD_temp[:,1:3] += M1_MSD; M2_MSD_temp[:,1:3] += M2_MSD
        RecordEcad(t, M1_MSD_temp, M2_MSD_temp, file_M1_MSD_s, file_M2_MSD_s, file_record_M1, file_record_M2, MSD = True)
        
    M1_prev, M2_prev = np.copy(M1),np.copy(M2)    
    
file_record_M1.close()
file_record_M2.close()
file_M1.close()
file_M2.close()
file_M1_MSD.close()
file_M2_MSD.close()
file_Actin1.close()
file_Actin2.close()

file_record_M1_s.close()
file_record_M2_s.close()
file_M1_s.close()
file_M2_s.close()
file_M1_MSD_s.close()
file_M2_MSD_s.close()
file_M1_cis_LT.close()
file_M2_cis_LT.close()
file_M1_trans_LT.close()
file_M2_trans_LT.close()
print("total time: ", time.time()-start_t)

