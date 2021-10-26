# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 13:12:31 2020

@author: qiliny
"""

import numpy as np
import sys
import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from inputset import num_fila, num_seg_fila, len_seg, side_len_domain, start_index 

def fixmissing(file):
    file_1 = []
    
    for line in file:
        if len(line.split(' ')) == 12:
            line = line+' 0'
            file_1.append(list(map(float, line.rstrip().split(' '))))
        else:
            file_1.append(list(map(float, line.rstrip().split('    '))))
            
    file_1 = np.array(file_1)
    return file_1

# load data
Actin1_data = np.loadtxt('Actin1.txt')[1:num_fila+1,:]
Actin2_data = np.loadtxt('Actin2.txt')[1:num_fila+1,:]

M1_data = np.loadtxt('M1Info.txt')
M2_data = np.loadtxt('M2Info.txt')
# M1_data = open('M1Info.txt','r')
# M2_data = open('M2Info.txt','r')
# M1_data = fixmissing(M1_data)
# M2_data = fixmissing(M2_data)


frame_side_len = side_len_domain
num_par = M1_data[0,0]

num_frames= int(len(M1_data)/(num_par+1))
num_par_index = 0

def ShiftXY(x):
    if x > side_len_domain:
        x -= side_len_domain
    elif x < 0:
        x += side_len_domain 
    return x

def PlotActin(ax, Actin, color):
    # loop over actin filament
    for i in range(0,len(Actin)):
        point1 = Actin[i,1:3]
        point2 = [0,0]
        # loop over segment in the filament
        for j in range(num_seg_fila):
            if j < start_index:
                color = 'b'
            else:
                color = 'k'
        
            point2[0] = point1[0]+np.cos(Actin[i,5])*len_seg
            point2[1] = point1[1]+np.sin(Actin[i,5])*len_seg
            ax.plot([point1[0],point2[0]],[point1[1],point2[1]],c = color)   

            point1 = [ShiftXY(point2[0]),ShiftXY(point2[1])]
            
######################################################
#  movie initialization 
######################################################
dpi=100     # dots per inch, which is the resolution
fig = plt.figure()
ax = fig.add_subplot(111)
ax.axis('off')      

PlotActin(ax,Actin1_data,'b')
PlotActin(ax,Actin2_data,'k')
                        
M1pos = ax.scatter([],[], s=10, c='b')   # M1 monomers
M2pos = ax.scatter([],[], s=10, c='k')   # M2 monomers
transpos = ax.scatter([],[], s=10, c='r')   # Trans dimers

#E_trans_actindata = ax.scatter([],[], s=30, c='b')   # Trans dimers

ax.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
ax.tick_params(axis='y', which='both', bottom=False, top=False, labelleft=False)
fig.set_size_inches(15, 15)   
ax.axis([0, side_len_domain, 0, side_len_domain])

def update_img(n,M1,M2):
    global num_par_index
    print(n)
    num_par = int(M1_data[num_par_index,0])
    
    M1_temp = M1[num_par_index+1:num_par_index+num_par+1,:]
    M2_temp = M2[num_par_index+1:num_par_index+num_par+1,:]
    
    M1_mono_or_cis = (M1_temp[:,4]==1) | (M1_temp[:,7]!=-1) | (M1_temp[:,8]!=-1)
    M2_mono_or_cis = (M2_temp[:,4]==1) | (M2_temp[:,7]!=-1) | (M2_temp[:,8]!=-1)
    
    M1_trans = (M1_temp[:,5]==1) | (M1_temp[:,6]==1)
    
    M1pos.set_offsets(np.c_[M1_temp[np.where(M1_mono_or_cis),1].ravel(),M1_temp[np.where(M1_mono_or_cis),2].ravel()])
    
    M2pos.set_offsets(np.c_[M2_temp[np.where(M2_mono_or_cis),1].ravel(),M2_temp[np.where(M2_mono_or_cis),2].ravel()])
    
    transpos.set_offsets(np.c_[M1_temp[np.where(M1_trans),1].ravel(),M1_temp[np.where(M1_trans),2].ravel()])
    
    num_par_index += num_par+1

    return transpos,M1pos,M2pos


ani = animation.FuncAnimation(fig,update_img,frames=num_frames,fargs=(M1_data,M2_data))
plt.show()
writer = animation.writers['ffmpeg'](fps=10)

ani.save('EcadActinMov.mp4',writer=writer,dpi=dpi)

