# -*- coding: utf-8 -*-
"""
Created on Fri Mar  6 08:27:46 2020

@author: qiliny
"""
import numpy as np

###############################################################################
# INPUTS (Constants & Global Variables)
###############################################################################
Total_time      = 5                                 # s, Time
deltat          = 0.00001                           # s, A single time step
Total_steps     = int(Total_time/deltat)+2          # Total time steps
Record_freq     = 1000                             # Recording frequency in time step

diffusivity     = 0.028*10**6                       # nm^2/s, Diffusivity of the E-cad monomer and trans dimer nm^2/s
steplength      = np.sqrt(4*diffusivity*deltat)   # nm, length of diffusionin each time step nm

num_par         = 1200                               # #, Number of Particles
side_len_domain = 1000                              # nm, side length of the domain
rad_ecad        = 2.5                               # nm, Radius of the E-cad cylinder 

###############################################################################
# Distances & Angle Conditions
###############################################################################

trans_thre      = 1.5                               # nm, mini distance for trans to happen
cis_thre        = 3.0                               # nm, mini distance for ciss to happen
EcadActin_bind_thre  = 3.5                          # nm, mini distance for ciss to happen
collision_thre  = 8

rad_ecad_sq        = rad_ecad**2                    # nm, Radius of the E-cad cylinder 
two_rad_ecad_sq    = (2*rad_ecad)**2                    # nm, Radius of the E-cad cylinder
trans_thre_sq      = trans_thre**2                  # nm, mini distance for trans to happen
cis_thre_sq        = cis_thre**2                    # nm, mini distance for ciss to happen
EcadActin_bind_thre_sq  = EcadActin_bind_thre**2              # nm, mini distance for ciss to happen
collision_thre_sq  = collision_thre**2

PBC_rad_ecad        = side_len_domain-rad_ecad*2                    # nm, Radius of the E-cad cylinder 
PBC_trans_thre      = side_len_domain-trans_thre*2                  # nm, mini distance for trans to happen
PBC_cis_thre        = side_len_domain-cis_thre*2                    # nm, mini distance for ciss to happen
PBC_EcadActin_bind_thre  = side_len_domain-EcadActin_bind_thre*2              # nm, mini distance for ciss to happen
PBC_collision_thre  = side_len_domain-collision_thre*2              # nm, mini distance for ciss to happen

angle_thre_s    = np.pi/6                           # rad, mini angle difference for interactions to happen
angle_thre_l    = angle_thre_s*11                   # rad, mini angle difference for interactions to happen

part_scale_factor = 2                               # Partition scale factor
hop_p           = 0                             # Hopping possibility

###############################################################################
# Actin Constants initialization
###############################################################################
cort_thick = 50                             # (nm) Thickness of cortical actin 
actin_mean  = 80.57; actin_std  = 16.86     # mean and std of the z-pos of alpha-catenin C-terminus
total_actin = 2000                          # Total number of G-actin

actin_angle = 75                            # angle relative to the surface of membrane, unit: degree
actin_bindingeff_ratio = 20/100               # effective ratio of actin binding region

len_fila_2d = cort_thick/np.tan(np.deg2rad(actin_angle))    # Length of F-actin in 2D
len_seg = 2.7*np.cos(np.deg2rad(actin_angle))               # Length of each actin segment
num_seg_fila = int(np.round(len_fila_2d/len_seg))              # Number of segment in each filament
num_fila = int(np.round(total_actin/num_seg_fila))          # Total Number of actin Filaments

start_index = num_seg_fila*actin_bindingeff_ratio

TO_r = 0.2                                 # Actin turnover rate
p_TO_r = TO_r*deltat                        # turnover probability each time step

memb_z   = 35
fencing = 1                                 # 1, actin acts as both fences and binding sites
                                            # 0, actin acts only as binding sites
###############################################################################
# Distribution constants of alpha-catenin
###############################################################################
alpha_cat_mean = 53.83; alpha_cat_std = 8.64    # mean and std of the z-pos of alpha-catenin C-terminus
vinculin_mean  = 87.23; vinculin_std  = 7.95   # mean and std of the z-pos of vinculin C-terminus


linker = 'vinculin'
if linker == 'alpha_cat':
    cad_z_mean, cad_z_std = alpha_cat_mean, alpha_cat_std
elif linker == 'vinculin':
    cad_z_mean, cad_z_std = vinculin_mean, vinculin_std

###############################################################################
# Parameters for position update
###############################################################################
mu = 0; sigma = 1

###############################################################################
# Kinetic Constants, followed by the association and dissociation probability
###############################################################################

k_mx = 38000;           p_k_mx          = 1 - np.exp(-1 * k_mx * deltat)

k_xm = 1840;            p_k_xm          = 1 - np.exp(-1 * k_xm * deltat)

k_xs = 86;              p_k_xs          = 1 - np.exp(-1 * k_xs * deltat)

k_sx = 0.8;             p_k_sx          = 1 - np.exp(-1 * k_sx * deltat)

k_ms = 0.31;            p_k_ms          = 1 - np.exp(-1 * k_ms * deltat)

k_sm = 0.000127;        p_k_sm          = 1 - np.exp(-1 * k_sm * deltat)

k_scis = 100000;          p_k_scis        = 1 - np.exp(-1 * k_scis * deltat)

k_ciss = 1000;             p_k_ciss        = 1 - np.exp(-1 * k_ciss * deltat)

k_cadactin_a = 100000;p_k_cadactin_a  = 1 - np.exp(-1 * k_cadactin_a * deltat)

k_cadactin_d = 0.1;      p_k_cadactin_d  = 1 - np.exp(-1 * k_cadactin_d * deltat)

###############################################################################
# Toggles
###############################################################################
HaveActin_toggle                        = True
ThreeDimensionalActin_toggle            = False     # (False = 2DActin)
ActinTurnover_toggle                    = True
ActinBinding_toggle                     = True  
Rotation_withActinBinding               = True

###############################################################################
# Actin Constants initialization
###############################################################################
# Create Hash tables containing index of M1 and M2 monomers, which were assigned
# to corresponding partition
side_len_par_trans = trans_thre
num_part_side_trans = int(np.floor(side_len_domain/side_len_par_trans))
side_len_par_trans = side_len_domain/num_part_side_trans

side_len_par_cis = cis_thre
num_part_side_cis = int(np.floor(side_len_domain/side_len_par_cis))
side_len_par_cis = side_len_domain/num_part_side_cis

side_len_par_collision = collision_thre
num_part_side_collision = int(np.floor(side_len_domain/side_len_par_collision))
side_len_par_collision = side_len_domain/num_part_side_collision

side_len_par_EcadActin = collision_thre
num_part_side_EcadActin = int(np.floor(side_len_domain/side_len_par_EcadActin))
side_len_par_EcadActin = side_len_domain/num_part_side_EcadActin

def ori_index_hash():
    index_hash = dict()

    for i in range(num_par):
        index_hash[i] = [i]
    return index_hash

par_index_hash = ori_index_hash()



