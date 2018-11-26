import math
import numpy as np
import random
import time
import matplotlib.pyplot as plt
import matplotlib.animation as animation

start = time.time()
t_total_step = 18000    # Time steps(Real time=t_total_step * deltat)
Ecad_conc = 0.16       # Concentration of the E-cad
side_len = 7.9          # Side length of a lattice(nm)
diffusivity = 0.028 * 10 ** 6                # Diffusivity of the E - cad monomer and trans dimer nm ^ 2 / s
deltat = side_len ** 2 / (4 * diffusivity)   # A single time step

x_len = 25              # The number of lattices alone the x axis
y_len = 25              # The number of lattices alone the y axis

total_side_len = side_len*x_len         # Side length of the square
radius = total_side_len/2.0             # radius of the circle
num_sites = x_len * y_len                   # The number of lattices on the surface

g_trans = 6             # Free energy change for trans dimer
g_cis = 4               # Free energy change for cis

gamma_t = math.exp(g_trans)
ass_p_t = gamma_t / (1 + gamma_t)   # Association probability trans
dis_p_t = 1 / (1 + gamma_t)         # Dissociation probability trans

gamma_c = math.exp(g_cis)
ass_p_c = gamma_c / (1 + gamma_c)   # Association probability cis
dis_p_c = 1 / (1 + gamma_c)         # Dissociation probability cis

num_mono = np.zeros((t_total_step / 100, 1), dtype=int)
num_trans = np.zeros((t_total_step / 100, 1), dtype=int)
num_cis = np.zeros((t_total_step / 100, 1), dtype=int)

#----------------------------------------------------------------------------------
# User Defined Functions
#----------------------------------------------------------------------------------
def select(data):
    if data != []:
        elem = random.choice(data)
        data.remove(elem)
        return elem
    else:
        return None

"""
Initialize position and orientations of M1 and M2
col_1: index, col_2: monomer, col_3: trans, col_4: cis
col_5: orientation(N: 5, E: 6, S: 7, W: 8),
col_6: x, col_7: y, col_8: occupied.
"""
# M1[:, 0] = range(num_sites) + np.ones(num_sites)

file_M1 = open("M1_data.txt", "w")
file_M2 = open("M2_data.txt", "w")

M1 = [[0 for x in range(num_sites)] for y in range(8)]  # Initialize M1 matrix
M2 = [[0 for x in range(num_sites)] for y in range(8)]  # Initialize M2 matrix

# assign index and x,y positions
for i in range(x_len*y_len):
    M1[0][i] = i
    M1[5][i] = i % x_len
    M1[6][i] = int(math.floor(i/x_len))

M2[0][:], M2[5][:], M2[6][:] = M1[0][:], M1[5][:], M1[6][:]

# Find all lattice points within the circle
par_x, par_y = [0 for x in range(num_sites)],[0 for x in range(num_sites)]
if x_len%2 != 0:
    par_x = (np.array(M1[5][:])-x_len/2)*side_len
    par_y = (np.array(M1[6][:])-y_len/2)*side_len
else:
    par_x = (np.array(M1[5][:])-(x_len/2-0.5))*side_len
    par_y = (np.array(M1[6][:])-(y_len/2-0.5))*side_len
dist_center = np.sqrt(par_x**2+par_y**2)
incircle = np.array(np.where(np.sqrt(par_x**2+par_y**2)<radius))        # Indexes of lattice points within the circle

# edit properties of occupied sites
num_par = int(Ecad_conc*incircle.size)
incircle_list = map(int,incircle.ravel().tolist())                  # turn it into a list
sites_index = range(x_len*y_len)
occupied_site_M1 = [select(incircle_list) for x in range(num_par)]  # randomly select num_par occupied sites,M1
sites_index = range(x_len*y_len)
occupied_site_M2 = [select(incircle_list) for x in range(num_par)]  # randomly select num_par occupied sites,M2
for i in range(num_par):
    M1[1][occupied_site_M1[i]] = 1
    M1[4][occupied_site_M1[i]] = random.randint(5, 8)
    M1[7][occupied_site_M1[i]] = 1
    
    M2[1][occupied_site_M2[i]] = 1
    M2[4][occupied_site_M2[i]] = random.randint(5, 8)
    M2[7][occupied_site_M2[i]] = 1

M1 = np.transpose(np.array(M1))
M2 = np.transpose(np.array(M2))

# Initialize the arrays
all_bond_trans_old = np.empty((0,3),int)
all_bond_cis_M1_old = np.empty((0,3),int)
all_bond_cis_M2_old = np.empty((0,3),int)

onetwo = np.array([0, 1])
threefour = np.array([2, 3])
col_index = np.array([1, 2, 3, 4, 7])

for t in range(t_total_step):
    entropy_diff = 0;
    attempt = 0;
    
    # Check entropy
    while entropy_diff < random.random() and attempt<5:
        attempt = attempt+1;
        # Filter occupied site to simplify the calculation   
        M1_occupied = M1[np.array(np.where(M1[:,7]==1)).ravel(),:]
        M2_occupied = M2[np.array(np.where(M2[:,7]==1)).ravel(),:]
        
    #######################################################################
    # Formation of dissociation of cis and trans
    # (1: monomer, 2: trans, 3: Cis, 4: Combo)
    #######################################################################
    
    # Formation of trans dimer
        # monomer to trans (1 to 2)
        
        meet_mono_trans = np.array(np.where((M1[:,7]==1)&(M2[:,7]==1)&((M1[:,4]-M2[:,4]==3)|(M1[:,4]-M2[:,5]==-1))))
        new_trans = np.setdiff1d(meet_mono_trans,all_bond_trans_old);        # index of new trans
        new_trans = new_trans[np.random.rand(len(new_trans)) < ass_p_t]
        
    # Formation of cis bond from monomers M1  
    # Formation of cis bond 
        # monomer to cis 
#        trans_trans_M1 = np.empty((0,2))
        new_bond_cis_M1 = np.empty((0,2),int)
        succ_cis_M1 = np.empty((0,1))
        sum_part1part2 = np.array([])
        
        # row 0: site on the left; row 1: site on the right; row 2: site on the bottom; row 3: site on the top 
        check_row_M1 = np.transpose(np.array([M1_occupied[:,0]-1, M1_occupied[:,0]+1, M1_occupied[:,0]-x_len, M1_occupied[:,0]+x_len]))   
#        # check for periodic boundary conditions
        check_row_M1[M1_occupied[:,0]%x_len==0,0] = check_row_M1[M1_occupied[:,0]%x_len==0,0]+x_len
        check_row_M1[M1_occupied[:,0]%x_len==24,1] = check_row_M1[M1_occupied[:,0]%x_len==24,1]-x_len
        check_row_M1[M1_occupied[:,0]<(x_len),2] = check_row_M1[M1_occupied[:,0]<(x_len),2]+x_len*y_len
        check_row_M1[M1_occupied[:,0]>(x_len*y_len-x_len-1),3] = check_row_M1[M1_occupied[:,0]>(x_len*y_len-x_len-1),3]-x_len*y_len

        for i in range(num_par):
            cis_bind_M1 = np.array([])
            if M1_occupied[i,4] in [5, 7]:
                cis_bind_M1 = threefour[M1[check_row_M1[i,[2, 3]],4] == M1_occupied[i,4]]
            elif M1_occupied[i,4] in [6, 8]:
                cis_bind_M1 = onetwo[M1[check_row_M1[i,[0, 1]],4] == M1_occupied[i,4]]
            
#             cis_bind_M1 = find(M1(check_row_M1(i,:),5) == M1_occupied(i,5));    % index of particle next to the particle i
            if cis_bind_M1.size != 0:  
                for j in range(cis_bind_M1.size):
                    part1 = M1_occupied[i,0]
                    part2 = check_row_M1[i,cis_bind_M1[j]]   # index of part 1 and part 2

                    if (all_bond_cis_M1_old.size != 0) and (part1+part2 in all_bond_cis_M1_old[:,2]):
                        continue
                    
                    if (random.random()<ass_p_c) and (part1+part2 not in sum_part1part2):
                        new_bond_cis_M1 = np.append(new_bond_cis_M1, np.array([[part1, part2]]), axis = 0)   # all cis bond on M1 in this time step 
                        sum_part1part2 = np.append(sum_part1part2, part1+part2)
                # append chosen monomer index to the successful binding list        
#                succ_cis_M1 = np.append(succ_cis_M1, [check_row_M1[i,cis_bind_M1]], axis = 0)                # Index of selected Cis
#                succ_cis_M1 = np.append(succ_cis_M1, np.array([[M1_occupied[i,1]]]), axis = 0)
#        
#        [~,ia_M1,~] = unique(sum_part1part2);
#        new_bond_cis_M1 = new_bond_cis_M1(ia_M1,:);     # index of particles in M2 cis bonds
        
    # Formation of cis bond from monomers M2 
        # Formation of cis bond 
        # monomer to cis 
#        trans_trans_M2 = np.empty((0,2))
        new_bond_cis_M2 = np.empty((0,2),int)
        succ_cis_M2 = np.empty((0,1))
        sum_part1part2 = np.empty((0,1))
        
        # row 0: site on the left; row 1: site on the right; row 2: site on the bottom; row 3: site on the top 
        check_row_M2 = np.transpose(np.array([M2_occupied[:,0]-1, M2_occupied[:,0]+1, M2_occupied[:,0]-x_len, M2_occupied[:,0]+x_len]))   
#        # check for periodic boundary conditions
        check_row_M2[M2_occupied[:,0]%x_len==0,0] = check_row_M2[M2_occupied[:,0]%x_len==0,0]+x_len
        check_row_M2[M2_occupied[:,0]%x_len==24,1] = check_row_M2[M2_occupied[:,0]%x_len==24,1]-x_len
        check_row_M2[M2_occupied[:,0]<(x_len),2] = check_row_M2[M2_occupied[:,0]<(x_len),2]+x_len*y_len
        check_row_M2[M2_occupied[:,0]>(x_len*y_len-x_len-1),3] = check_row_M2[M2_occupied[:,0]>(x_len*y_len-x_len-1),3]-x_len*y_len

        for i in range(num_par):
            cis_bind_M2 = np.array([])
            if M2_occupied[i,4] in [5, 7]:
                cis_bind_M2 = threefour[M2[check_row_M2[i,[2, 3]],4] == M2_occupied[i,4]]
            elif M2_occupied[i,4] in [6, 8]:
                cis_bind_M2 = onetwo[M2[check_row_M2[i,[0, 1]],4] == M2_occupied[i,4]]
#             cis_bind_M2 = find(M2(check_row_M2(i,:),5) == M2_occupied(i,5));    % index of particle next to the particle i
            if cis_bind_M2.size != 0:  
                for j in range(cis_bind_M2.size):
                    part1 = M2_occupied[i,0]
                    part2 = check_row_M2[i,cis_bind_M2[j]]   # index of part 1 and part 2

                    if (all_bond_cis_M2_old.size != 0) and (part1+part2 in all_bond_cis_M2_old[:,2]):
                        continue
                    
                    if (random.random()<ass_p_c) and (part1+part2 not in sum_part1part2):
                        new_bond_cis_M2 = np.append(new_bond_cis_M2, np.array([[part1, part2]]), axis = 0)   # all cis bond on M2 in this time step 
                        sum_part1part2 = np.append(sum_part1part2, part1+part2)
                # append chosen monomer index to the successful binding list        
#                succ_cis_M2 = np.append(succ_cis_M2, [check_row_M2[i,cis_bind_M2]], axis = 0)                # Index of selected Cis
#                succ_cis_M2 = np.append(succ_cis_M2, np.array([[M2_occupied[i,1]]]), axis = 0)
#        
#        [~,ia_M2,~] = unique(sum_part1part2);
#        new_bond_cis_M2 = new_bond_cis_M1(ia_M1,:);     # index of particles in M1 cis bonds
       

    # Dissociation of Trans
        trans_index = np.array(np.where(M1[:,2]==1))
        possibility_trans = np.array(np.where([random.random() for x in range(trans_index.size)] < dis_p_t))
        trans_dis = trans_index[possibility_trans] # index of dissociated trans
#        trans_index[np.array([random.random() for x in range(trans_index.size)]) < dis_p_t]
        cis_dis_M1 = np.array([])
        cis_dis_M2 = np.array([])
    # Dissociation of Cis on M1
        if all_bond_cis_M1_old.size != 0:
            cis_dis_M1 = all_bond_cis_M1_old[np.array([random.random() for x in range(all_bond_cis_M1_old.size/3)]) < dis_p_c,0:2]    # index of cis bond that will be dissociated on M1

    # Dissociation of Cis on M2
        if all_bond_cis_M2_old.size != 0:
            cis_dis_M2 = all_bond_cis_M2_old[np.array([random.random() for x in range(all_bond_cis_M2_old.size/3)]) < dis_p_c,0:2]    # index of cis bond that will be dissociated on M1

        entropy_diff = math.exp((new_trans.size-trans_dis.size)*g_trans+(new_bond_cis_M1.size-cis_dis_M1.size+new_bond_cis_M2.size-cis_dis_M2.size)/2*g_cis);
        
    """
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % update positions states 'moving direction' orientations 
    % monomer: position states 'moving direction' orientation
    % trans: position states 'moving direction' 'both orientation'
    % cis: states
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    """
    
    # update association 
    # update trans formation
    M1[new_trans,2] = 1
    M2[new_trans,2] = 1
    
    # update cis formation
    uni_cis_bond_M1, count_M1_cis_for = np.unique(new_bond_cis_M1, return_counts=True)
    M1[uni_cis_bond_M1.astype(int),3] = M1[uni_cis_bond_M1.astype(int),3]+count_M1_cis_for
    uni_cis_bond_M2, count_M2_cis_for = np.unique(new_bond_cis_M2, return_counts=True)
    M2[uni_cis_bond_M2.astype(int),3] = M2[uni_cis_bond_M2.astype(int),3]+count_M2_cis_for
    
    # update trans dissociation
    M1[trans_dis,3] = 0
    M2[trans_dis,3] = 0
    
    # update cis dissociation
    uni_cis_dis_M1, count_M1_cis_dis = np.unique(cis_dis_M1, return_counts=True)
    M1[uni_cis_dis_M1.astype(int),3] = M1[uni_cis_dis_M1.astype(int),3]-count_M1_cis_dis
    uni_cis_dis_M2, count_M2_cis_dis = np.unique(cis_dis_M2, return_counts=True)
    M2[uni_cis_dis_M2.astype(int),3] = M2[uni_cis_dis_M2.astype(int),3]-count_M2_cis_dis

    # code for debug
    if True in (M2[:,3]>2)|(M2[:,3]<0)|(M1[:,3]>2)|(M1[:,3]<0):
        print "Error in row 252\n"
        
        break

    # update monomer
    # Dissociated trans or cis. if a site is occupied but is not the state is not trans or cis, it must be a mono
    M1[(M1[:,2]==0)&(M1[:,3]==0)&(M1[:,7]==1),1] = 1
    M2[(M2[:,2]==0)&(M2[:,3]==0)&(M2[:,7]==1),1] = 1
    # Associated Monomers. if the state of a site is trans or cis, it must not be a mono
    M1[(M1[:,3]!=0)|(M1[:,2]!=0),1] = 0
    M2[(M2[:,3]!=0)|(M2[:,2]!=0),1] = 0

    # update cis bond info, the third column equals to the sum of the first two columns 
    new_bond_cis_M1 = np.c_[new_bond_cis_M1, new_bond_cis_M1[:,0]+new_bond_cis_M1[:,1]]

    new_bond_cis_M2 = np.c_[new_bond_cis_M2, new_bond_cis_M2[:,0]+new_bond_cis_M2[:,1]]

    # Remove dissociated cis bonds from the list
    if cis_dis_M1.size != 0:
        cis_dis_index_M1 = np.array(np.intersect1d((cis_dis_M1[:,0]+cis_dis_M1[:,1]),all_bond_cis_M1_old[:,2], return_indices=True))
        all_bond_cis_M1_old = np.delete(all_bond_cis_M1_old, cis_dis_index_M1[2], axis=0)
 
    if cis_dis_M2.size != 0:
        cis_dis_index_M2 = np.array(np.intersect1d((cis_dis_M2[:,0]+cis_dis_M2[:,1]),all_bond_cis_M2_old[:,2], return_indices=True))
        all_bond_cis_M2_old = np.delete(all_bond_cis_M2_old, cis_dis_index_M2[2], axis=0)


    # Add new cis bonds to the list
    # append newly formed cis bond to the cis bond information list
    all_bond_cis_M1_old = np.append(all_bond_cis_M1_old, new_bond_cis_M1, axis = 0)
    all_bond_cis_M2_old = np.append(all_bond_cis_M2_old, new_bond_cis_M2, axis = 0)
    
    #----------------------------------------------------------------------
    # update M1 position
    #----------------------------------------------------------------------
    M1_occupied_new = np.array(np.where(M1[:,7]==1))
    mono_index_M1 = np.transpose(np.array(np.where(M1[:,1]==1))).ravel()    # index of monomers
    moving_dir_M1 = [random.randint(5, 8) for x in range(mono_index_M1.size)]
    
    M1_pos = np.c_[mono_index_M1-1, mono_index_M1+1, mono_index_M1-x_len, mono_index_M1+x_len]

    # check for periodic boundary conditions
    M1_pos[mono_index_M1%x_len==0,0] = M1_pos[mono_index_M1%x_len==0,0]+x_len
    M1_pos[mono_index_M1%x_len==24,1] = M1_pos[mono_index_M1%x_len==24,1]-x_len
    M1_pos[mono_index_M1<(x_len),2] = M1_pos[mono_index_M1<(x_len),2]+x_len*y_len
    M1_pos[mono_index_M1>(x_len*y_len-x_len-1),3] = M1_pos[mono_index_M1>(x_len*y_len-x_len-1),3]-x_len*y_len

    for i in range(mono_index_M1.size):
        conditions = np.logical_and(np.isin(M1_pos[i,:],M1_occupied_new,invert=True),np.isin(M1_pos[i,:],incircle))
        avai_pos = M1_pos[i, conditions]
        if avai_pos.size != 0:
            desti = avai_pos[random.randint(0, avai_pos.size-1)] 
            
            M1[desti,col_index] = M1[mono_index_M1[i],col_index]
            M1[mono_index_M1[i],col_index] = 0
            M1_occupied_new[M1_occupied_new==mono_index_M1[i]] = desti       # update occupied sites
    
    M1[M1[:,1]==1,4] = [random.randint(5, 8) for x in range(mono_index_M1.size)]          # update orientation

    #----------------------------------------------------------------------
    # update M2 position
    #----------------------------------------------------------------------
    M2_occupied_new = np.array(np.where(M2[:,7]==1))
    mono_index_M2 = np.transpose(np.array(np.where(M2[:,1]==1))).ravel()
    moving_dir_M2 = [random.randint(5, 8) for x in range(mono_index_M2.size)]
    
    M2_pos = np.c_[mono_index_M2-1, mono_index_M2+1, mono_index_M2-x_len, mono_index_M2+x_len]

    # check for periodic boundary conditions
    M2_pos[mono_index_M2%x_len==0,0] = M2_pos[mono_index_M2%x_len==0,0]+x_len
    M2_pos[mono_index_M2%x_len==24,1] = M2_pos[mono_index_M2%x_len==24,1]-x_len
    M2_pos[mono_index_M2<(x_len),2] = M2_pos[mono_index_M2<(x_len),2]+x_len*y_len
    M2_pos[mono_index_M2>(x_len*y_len-x_len-1),3] = M2_pos[mono_index_M2>(x_len*y_len-x_len-1),3]-x_len*y_len


    for i in range(mono_index_M2.size):
        conditions = np.logical_and(np.isin(M2_pos[i,:],M2_occupied_new,invert=True),np.isin(M2_pos[i,:],incircle))
        avai_pos = M2_pos[i, conditions]
        if avai_pos.size != 0:
            desti = avai_pos[random.randint(0, avai_pos.size-1)] 
            
            M2[desti,col_index] = M2[mono_index_M2[i],col_index]
            M2[mono_index_M2[i],col_index] = 0
            M2_occupied_new[M2_occupied_new==mono_index_M2[i]] = desti       # update occupied sites
    
    M2[M2[:,1]==1,4] = [random.randint(5, 8) for x in range(mono_index_M2.size)]          # update orientation

    #----------------------------------------------------------------------
    # update trans position
    #----------------------------------------------------------------------
    all_occupied_new = np.unique(np.append(M1_occupied_new, M2_occupied_new))
    trans_index_movable = np.transpose(np.array(np.where((M1[:,2]==1)&(M2[:,2]==1)&(M1[:,3]==0)&(M2[:,3]==0)))).ravel()
    moving_dir_trans = [random.randint(5, 8) for x in range(trans_index_movable.size)]
    
    trans_pos = np.c_[trans_index_movable-1, trans_index_movable+1, trans_index_movable-x_len, trans_index_movable+x_len]

    # check for periodic boundary conditions
    trans_pos[trans_index_movable%x_len==0,0] = trans_pos[trans_index_movable%x_len==0,0]+x_len
    trans_pos[trans_index_movable%x_len==24,1] = trans_pos[trans_index_movable%x_len==24,1]-x_len
    trans_pos[trans_index_movable<(x_len),2] = trans_pos[trans_index_movable<(x_len),2]+x_len*y_len
    trans_pos[trans_index_movable>(x_len*y_len-x_len-1),3] = trans_pos[trans_index_movable>(x_len*y_len-x_len-1),3]-x_len*y_len

    for i in range(trans_index_movable.size):
        conditions = np.logical_and(np.isin(trans_pos[i,:],all_occupied_new,invert=True),np.isin(trans_pos[i,:],incircle))
        avai_pos = trans_pos[i, conditions]
        if avai_pos.size != 0:
            desti = avai_pos[random.randint(0, avai_pos.size-1)]  
            M1[desti,col_index] = M1[trans_index_movable[i],col_index]
            M1[trans_index_movable[i],col_index] = 0
            M2[desti,col_index] = M2[trans_index_movable[i],col_index]
            M2[trans_index_movable[i],col_index] = 0
            all_occupied_new[all_occupied_new==trans_index_movable[i]] = desti
            
    M1[trans_index_movable,4] = [random.randint(5, 8) for x in range(trans_index_movable.size)]          # update orientation
    
    M2[trans_index_movable,4] = [random.randint(5, 8) for x in range(trans_index_movable.size)]          # update orientation
    M2[M2[:,4]==9,4] = M2[M2[:,4]==9,4]-4
    
    for i in col_index:
        M1[M1[:,7] == 0, i] = 0
        M2[M2[:,7] == 0, i] = 0
        
    if (M1[:,2]!=M2[:,2]).any():
        print 'Error: M1 trans ~= M2 trans\n'

    if t%100==0:
        num_mono[t/100] = np.array(np.where(M1[:,1]==1)).size+np.array(np.where(M2[:,1]==1)).size
        num_trans[t/100] = np.array(np.where(M1[:,2]==1)).size+np.array(np.where(M2[:,2]==1)).size
        num_cis[t/100] = np.array(np.where(M1[:,3]==1)).size+np.array(np.where(M2[:,3]==1)).size
        
#        fig = plt.figure()
#        
#        plt.scatter(M1[:,5]+1,M1[:,6]+1, s=3, c='k')                                    # all Lattices
#        
#        plt.scatter(M1[np.where(M1[:,7]),5]+1,M1[np.where(M1[:,7]),6]+1, s=15, c='r')   # M1 monomers
#         
#        plt.scatter(M2[np.where(M2[:,7]),5]+1,M2[np.where(M2[:,7]),6]+1, s=15, c='g')   # M2 monomers
#         
#        plt.scatter(M1[np.where(M1[:,2]),5]+1,M1[np.where(M1[:,2]),6]+1, s=30, c='b')   # Trans dimers
#    
#        plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
#        plt.tick_params(axis='y', which='both', bottom=False, top=False, labelleft=False)
#        fig.set_size_inches(3, 3)   
#        plt.axis([-0.3, x_len+1, -0.3, y_len+1])
##        imgname = 'scatterplot_%d.avi' % (t/100)
#
#        imgname = 'D:\\2018_02 Year 1\\ResearchProjects\\E-cad dynamics\\latticebasedmodel\\scatterplot\\scatterplot_%d.png' % (t/100)
#
## Save a single image                
#        plt.savefig(imgname, bbox_inches='tight',dpi=200)
#       
#        plt.show()
        
        # backslash "\" can escape itself
        
    

##hold off
##close(v); % Saves the movie.
#num_mono = np.append(num_par*2,num_mono)
#num_trans = np.append(0,num_trans)
#num_cis = np.append(0,num_cis)
##% save(['Variables_' num2str(np) '.mat'])
#
#fig1 = plt.figure()
#ax = fig1.add_axes([0.1, 0.1, 0.8, 0.8])
##set(gca,'Position',[0.17 0.25 0.525 0.7]);
#plt.plot(np.linspace(0,t_total_step*deltat,t_total_step/100+1),num_mono, label='mono')
#plt.plot(np.linspace(0,t_total_step*deltat,t_total_step/100+1),num_trans, label='trans')
#plt.plot(np.linspace(0,t_total_step*deltat,t_total_step/100+1),num_cis, label='cis')
#plt.xlabel('Time (s)')
#plt.ylabel('Number of E-cad')
#fig = plt.gcf()
#fig.set_size_inches(3, 3)
#plt.axis([0, t_total_step*deltat, 0, num_par*2])
#plt.legend()
#
#plt.savefig('Number of E-cad.png', bbox_inches='tight',dpi=800)
#
#plt.show()


#
#figure(2)
#set(gca,'Position',[0.17 0.25 0.525 0.7]);
#plot(0:deltat*100:t_total_step*deltat,num_trans,'linewidth',2)
#xlabel('Time (s)')
#ylabel('Number of Trans')
#b=gca;
#b.FontSize=20;
#fig = gcf;
#fig.PaperPositionMode = 'auto';
#axis([0 t_total_step*deltat 0 np*2])
#print(['NumberofTrans_' num2str(np)],'-dpng','-r0')
#
#figure(3)
#set(gca,'Position',[0.17 0.25 0.525 0.7]);
#plot(0:deltat*100:t_total_step*deltat,num_cis,'linewidth',2)
#xlabel('Time (s)')
#ylabel('Number of Cis')
#c=gca;
#c.FontSize=20;
#fig = gcf;
#fig.PaperPositionMode = 'auto';
#axis([0 t_total_step*deltat 0 np*2])
#print(['NumberofCis_' num2str(np)],'-dpng','-r0')


























end = time.time()
print end-start
            
        








