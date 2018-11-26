# -*- coding: utf-8 -*-
"""
Created on Wed Nov 21 13:48:21 2018

@author: qiliny

"""

import numpy as np
import math
import random


side_len = 8.0          # Side length of a lattice(nm)
x_len = 125              # The number of lattices alone the x axis
y_len = 125              # The number of lattices alone the y axis
Ecad_conc = 0.16       # Concentration of the E-cad

total_side_len = side_len*x_len         # Side length of the square
radius = total_side_len/2.0             # radius of the circle
num_sites = x_len * y_len                   # The number of lattices on the surface

boundTol = 1e-10
num_par = 100                # Number of Particles

matchedlist = []

def select(data):
    if data != []:
        elem = random.choice(data)
        data.remove(elem)
        return elem
    else:
        return None

#------------------------------------------------------------------------------
# Read the node file
#------------------------------------------------------------------------------
node_file = open('ellipse.1.node.txt', 'r')         # Inputing node file

num_nodes, num_coords, num_attributes, boundary_markers = node_file.readline().split()      # Reading the values of the first line

# Converting the numbers to integers
num_nodes = int(num_nodes)
num_coords = int(num_coords)
num_attributes = int(num_attributes)
boundary_marker = int(boundary_markers)

NodeNums = [[0 for m in range(2)] for n in range(num_nodes)]                # Creating variables to store node number & boundary marker
NodeCoords = [[0 for m in range(num_coords)] for n in range(num_nodes)]     # Creating variable to store x, y and z coordinates
# Reading data from nodefiler
for i in range(num_nodes):
    NodeNums[i][0], NodeCoords[i][0], NodeCoords[i][1], NodeNums[i][1] = node_file.readline().split()
    # Converting from string to int/float
    NodeNums[i][0] = int(NodeNums[i][0])
    NodeNums[i][1] = int(NodeNums[i][1])
    NodeCoords[i][0] = float(NodeCoords[i][0])
    NodeCoords[i][1] = float(NodeCoords[i][1])

NodeNums = np.array(NodeNums)
NodeCoords = np.array(NodeCoords)
Num_BoundaryNodes = [row[1] for row in NodeNums].count(1)# The number of Boundary Nodes

node_file.close()# Closing the file

#------------------------------------------------------------------------------
# Read the element file
#------------------------------------------------------------------------------
elem_file = open('ellipse.1.ele.txt', 'r')

# Reading the values of the first line
num_elements, nodes_per_ele, ele_attributes = elem_file.readline().split()

# Converting the numbers to integers
num_elements = int(num_elements)
nodes_per_ele = int(nodes_per_ele)
ele_attributes = int(ele_attributes)

ElemMap = [[0 for x in range(3)] for y in range(num_elements)]          # Creating variable to store element map
Elemindex = [[0 for m in range(1)] for n in range(num_elements)]
# Reading data from elemfile
for i in range(num_elements):

    Elemindex[i][0], ElemMap[i][0], ElemMap[i][1], ElemMap[i][2] = elem_file.readline().split()

    # Converting from string to int
    Elemindex[i][0] = int(Elemindex[i][0])
    ElemMap[i][0] = int(ElemMap[i][0])
    ElemMap[i][1] = int(ElemMap[i][1])
    ElemMap[i][2] = int(ElemMap[i][2])

elem_file.close()   # Closing the file

#------------------------------------------------------------------------------
# M1 M2 initialization
#------------------------------------------------------------------------------

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

# Area Calculation
def area(x1, y1, x2, y2, x3, y3):
    return abs((x1*(y2-y3)+x2*(y3-y1)+x3*(y1-y2))/2.0)

matchedlist = []    # Elements with corresponding lattice points      

matchedfile = open('MatchLatticeElement.txt','w')

# loop through all elements
for m in range(num_elements):
    matchedfile.write(str(m+1)+'    ')
    templist = []
    # loop through all lattice points
    for n in range(num_sites):
        lat_x = par_x[n]/1000.0     # x of the lattice point
        lat_y = par_y[n]/1000.0     # y of the lattice point
        
        # assign x,y of the three nodes of the element to new variables 
        node1_x, node1_y = NodeCoords[ElemMap[m][0]-1][0], NodeCoords[ElemMap[m][0]-1][1]
        node2_x, node2_y = NodeCoords[ElemMap[m][1]-1][0], NodeCoords[ElemMap[m][1]-1][1]
        node3_x, node3_y = NodeCoords[ElemMap[m][2]-1][0], NodeCoords[ElemMap[m][2]-1][1]
        
        # Area Calculation
        Area = area(node1_x, node1_y, node2_x, node2_y, node3_x, node3_y)   # Area of the element
        Area1 = area(lat_x, lat_y, node2_x, node2_y, node3_x, node3_y)      # Area of the triangle with vertexs P,B,C
        Area2 = area(node1_x, node1_y, lat_x, lat_y, node3_x, node3_y)      # Area of the triangle with vertexs A,P,C
        Area3 = area(node1_x, node1_y, node2_x, node2_y, lat_x, lat_y)      # Area of the triangle with vertexs A,B,P
        
        # if P is in the Element
        if abs(Area-Area1-Area2-Area3)<boundTol:
            templist.append(n) 
    if templist != []:
        matchedfile.write('   '.join(map(str,templist)))
    matchedfile.write('\n')             # file written
    matchedlist.append(templist)        # list stored
matchedfile.close()                     # close file

# Check duplication
# result shows that the program successfully match lattice with nodes except when 
# the position of a lattice point is exactly the same as the position of a node 
seen = set()
duplicated = []
for x in matchedlist:
    for y in x:
        if y not in seen:
            seen.add(y)
        else:
            duplicated.append(y)
            
            
#------------------------------------------------------------------------------
# generate a list, find nearby elements of a node
#------------------------------------------------------------------------------
ElemMap = np.array(ElemMap)
nodenearbyelem = []
nodenearbyelemfile = open('MatchNodeElement.txt','w')
for i in range(num_nodes):
    nodenearbyelemfile.write(str(i+1)+'    ')
    node_elem = np.array(np.where(ElemMap==(i+1)))      # find all elements nearby the node
    nodenearbyelem.append(node_elem[0,:]+1)             # append the elements to the list
    
    nodenearbyelemfile.write('   '.join(map(str,node_elem[0,:]+1)))
    nodenearbyelemfile.write('\n')

nodenearbyelemfile.close()










