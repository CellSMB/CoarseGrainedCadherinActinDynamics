#!/usr/bin/env python

import re
import numpy as np
import math
import random
import os
# load up the opencmiss library of functions into an object called iron.
from opencmiss.iron import iron

#_________________________________________________________________________________________________
# INPUTS
#_________________________________________________________________________________________________

#-------------------------------------------------------------------------------------------------
# READING INPUT FILE
#-------------------------------------------------------------------------------------------------

#Open the input file
input_file = open("./input_3.txt","r")


#Create variable that extracts each line from text. Also shifts starting index from 0 to 1.
input_lines = input_file.readlines()
input_lines.insert(0,' ')
#Dimension of the model
number_D = 2

#Geometric Inputs
nodeinput, eleminput                           = re.split(',|\n',input_lines[6])[0:2] #Link to nodes and elements created via triangle

#Simulation runtime inputs
startT, endT, Tstep, ODE_TIME_STEP             	= input_lines[13].split(',')           #Simulation Time Parameters
load_steps                                     	= input_lines[16]                      #Output Frequence of Results
outputfreq                                     	= input_lines[19]                      #Output Frequence of Results

#Reaction Diffusion Inputs
	#E-cad dynamics model

EcadDynamicsModel                               = input_lines[26].split('\n')[0]      #Link to Cellml model
init_M1, Dx_M1, Dy_M1                           = input_lines[29].split(',')		  #M1 inputs
init_M2, Dx_M2, Dy_M2                           = input_lines[32].split(',')		  #M2 inputs
init_MM, Dx_MM, Dy_MM 	   					   	= input_lines[35].split(',')		  #MM inputs
init_X, Dx_X, Dy_X 	 	  					   	= input_lines[38].split(',')		  #X inputs
init_S, Dx_S, Dy_S 	 	  					  	= input_lines[41].split(',')		  #S inputs
init_SS, Dx_SS, Dy_SS 	 	  					= input_lines[44].split(',')		  #SS inputs
init_Cis_star, Dx_Cis_star, Dy_Cis_star 	 	= input_lines[47].split(',')		  #Cis_star inputs
init_Cis, Dx_Cis, Dy_Cis 	 	  			    = input_lines[50].split(',')		  #Cis inputs
k1, k_xp, k_x, k_sp, k_s, k_xs, k_sx 		   	= input_lines[54].split(',')		  #E-cad monomers to dimers kinetic constants
k2, k_cisp, k_cis, k_sa, k_cisp_m, k_cispstar_m, alpha, beta, jon, joff            = input_lines[57].split(',')		  #E-cad dimers to clusters kinetic constants

#Mechanics Inputs
C10, C01, kappa                               	= input_lines[64].split(',')          #Mooney-Rivlin Constants
f_crit										 	= input_lines[67]                     #Critical force for each cadherin dimer to break	
aimed_force							   		   	= input_lines[70]		              #Aimed micropipette_force per um^2 (pN)
aimed_t							   		   	   	= input_lines[73]		              #Time to reach the aimed force (s)

#debug_mode
debug_mode                                      = input_lines[80]

#break_mode
break_mode                                      = input_lines[87]                     #Break Mode				
#Closing File
input_file.close()


#-------------------------------------------------------------------------------------------------
# CONVERTING INPUTS FROM STRING TO APPROPRIATE DATA TYPES
#-------------------------------------------------------------------------------------------------

#Simulation Runtime Inputs
startT                     	= float(startT)
endT                       	= float(endT)
Tstep                      	= float(Tstep)
ODE_TIME_STEP              	= float(ODE_TIME_STEP)
load_steps                 	= int(load_steps)
outputfreq                 	= int(outputfreq)
debug_mode                  = int(debug_mode)
break_mode                  = int(break_mode)

#Reaction Diffusion Inputs
init_M1						= float(init_M1)		#M1 inputs
Dx_M1						= float(Dx_M1)
Dy_M1 	   					= float(Dy_M1)

init_M2						= float(init_M2)		#M2 inputs
Dx_M2						= float(Dx_M2)
Dy_M2 	   					= float(Dy_M2)
		  
init_MM						= float(init_MM)		#MM inputs
Dx_MM						= float(Dx_M2)
Dy_MM 	   					= float(Dy_M2)	
	
init_X						= float(init_X)			#X inputs
Dx_X						= float(Dx_X)
Dy_X 	 	  				= float(Dy_X)		

init_S						= float(init_S)			#S inputs
Dx_S						= float(Dx_S)
Dy_S 	 	  				= float(Dy_S)	

init_SS						= float(init_SS)		#MM inputs
Dx_SS						= float(Dx_SS)
Dy_SS 	   					= float(Dy_SS)	
	
init_Cis_star				= float(init_Cis_star)	#X inputs
Dx_Cis_star					= float(Dx_Cis_star)
Dy_Cis_star 	 	  		= float(Dy_Cis_star)		

init_Cis					= float(init_Cis)	    #S inputs
Dx_Cis						= float(Dx_Cis)
Dy_Cis 	 	  				= float(Dy_Cis)	
	
k1							= float(k1)				#Kinetic constants for monomers and dimers 
k_xp						= float(k_xp)
k_x							= float(k_x)
k_sp						= float(k_sp)
k_s							= float(k_s)
k_xs						= float(k_xs)
k_sx						= float(k_sx)

k2							= float(k2)				#Kinetic constants for dimers and clusters
k_cisp						= float(k_cisp)
k_cis						= float(k_cis)
k_sa						= float(k_sa)
k_cisp_m                    = float(k_cisp_m)
k_cispstar_m                = float(k_cispstar_m)
alpha						= float(alpha)
beta						= float(beta)
jon						    = float(jon)
joff					    = float(joff)
#Mechanics Inputs
C10							= float(C10)
C01							= float(C01)
kappa						= float(kappa)
f_crit						= float(f_crit)
aimed_force					= float(aimed_force)
aimed_t						= float(aimed_t)



#-------------------------------------------------------------------------------------------------
# LOGGING
#-------------------------------------------------------------------------------------------------

if debug_mode == 1:
    #Printing the input types and values to confirm a successful reading
    
    #Geometric Inputs
    print "nodeinput: ",                  type(nodeinput),                 nodeinput
    print "eleminput: ",                  type(eleminput),                 eleminput

    #Simulation Runtime Inputs
    print "startT: ",                     type(startT),                    startT
    print "endT: ",                       type(endT),                      endT
    print "Tstep: ",                      type(Tstep),                     Tstep
    print "ODE_TIME_STEP: ",              type(ODE_TIME_STEP),             ODE_TIME_STEP
    print "load_steps: ",                 type(load_steps),                load_steps
    print "outputfreq: ",                 type(outputfreq),                outputfreq

    #Reaction Diffusion Inputs
    print "EcadDynamicsModel: ",      	  type(EcadDynamicsModel),     	   EcadDynamicsModel
    print "init_M1: ",             		  type(init_M1),            	   init_M1
    print "Dx_M1: ",               		  type(Dx_M1),              	   Dx_M1
    print "Dy_M1: ",               		  type(Dy_M1),              	   Dy_M1
    print "init_M2: ",             		  type(init_M2),            	   init_M2
    print "Dx_M2: ",               		  type(Dx_M2),              	   Dx_M2
    print "Dy_M2: ",               		  type(Dy_M2),              	   Dy_M2
    print "init_MM: ",             		  type(init_MM),            	   init_MM
    print "Dx_MM: ",               		  type(Dx_MM),              	   Dx_MM
    print "Dy_MM: ",               		  type(Dy_MM),              	   Dy_MM
    print "init_X: ",             		  type(init_X),            	   	   init_X
    print "Dx_X: ",               		  type(Dx_X),              	       Dx_X
    print "Dy_X: ",               		  type(Dy_X),              	   	   Dy_X
    print "init_S: ",             		  type(init_S),            	       init_S
    print "Dx_S: ",               		  type(Dx_S),              	       Dx_S
    print "Dy_S: ",               		  type(Dy_S),              	       Dy_S
    print "k1: ",                         type(k1),                        k1
    print "k_xp: ",                       type(k_xp),                      k_xp
    print "k_x: ",                        type(k_x),                       k_x
    print "k_sp: ",                       type(k_sp),                      k_sp
    print "k_s: ",                        type(k_s),                  	   k_s
    print "k_xs: ",                       type(k_xs),                      k_xs
    print "k_sx: ",                       type(k_sx),                      k_sx
    print "k2: ",						  type(k2),                        k2
    print "k_cisp: ",				      type(k_cisp),                    k_cisp
    print "k_cis: ",					  type(k_cis),                     k_cis
    print "k_sa: ",						  type(k_sa),                      k_sa
    print "k_cisp_m: ",				      type(k_cisp_m),                  k_cisp_m
    print "k_cispstar_m: ",				  type(k_cispstar_m),              k_cispstar_m
    print "alpha: ",					  type(alpha),                     alpha
    print "beta: ",						  type(beta),                      beta
    print "jon: ",					      type(jon),                       jon
    print "joff: ",						  type(joff),                      joff
    
    #Mechanics Inputs
    print "C10: ",                        type(C10),                       C10
    print "C01: ",                        type(C01),                       C01
    print "kappa: ",                      type(kappa),                     kappa
    print "f_crit: ",                     type(f_crit),                    f_crit
    print "aimed_force: ",                type(aimed_force),               aimed_force
    print "aimed_t: ",                    type(aimed_t),                   aimed_t

#Defining Storage Coefficient
store_coeff = 1.0

numberOfXi = 2

CoordinateSystemUserNumber = 10
RegionUserNumber = 20
BasisUserNumber = 30
PressureBasisUserNumber = 40
GeneratedMeshUserNumber = 50
MeshUserNumber = 60
DecompositionUserNumber = 70
GeometricFieldUserNumber = 80

Mech_FibreFieldUserNumber = 90
Mech_MaterialsFieldUserNumber = 100
Mech_DependentFieldUserNumber = 110
Mech_EquationsSetUserNumber = 120
Mech_EquationsSetFieldUserNumber = 130

Mech_problemUserNumber				= 140
rd_problemUserNumber				= 150

M1_EquationsSetFieldUserNumber		= 160
M1_EquationsSetUserNumber			= 170
M1_FieldUserNumber					= 180
M1_MaterialsFieldUserNumber			= 190
M1_SourceFieldUserNumber			= 200

M2_EquationsSetFieldUserNumber		= 210
M2_EquationsSetUserNumber			= 220
M2_FieldUserNumber					= 230
M2_MaterialsFieldUserNumber			= 240
M2_SourceFieldUserNumber			= 250

MM_EquationsSetFieldUserNumber		= 260
MM_EquationsSetUserNumber			= 270
MM_FieldUserNumber					= 280
MM_MaterialsFieldUserNumber			= 290
MM_SourceFieldUserNumber			= 300

X_EquationsSetFieldUserNumber		= 310
X_EquationsSetUserNumber			= 320
X_FieldUserNumber					= 330
X_MaterialsFieldUserNumber			= 340
X_SourceFieldUserNumber				= 350

S_EquationsSetFieldUserNumber		= 360
S_EquationsSetUserNumber			= 370
S_FieldUserNumber					= 380
S_MaterialsFieldUserNumber			= 390
S_SourceFieldUserNumber				= 400

SS_EquationsSetFieldUserNumber		= 410
SS_EquationsSetUserNumber			= 420
SS_FieldUserNumber					= 430
SS_MaterialsFieldUserNumber			= 440
SS_SourceFieldUserNumber			= 450

Cis_star_EquationsSetFieldUserNumber= 460
Cis_star_EquationsSetUserNumber		= 470
Cis_star_FieldUserNumber			= 480
Cis_star_MaterialsFieldUserNumber	= 490
Cis_star_SourceFieldUserNumber		= 500

Cis_EquationsSetFieldUserNumber		= 510
Cis_EquationsSetUserNumber			= 520
Cis_FieldUserNumber					= 530
Cis_MaterialsFieldUserNumber	    = 540
Cis_SourceFieldUserNumber			= 550

k1_FieldUserNumber					= 560
k_xp_FieldUserNumber				= 570
k_x_FieldUserNumber					= 580
k_sp_FieldUserNumber				= 590
k_s_FieldUserNumber					= 600
k_xs_FieldUserNumber				= 610
k_sx_FieldUserNumber				= 620

k2_FieldUserNumber					= 630
k_cisp_FieldUserNumber				= 640
k_cis_FieldUserNumber				= 650
k_sa_FieldUserNumber				= 660
k_cisp_m_FieldUserNumber            = 670
k_cispstar_m_FieldUserNumber        = 680
alpha_FieldUserNumber				= 690
beta_FieldUserNumber				= 700
jon_FieldUserNumber 				= 710
joff_FieldUserNumber				= 720

cellmlUserNumber					= 730
cellmlModelsFieldUserNumber         = 740
cellmlStateFieldUserNumber          = 750
cellmlIntermediateFieldUserNumber   = 760
cellmlParametersFieldUserNumber     = 770

# Get the computational nodes information
numberOfComputationalNodes = iron.ComputationalNumberOfNodesGet()
computationalNodeNumber = iron.ComputationalNodeNumberGet()

# -------------------------------------------------------------------------------------------------
# COORDINATE SYSTEM
# -------------------------------------------------------------------------------------------------

# Two Dimensional Coordinate System
coordinateSystem = iron.CoordinateSystem()
coordinateSystem.CreateStart(CoordinateSystemUserNumber)
coordinateSystem.dimension = 2
coordinateSystem.CreateFinish()

# -------------------------------------------------------------------------------------------------
# REGION
# -------------------------------------------------------------------------------------------------

# Start Region
region = iron.Region()
region.CreateStart(RegionUserNumber, iron.WorldRegion)
region.label = "membrane"
region.coordinateSystem = coordinateSystem
region.CreateFinish()

# -------------------------------------------------------------------------------------------------
# BASIS
# -------------------------------------------------------------------------------------------------

# Simplex Basis
basis = iron.Basis()
basis.CreateStart(BasisUserNumber)
basis.type = iron.BasisTypes.SIMPLEX
basis.numberOfXi = numberOfXi
basis.interpolationXi = [iron.BasisInterpolationSpecifications.LINEAR_SIMPLEX] * numberOfXi
# basis.quadratureNumberOfGaussXi = [2]*2
basis.CreateFinish()

# -------------------------------------------------------------------------------------------------
# MESH
# -------------------------------------------------------------------------------------------------

	# ___________________________________________________________________________________
	#                              Initializing Node File

# Inputing node file
node_file = open('ellipse.1.node', 'r')

# Reading the values of the first line
num_nodes, num_coords, num_attributes, boundary_markers = node_file.readline().split()

# Converting the numbers to integers
num_nodes = int(num_nodes)
num_coords = int(num_coords)
num_attributes = int(num_attributes)
boundary_marker = int(boundary_markers)

# Creating variables to store node number & boundary marker
NodeNums = [[0 for m in range(2)] for n in range(num_nodes)]

# Creating variable to store x, y and z coordinates
NodeCoords = [[0 for m in range(num_coords)] for n in range(num_nodes)]

# Reading data from nodefiler
for i in range(num_nodes):
    NodeNums[i][0], NodeCoords[i][0], NodeCoords[i][1], NodeNums[i][1] = node_file.readline().split()
    # Converting from string to int/float
    NodeNums[i][0] = int(NodeNums[i][0])
    NodeNums[i][1] = int(NodeNums[i][1])
    NodeCoords[i][0] = float(NodeCoords[i][0])
    NodeCoords[i][1] = float(NodeCoords[i][1])
    
# The number of Boundary Nodes
Num_BoundaryNodes = [row[1] for row in NodeNums].count(1)

# Closing the file
node_file.close()

# ___________________________________________________________________________________
#                              Initializing Element File

# Inputing element file
elem_file = open('ellipse.1.ele', 'r')

# Reading the values of the first line
num_elements, nodes_per_ele, ele_attributes = elem_file.readline().split()

# Converting the numbers to integers
num_elements = int(num_elements)
nodes_per_ele = int(nodes_per_ele)
ele_attributes = int(ele_attributes)

# Creating variable to store element map
ElemMap = [[0 for x in range(3)] for y in range(num_elements)]
Elemindex = [[0 for m in range(1)] for n in range(num_elements)]
# Reading data from elemfile
for i in range(num_elements):

    Elemindex[i][0], ElemMap[i][0], ElemMap[i][1], ElemMap[i][2] = elem_file.readline().split()

    # Converting from string to int
    Elemindex[i][0] = int(Elemindex[i][0])
    ElemMap[i][0] = int(ElemMap[i][0])
    ElemMap[i][1] = int(ElemMap[i][1])
    ElemMap[i][2] = int(ElemMap[i][2])

# Closing the file
elem_file.close()

# Initialise Nodes
nodes = iron.Nodes()
nodes.CreateStart(region, num_nodes)
nodes.CreateFinish()

# Initialise Mesh
mesh = iron.Mesh()
mesh.CreateStart(MeshUserNumber, region, num_coords)
mesh.NumberOfElementsSet(num_elements)
mesh.NumberOfComponentsSet(1)

# Initialise Elements
meshElements = iron.MeshElements()
meshElements.CreateStart(mesh, 1, basis)
for i in range(num_elements):
    element = Elemindex[i][0]
    meshElements.NodesSet(element, [ElemMap[i][0], ElemMap[i][1], ElemMap[i][2]])
meshElements.CreateFinish()

# Finilise Mesh
mesh.CreateFinish()

# -------------------------------------------------------------------------------------------------
# MESH DECOMPOSITION
# -------------------------------------------------------------------------------------------------

# Parallelization
decomposition = iron.Decomposition()
decomposition.CreateStart(DecompositionUserNumber, mesh)
decomposition.type = iron.DecompositionTypes.CALCULATED
decomposition.numberOfDomains = numberOfComputationalNodes
decomposition.CreateFinish()

# -------------------------------------------------------------------------------------------------
# GEOMETRIC FIELD
# -------------------------------------------------------------------------------------------------

# Geometric Field
geometricField = iron.Field()
geometricField.CreateStart(GeometricFieldUserNumber, region)
geometricField.MeshDecompositionSet(decomposition)
geometricField.TypeSet(iron.FieldTypes.GEOMETRIC)

geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 1, 1)
geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 2, 1)
geometricField.CreateFinish()

# Update Geometric Field from customized mesh
for i in range(num_nodes):
    node = NodeNums[i][0]
    nodeDomain = decomposition.NodeDomainGet(node, 1)
    if nodeDomain == computationalNodeNumber:
        nodex = NodeCoords[i][0]
        nodey = NodeCoords[i][1]
        geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,
                                                iron.FieldParameterSetTypes.VALUES,
                                                1, 1, node, 1, nodex)
        geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,
                                                iron.FieldParameterSetTypes.VALUES,
                                                1, 1, node, 2, nodey)


# Update Geometric Field
geometricField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,
                                       iron.FieldParameterSetTypes.VALUES)
geometricField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,
                                        iron.FieldParameterSetTypes.VALUES)

#_________________________________________________________________________________________________
# MECHANICS 
#_________________________________________________________________________________________________


	# -------------------------------------------------------------------------------------------------
	# FIBRE FIELD
	# -------------------------------------------------------------------------------------------------

# Setting up the fibre field for reaction diffusion equation set.
Mech_FibreField = iron.Field()
Mech_FibreField.CreateStart(Mech_FibreFieldUserNumber, region)
Mech_FibreField.TypeSet(iron.FieldTypes.FIBRE)
Mech_FibreField.MeshDecompositionSet(decomposition)
Mech_FibreField.GeometricFieldSet(geometricField)
Mech_FibreField.NumberOfVariablesSet(1)
Mech_FibreField.NumberOfComponentsSet(iron.FieldVariableTypes.U, 2)
Mech_FibreField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 1, 1)
Mech_FibreField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 2, 1)
Mech_FibreField.VariableLabelSet(iron.FieldVariableTypes.U, "Fibre Field")
Mech_FibreField.CreateFinish()

	# -------------------------------------------------------------------------------------------------
	# DEPENDENT FIELD
	# -------------------------------------------------------------------------------------------------

# Start and label deformation Field Dependent Field
Mech_DependentField = iron.Field()
Mech_DependentField.CreateStart(Mech_DependentFieldUserNumber, region)
Mech_DependentField.TypeSet(iron.FieldTypes.GEOMETRIC_GENERAL)
Mech_DependentField.MeshDecompositionSet(decomposition)
Mech_DependentField.GeometricFieldSet(geometricField)
Mech_DependentField.DependentTypeSet(iron.FieldDependentTypes.DEPENDENT)

Mech_DependentField.NumberOfVariablesSet(2)
Mech_DependentField.NumberOfComponentsSet(iron.FieldVariableTypes.U, 2)
Mech_DependentField.NumberOfComponentsSet(iron.FieldVariableTypes.DELUDELN, 2)
Mech_DependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 1, 1)
Mech_DependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 2, 1)
Mech_DependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.DELUDELN, 1, 1)
Mech_DependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.DELUDELN, 2, 1)
Mech_DependentField.VariableLabelSet(iron.FieldVariableTypes.U, "deformationField")
Mech_DependentField.VariableLabelSet(iron.FieldVariableTypes.DELUDELN, "reaction force Field")

# Finish deformation Field Dependent Field
Mech_DependentField.CreateFinish()

Mech_DependentField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                           iron.FieldParameterSetTypes.VALUES,
                                           1, 0.0)
Mech_DependentField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                           iron.FieldParameterSetTypes.VALUES,
                                           2, 0.0)

	# --------------------------------------------------------------------------------------------
	# MATERIAL FIELD - DEFORMATION
	# --------------------------------------------------------------------------------------------

Mech_MaterialsField = iron.Field()
Mech_MaterialsField.CreateStart(Mech_MaterialsFieldUserNumber, region)
Mech_MaterialsField.TypeSet(iron.FieldTypes.MATERIAL)
Mech_MaterialsField.MeshDecompositionSet(decomposition)
Mech_MaterialsField.GeometricFieldSet(geometricField)
Mech_MaterialsField.NumberOfVariablesSet(1)
Mech_MaterialsField.NumberOfComponentsSet(iron.FieldVariableTypes.U, 3)
Mech_MaterialsField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 1, 1)
Mech_MaterialsField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 2, 1)
Mech_MaterialsField.VariableLabelSet(iron.FieldVariableTypes.U, "deformation Materials Field")
Mech_MaterialsField.CreateFinish()

# Set up the constants in the equation
# Mooney-Rivlin constant C10
Mech_MaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                           iron.FieldParameterSetTypes.VALUES,
                                           1, C10)
# Mooney-Rivlin constant C01
Mech_MaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                           iron.FieldParameterSetTypes.VALUES,
                                           2, C01)
# Mooney-Rivlin constant kappa
Mech_MaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                           iron.FieldParameterSetTypes.VALUES,
                                           3, kappa)
	# -------------------------------------------------------------------------------------------------
	# EQUATIONS SET
	# -------------------------------------------------------------------------------------------------

Mech_EquationsSet = iron.EquationsSet()
Mech_EquationsSetField = iron.Field()
Mech_EquationsSetSpecification = [iron.EquationsSetClasses.ELASTICITY,
                             iron.EquationsSetTypes.FINITE_ELASTICITY,
                             iron.EquationsSetSubtypes.NEARLY_INCOMPRESSIBLE_MOONEY_RIVLIN]


Mech_EquationsSet.CreateStart(Mech_EquationsSetUserNumber, region,
                         Mech_FibreField, Mech_EquationsSetSpecification,
                         Mech_EquationsSetFieldUserNumber, Mech_EquationsSetField)

Mech_EquationsSet.CreateFinish()

Mech_EquationsSet.DependentCreateStart(Mech_DependentFieldUserNumber, Mech_DependentField)
Mech_EquationsSet.DependentCreateFinish()

Mech_EquationsSet.MaterialsCreateStart(Mech_MaterialsFieldUserNumber, Mech_MaterialsField)
Mech_EquationsSet.MaterialsCreateFinish()

	# -------------------------------------------------------------------------------------------------
	# Create equations
	# -------------------------------------------------------------------------------------------------

Mech_Equations = iron.Equations()
Mech_EquationsSet.EquationsCreateStart(Mech_Equations)
Mech_Equations.sparsityType = iron.EquationsSparsityTypes.SPARSE
Mech_Equations.outputType = iron.EquationsOutputTypes.NONE
Mech_EquationsSet.EquationsCreateFinish()

iron.Field.ParametersToFieldParametersComponentCopy(
    geometricField, iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1,
    Mech_DependentField, iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1)

iron.Field.ParametersToFieldParametersComponentCopy(
    geometricField, iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 2,
    Mech_DependentField, iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 2)




#_________________________________________________________________________________________________
# REACTION DIFFUSION
#_________________________________________________________________________________________________

#-------------------------------------------------------------------------------------------------
# M1 EQUATIONS (THE FIRST PROTEIN)
#-------------------------------------------------------------------------------------------------

     #--------------------------------------------------------------------------------------------
     # EQUATION SETS - M1
     #--------------------------------------------------------------------------------------------

# Create standard Laplace equations set
M1_EquationsSetField         = iron.Field()
M1_EquationsSet              = iron.EquationsSet()
M1_EquationsSetSpecification = [iron.EquationsSetClasses.CLASSICAL_FIELD,
                                       iron.EquationsSetTypes.REACTION_DIFFUSION_EQUATION,
                                       iron.EquationsSetSubtypes.CELLML_REAC_SPLIT_REAC_DIFF]
M1_EquationsSet.CreateStart(M1_EquationsSetUserNumber,region,
                                   geometricField,M1_EquationsSetSpecification,
                                   M1_EquationsSetFieldUserNumber,M1_EquationsSetField)
M1_EquationsSet.CreateFinish()



     #--------------------------------------------------------------------------------------------
     # DEPENDENT FIELD - M1
     #--------------------------------------------------------------------------------------------

#Start and label M1 Field Dependent Field
M1_Field = iron.Field()
M1_EquationsSet.DependentCreateStart(M1_FieldUserNumber, M1_Field)
M1_Field.VariableLabelSet(iron.FieldVariableTypes.U,"M1 Field")
M1_Field.VariableLabelSet(iron.FieldVariableTypes.DELUDELN,"M1 Field DELUDELN")

#Finish M1 Field Dependent Field
M1_EquationsSet.DependentCreateFinish()

#Initialise M1 Field Dependent field
M1_Field.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,init_M1)


     #--------------------------------------------------------------------------------------------
     # MATERIAL FIELD - M1
     #--------------------------------------------------------------------------------------------

M1_MaterialsField = iron.Field()
M1_EquationsSet.MaterialsCreateStart(M1_MaterialsFieldUserNumber,M1_MaterialsField)

M1_MaterialsField.VariableLabelSet(iron.FieldVariableTypes.U,"M1 Materials Field")
M1_MaterialsField.ComponentInterpolationSet(iron.FieldVariableTypes.U,1,iron.FieldInterpolationTypes.NODE_BASED)
M1_MaterialsField.ComponentInterpolationSet(iron.FieldVariableTypes.U,2,iron.FieldInterpolationTypes.NODE_BASED)
M1_MaterialsField.ComponentInterpolationSet(iron.FieldVariableTypes.U,3,iron.FieldInterpolationTypes.NODE_BASED)
M1_EquationsSet.MaterialsCreateFinish()


#Diffusion Coefficient in X
M1_MaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                     iron.FieldParameterSetTypes.VALUES,
                                                     1,Dx_M1)

#Diffusion Coefficient in Y
M1_MaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                     iron.FieldParameterSetTypes.VALUES,
                                                     2,Dy_M1)

#Storage Coefficient
M1_MaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                     iron.FieldParameterSetTypes.VALUES,
                                                     3,store_coeff)


     #--------------------------------------------------------------------------------------------
     # SOURCE FIELD - M1
     #--------------------------------------------------------------------------------------------

#Setting up the source field for reaction diffusion equation set.
#For split problem subtype, source field is not used
M1_SourceField = iron.Field()
M1_EquationsSet.SourceCreateStart(M1_SourceFieldUserNumber, M1_SourceField)
M1_SourceField.VariableLabelSet(iron.FieldVariableTypes.U,"M1 Source Field")
M1_EquationsSet.SourceCreateFinish()

#Initialise to 0.0
M1_SourceField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                  iron.FieldParameterSetTypes.VALUES,
                                                  1,0.0)


#Update M1 Source Field
M1_SourceField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,
                                              iron.FieldParameterSetTypes.VALUES)
M1_SourceField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,
                                               iron.FieldParameterSetTypes.VALUES)

#Update M1 Dependent Field
M1_Field.ParameterSetUpdateStart(iron.FieldVariableTypes.U,
                                        iron.FieldParameterSetTypes.VALUES)
M1_Field.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,
                                         iron.FieldParameterSetTypes.VALUES)


#-------------------------------------------------------------------------------------------------
# M2 EQUATIONS (THE SECOND PROTEIN)
#-------------------------------------------------------------------------------------------------

     #--------------------------------------------------------------------------------------------
     # EQUATION SETS - M2
     #--------------------------------------------------------------------------------------------

# Create standard Laplace equations set
M2_EquationsSetField         = iron.Field()
M2_EquationsSet              = iron.EquationsSet()
M2_EquationsSetSpecification = [iron.EquationsSetClasses.CLASSICAL_FIELD,
                                       iron.EquationsSetTypes.REACTION_DIFFUSION_EQUATION,
                                       iron.EquationsSetSubtypes.CELLML_REAC_SPLIT_REAC_DIFF]
M2_EquationsSet.CreateStart(M2_EquationsSetUserNumber,region,
                                   geometricField,M2_EquationsSetSpecification,
                                   M2_EquationsSetFieldUserNumber,M2_EquationsSetField)
M2_EquationsSet.CreateFinish()



     #--------------------------------------------------------------------------------------------
     # DEPENDENT FIELD - M2
     #--------------------------------------------------------------------------------------------

#Start and label M2 Field Dependent Field
M2_Field = iron.Field()
M2_EquationsSet.DependentCreateStart(M2_FieldUserNumber, M2_Field)
M2_Field.VariableLabelSet(iron.FieldVariableTypes.U,"M2 Field")
M2_Field.VariableLabelSet(iron.FieldVariableTypes.DELUDELN,"M2 Field DELUDELN")

#Finish M2 Field Dependent Field
M2_EquationsSet.DependentCreateFinish()

#Initialise M2 Field Dependent field
M2_Field.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,init_M2)


     #--------------------------------------------------------------------------------------------
     # MATERIAL FIELD - M2
     #--------------------------------------------------------------------------------------------

M2_MaterialsField = iron.Field()
M2_EquationsSet.MaterialsCreateStart(M2_MaterialsFieldUserNumber,M2_MaterialsField)

M2_MaterialsField.VariableLabelSet(iron.FieldVariableTypes.U,"M2 Materials Field")
M2_MaterialsField.ComponentInterpolationSet(iron.FieldVariableTypes.U,1,iron.FieldInterpolationTypes.NODE_BASED)
M2_MaterialsField.ComponentInterpolationSet(iron.FieldVariableTypes.U,2,iron.FieldInterpolationTypes.NODE_BASED)
M2_MaterialsField.ComponentInterpolationSet(iron.FieldVariableTypes.U,3,iron.FieldInterpolationTypes.NODE_BASED)
M2_EquationsSet.MaterialsCreateFinish()


#Diffusion Coefficient in X
M2_MaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                     iron.FieldParameterSetTypes.VALUES,
                                                     1,Dx_M2)

#Diffusion Coefficient in Y
M2_MaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                     iron.FieldParameterSetTypes.VALUES,
                                                     2,Dy_M2)

#Storage Coefficient
M2_MaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                     iron.FieldParameterSetTypes.VALUES,
                                                     3,store_coeff)


     #--------------------------------------------------------------------------------------------
     # SOURCE FIELD - M2
     #--------------------------------------------------------------------------------------------

#Setting up the source field for reaction diffusion equation set.
#For split problem subtype, source field is not used
M2_SourceField = iron.Field()
M2_EquationsSet.SourceCreateStart(M2_SourceFieldUserNumber, M2_SourceField)
M2_SourceField.VariableLabelSet(iron.FieldVariableTypes.U,"M2 Source Field")
M2_EquationsSet.SourceCreateFinish()

#Initialise to 0.0
M2_SourceField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                  iron.FieldParameterSetTypes.VALUES,
                                                  1,0.0)


#Update M2 Source Field
M2_SourceField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,
                                              iron.FieldParameterSetTypes.VALUES)
M2_SourceField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,
                                               iron.FieldParameterSetTypes.VALUES)

#Update M2 Dependent Field
M2_Field.ParameterSetUpdateStart(iron.FieldVariableTypes.U,
                                        iron.FieldParameterSetTypes.VALUES)
M2_Field.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,
                                         iron.FieldParameterSetTypes.VALUES)


#-------------------------------------------------------------------------------------------------
# MM EQUATIONS (M1 AND M2 MEET EACH OTHER)
#-------------------------------------------------------------------------------------------------

     #--------------------------------------------------------------------------------------------
     # EQUATION SETS - MM
     #--------------------------------------------------------------------------------------------

# Create standard Laplace equations set
MM_EquationsSetField         = iron.Field()
MM_EquationsSet              = iron.EquationsSet()
MM_EquationsSetSpecification = [iron.EquationsSetClasses.CLASSICAL_FIELD,
                                       iron.EquationsSetTypes.REACTION_DIFFUSION_EQUATION,
                                       iron.EquationsSetSubtypes.CELLML_REAC_SPLIT_REAC_DIFF]
MM_EquationsSet.CreateStart(MM_EquationsSetUserNumber,region,
                                   geometricField,MM_EquationsSetSpecification,
                                   MM_EquationsSetFieldUserNumber,MM_EquationsSetField)
MM_EquationsSet.CreateFinish()



     #--------------------------------------------------------------------------------------------
     # DEPENDENT FIELD - MM
     #--------------------------------------------------------------------------------------------

#Start and label MM Field Dependent Field
MM_Field = iron.Field()
MM_EquationsSet.DependentCreateStart(MM_FieldUserNumber, MM_Field)
MM_Field.VariableLabelSet(iron.FieldVariableTypes.U,"MM Field")
MM_Field.VariableLabelSet(iron.FieldVariableTypes.DELUDELN,"MM Field DELUDELN")

#Finish MM Field Dependent Field
MM_EquationsSet.DependentCreateFinish()

#Initialise MM Field Dependent field
MM_Field.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,init_MM)


     #--------------------------------------------------------------------------------------------
     # MATERIAL FIELD - MM
     #--------------------------------------------------------------------------------------------

MM_MaterialsField = iron.Field()
MM_EquationsSet.MaterialsCreateStart(MM_MaterialsFieldUserNumber,MM_MaterialsField)

MM_MaterialsField.VariableLabelSet(iron.FieldVariableTypes.U,"MM Materials Field")
MM_MaterialsField.ComponentInterpolationSet(iron.FieldVariableTypes.U,1,iron.FieldInterpolationTypes.NODE_BASED)
MM_MaterialsField.ComponentInterpolationSet(iron.FieldVariableTypes.U,2,iron.FieldInterpolationTypes.NODE_BASED)
MM_MaterialsField.ComponentInterpolationSet(iron.FieldVariableTypes.U,3,iron.FieldInterpolationTypes.NODE_BASED)
MM_EquationsSet.MaterialsCreateFinish()


#Diffusion Coefficient in X
MM_MaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                     iron.FieldParameterSetTypes.VALUES,
                                                     1,Dx_MM)

#Diffusion Coefficient in Y
MM_MaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                     iron.FieldParameterSetTypes.VALUES,
                                                     2,Dy_MM)

#Storage Coefficient
MM_MaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                     iron.FieldParameterSetTypes.VALUES,
                                                     3,store_coeff)


     #--------------------------------------------------------------------------------------------
     # SOURCE FIELD - MM
     #--------------------------------------------------------------------------------------------

#Setting up the source field for reaction diffusion equation set.
#For split problem subtype, source field is not used
MM_SourceField = iron.Field()
MM_EquationsSet.SourceCreateStart(MM_SourceFieldUserNumber, MM_SourceField)
MM_SourceField.VariableLabelSet(iron.FieldVariableTypes.U,"MM Source Field")
MM_EquationsSet.SourceCreateFinish()

#Initialise to 0.0
MM_SourceField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                  iron.FieldParameterSetTypes.VALUES,
                                                  1,0.0)


#Update MM Source Field
MM_SourceField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,
                                              iron.FieldParameterSetTypes.VALUES)
MM_SourceField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,
                                               iron.FieldParameterSetTypes.VALUES)

#Update MM Dependent Field
MM_Field.ParameterSetUpdateStart(iron.FieldVariableTypes.U,
                                        iron.FieldParameterSetTypes.VALUES)
MM_Field.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,
                                         iron.FieldParameterSetTypes.VALUES)


#-------------------------------------------------------------------------------------------------
# X EQUATIONS (X SHAPE DIMER OF M1 & M2)
#-------------------------------------------------------------------------------------------------

     #--------------------------------------------------------------------------------------------
     # EQUATION SETS - X
     #--------------------------------------------------------------------------------------------

# Create standard Laplace equations set
X_EquationsSetField         = iron.Field()
X_EquationsSet              = iron.EquationsSet()
X_EquationsSetSpecification = [iron.EquationsSetClasses.CLASSICAL_FIELD,
                                       iron.EquationsSetTypes.REACTION_DIFFUSION_EQUATION,
                                       iron.EquationsSetSubtypes.CELLML_REAC_SPLIT_REAC_DIFF]
X_EquationsSet.CreateStart(X_EquationsSetUserNumber,region,
                                   geometricField,X_EquationsSetSpecification,
                                   X_EquationsSetFieldUserNumber,X_EquationsSetField)
X_EquationsSet.CreateFinish()



     #--------------------------------------------------------------------------------------------
     # DEPENDENT FIELD - X
     #--------------------------------------------------------------------------------------------

#Start and label X Field Dependent Field
X_Field = iron.Field()
X_EquationsSet.DependentCreateStart(X_FieldUserNumber, X_Field)
X_Field.VariableLabelSet(iron.FieldVariableTypes.U,"X Field")
X_Field.VariableLabelSet(iron.FieldVariableTypes.DELUDELN,"X Field DELUDELN")

#Finish X Field Dependent Field
X_EquationsSet.DependentCreateFinish()

#Initialise X Field Dependent field
X_Field.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,init_X)


     #--------------------------------------------------------------------------------------------
     # MATERIAL FIELD - X
     #--------------------------------------------------------------------------------------------

X_MaterialsField = iron.Field()
X_EquationsSet.MaterialsCreateStart(X_MaterialsFieldUserNumber,X_MaterialsField)
X_MaterialsField.VariableLabelSet(iron.FieldVariableTypes.U,"X Materials Field")
X_MaterialsField.ComponentInterpolationSet(iron.FieldVariableTypes.U,1,iron.FieldInterpolationTypes.NODE_BASED)
X_MaterialsField.ComponentInterpolationSet(iron.FieldVariableTypes.U,2,iron.FieldInterpolationTypes.NODE_BASED)
X_MaterialsField.ComponentInterpolationSet(iron.FieldVariableTypes.U,3,iron.FieldInterpolationTypes.NODE_BASED)
X_EquationsSet.MaterialsCreateFinish()

#Diffusion Coefficient in X
X_MaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                     iron.FieldParameterSetTypes.VALUES,
                                                     1,Dx_X)

#Diffusion Coefficient in Y
X_MaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                     iron.FieldParameterSetTypes.VALUES,
                                                     2,Dy_X)

#Storage Coefficient
X_MaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                     iron.FieldParameterSetTypes.VALUES,
                                                     3,store_coeff)
#for i in range(num_nodes):
#    node = NodeNums[i][0]
#    #Diffusion Coefficient in X
#    X_MaterialsField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,node,1,Dx_X)

#    #Diffusion Coefficient in Y
#    X_MaterialsField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,node,2,Dy_X)

#    #Storage Coefficient
#    X_MaterialsField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,node,3,store_coeff)




     #--------------------------------------------------------------------------------------------
     # SOURCE FIELD - X
     #--------------------------------------------------------------------------------------------

#Setting up the source field for reaction diffusion equation set.
#For split problem subtype, source field is not used
X_SourceField = iron.Field()
X_EquationsSet.SourceCreateStart(X_SourceFieldUserNumber, X_SourceField)
X_SourceField.VariableLabelSet(iron.FieldVariableTypes.U,"X Source Field")
X_EquationsSet.SourceCreateFinish()

#Initialise to 0.0
X_SourceField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                  iron.FieldParameterSetTypes.VALUES,
                                                  1,0.0)


#Update X Source Field
X_SourceField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,
                                              iron.FieldParameterSetTypes.VALUES)
X_SourceField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,
                                               iron.FieldParameterSetTypes.VALUES)

#Update X Dependent Field
X_Field.ParameterSetUpdateStart(iron.FieldVariableTypes.U,
                                        iron.FieldParameterSetTypes.VALUES)
X_Field.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,
                                         iron.FieldParameterSetTypes.VALUES)


#-------------------------------------------------------------------------------------------------
# S EQUATIONS (STRAND-SWAPPED DIMER)
#-------------------------------------------------------------------------------------------------

     #--------------------------------------------------------------------------------------------
     # EQUATION SETS - S
     #--------------------------------------------------------------------------------------------

# Create standard Laplace equations set
S_EquationsSetField         = iron.Field()
S_EquationsSet              = iron.EquationsSet()
S_EquationsSetSpecification = [iron.EquationsSetClasses.CLASSICAL_FIELD,
                                       iron.EquationsSetTypes.REACTION_DIFFUSION_EQUATION,
                                       iron.EquationsSetSubtypes.CELLML_REAC_SPLIT_REAC_DIFF]
S_EquationsSet.CreateStart(S_EquationsSetUserNumber,region,
                                   geometricField,S_EquationsSetSpecification,
                                   S_EquationsSetFieldUserNumber,S_EquationsSetField)
S_EquationsSet.CreateFinish()



     #--------------------------------------------------------------------------------------------
     # DEPENDENT FIELD - S
     #--------------------------------------------------------------------------------------------

#Start and label S Field Dependent Field
S_Field = iron.Field()
S_EquationsSet.DependentCreateStart(S_FieldUserNumber, S_Field)
S_Field.VariableLabelSet(iron.FieldVariableTypes.U,"S Field")
S_Field.VariableLabelSet(iron.FieldVariableTypes.DELUDELN,"S Field DELUDELN")

#Finish S Field Dependent Field
S_EquationsSet.DependentCreateFinish()

#Initialise S Field Dependent field
S_Field.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,init_S)


     #--------------------------------------------------------------------------------------------
     # MATERIAL FIELD - S
     #--------------------------------------------------------------------------------------------

S_MaterialsField = iron.Field()
S_EquationsSet.MaterialsCreateStart(S_MaterialsFieldUserNumber,S_MaterialsField)
S_MaterialsField.VariableLabelSet(iron.FieldVariableTypes.U,"S Materials Field")
S_MaterialsField.ComponentInterpolationSet(iron.FieldVariableTypes.U,1,iron.FieldInterpolationTypes.NODE_BASED)
S_MaterialsField.ComponentInterpolationSet(iron.FieldVariableTypes.U,2,iron.FieldInterpolationTypes.NODE_BASED)
S_MaterialsField.ComponentInterpolationSet(iron.FieldVariableTypes.U,3,iron.FieldInterpolationTypes.NODE_BASED)
S_EquationsSet.MaterialsCreateFinish()


#Diffusion Coefficient in X
S_MaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                     iron.FieldParameterSetTypes.VALUES,
                                                     1,Dx_S)

#Diffusion Coefficient in Y
S_MaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                     iron.FieldParameterSetTypes.VALUES,
                                                     2,Dy_S)

#Storage Coefficient
S_MaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                     iron.FieldParameterSetTypes.VALUES,
                                                     3,store_coeff)
#for i in range(num_nodes):
#    node = NodeNums[i][0]
#    #Diffusion Coefficient in X
#    S_MaterialsField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,node,1,Dx_S)

#    #Diffusion Coefficient in Y
#    S_MaterialsField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,node,2,Dy_S)

#    #Storage Coefficient
#    S_MaterialsField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,node,3,store_coeff)

     #--------------------------------------------------------------------------------------------
     # SOURCE FIELD - S
     #--------------------------------------------------------------------------------------------

#Setting up the source field for reaction diffusion equation set.
#For split problem subtype, source field is not used
S_SourceField = iron.Field()
S_EquationsSet.SourceCreateStart(S_SourceFieldUserNumber, S_SourceField)
S_SourceField.VariableLabelSet(iron.FieldVariableTypes.U,"S Source Field")
S_EquationsSet.SourceCreateFinish()

#Initialise to 0.0
S_SourceField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                  iron.FieldParameterSetTypes.VALUES,
                                                  1,0.0)


#Update S Source Field
S_SourceField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,
                                              iron.FieldParameterSetTypes.VALUES)
S_SourceField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,
                                               iron.FieldParameterSetTypes.VALUES)

#Update S Dependent Field
S_Field.ParameterSetUpdateStart(iron.FieldVariableTypes.U,
                                        iron.FieldParameterSetTypes.VALUES)
S_Field.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,
                                         iron.FieldParameterSetTypes.VALUES)

#-------------------------------------------------------------------------------------------------
# SS EQUATIONS (THE FIRST PROTEIN)
#-------------------------------------------------------------------------------------------------

     #--------------------------------------------------------------------------------------------
     # EQUATION SETS - SS
     #--------------------------------------------------------------------------------------------

# Create standard Laplace equations set
SS_EquationsSetField         = iron.Field()
SS_EquationsSet              = iron.EquationsSet()
SS_EquationsSetSpecification = [iron.EquationsSetClasses.CLASSICAL_FIELD,
                                       iron.EquationsSetTypes.REACTION_DIFFUSION_EQUATION,
                                       iron.EquationsSetSubtypes.CELLML_REAC_SPLIT_REAC_DIFF]
SS_EquationsSet.CreateStart(SS_EquationsSetUserNumber,region,
                                   geometricField,SS_EquationsSetSpecification,
                                   SS_EquationsSetFieldUserNumber,SS_EquationsSetField)
SS_EquationsSet.CreateFinish()



     #--------------------------------------------------------------------------------------------
     # DEPENDENT FIELD - SS
     #--------------------------------------------------------------------------------------------

#Start and label SS Field Dependent Field
SS_Field = iron.Field()
SS_EquationsSet.DependentCreateStart(SS_FieldUserNumber, SS_Field)
SS_Field.VariableLabelSet(iron.FieldVariableTypes.U,"SS Field")
SS_Field.VariableLabelSet(iron.FieldVariableTypes.DELUDELN,"SS Field DELUDELN")

#Finish SS Field Dependent Field
SS_EquationsSet.DependentCreateFinish()

#Initialise SS Field Dependent field
SS_Field.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,init_SS)


     #--------------------------------------------------------------------------------------------
     # MATERIAL FIELD - SS
     #--------------------------------------------------------------------------------------------

SS_MaterialsField = iron.Field()
SS_EquationsSet.MaterialsCreateStart(SS_MaterialsFieldUserNumber,SS_MaterialsField)

SS_MaterialsField.VariableLabelSet(iron.FieldVariableTypes.U,"SS Materials Field")
SS_MaterialsField.ComponentInterpolationSet(iron.FieldVariableTypes.U,1,iron.FieldInterpolationTypes.NODE_BASED)
SS_MaterialsField.ComponentInterpolationSet(iron.FieldVariableTypes.U,2,iron.FieldInterpolationTypes.NODE_BASED)
SS_MaterialsField.ComponentInterpolationSet(iron.FieldVariableTypes.U,3,iron.FieldInterpolationTypes.NODE_BASED)
SS_EquationsSet.MaterialsCreateFinish()


#Diffusion Coefficient in X
SS_MaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                     iron.FieldParameterSetTypes.VALUES,
                                                     1,Dx_SS)

#Diffusion Coefficient in Y
SS_MaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                     iron.FieldParameterSetTypes.VALUES,
                                                     2,Dy_SS)

#Storage Coefficient
SS_MaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                     iron.FieldParameterSetTypes.VALUES,
                                                     3,store_coeff)


     #--------------------------------------------------------------------------------------------
     # SOURCE FIELD - SS
     #--------------------------------------------------------------------------------------------

#Setting up the source field for reaction diffusion equation set.
#For split problem subtype, source field is not used
SS_SourceField = iron.Field()
SS_EquationsSet.SourceCreateStart(SS_SourceFieldUserNumber, SS_SourceField)
SS_SourceField.VariableLabelSet(iron.FieldVariableTypes.U,"SS Source Field")
SS_EquationsSet.SourceCreateFinish()

#Initialise to 0.0
SS_SourceField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                  iron.FieldParameterSetTypes.VALUES,
                                                  1,0.0)


#Update SS Source Field
SS_SourceField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,
                                              iron.FieldParameterSetTypes.VALUES)
SS_SourceField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,
                                               iron.FieldParameterSetTypes.VALUES)

#Update SS Dependent Field
SS_Field.ParameterSetUpdateStart(iron.FieldVariableTypes.U,
                                        iron.FieldParameterSetTypes.VALUES)
SS_Field.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,
                                         iron.FieldParameterSetTypes.VALUES)

#-------------------------------------------------------------------------------------------------
# Cis_star EQUATIONS (THE FIRST PROTEIN)
#-------------------------------------------------------------------------------------------------

     #--------------------------------------------------------------------------------------------
     # EQUATION SETS - Cis_star
     #--------------------------------------------------------------------------------------------

# Create standard Laplace equations set
Cis_star_EquationsSetField         = iron.Field()
Cis_star_EquationsSet              = iron.EquationsSet()
Cis_star_EquationsSetSpecification = [iron.EquationsSetClasses.CLASSICAL_FIELD,
                                       iron.EquationsSetTypes.REACTION_DIFFUSION_EQUATION,
                                       iron.EquationsSetSubtypes.CELLML_REAC_SPLIT_REAC_DIFF]
Cis_star_EquationsSet.CreateStart(Cis_star_EquationsSetUserNumber,region,
                                   geometricField,Cis_star_EquationsSetSpecification,
                                   Cis_star_EquationsSetFieldUserNumber,Cis_star_EquationsSetField)
Cis_star_EquationsSet.CreateFinish()



     #--------------------------------------------------------------------------------------------
     # DEPENDENT FIELD - Cis_star
     #--------------------------------------------------------------------------------------------

#Start and label Cis_star Field Dependent Field
Cis_star_Field = iron.Field()
Cis_star_EquationsSet.DependentCreateStart(Cis_star_FieldUserNumber, Cis_star_Field)
Cis_star_Field.VariableLabelSet(iron.FieldVariableTypes.U,"Cis_star Field")
Cis_star_Field.VariableLabelSet(iron.FieldVariableTypes.DELUDELN,"Cis_star Field DELUDELN")

#Finish Cis_star Field Dependent Field
Cis_star_EquationsSet.DependentCreateFinish()

#Initialise Cis_star Field Dependent field
Cis_star_Field.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,init_Cis_star)


     #--------------------------------------------------------------------------------------------
     # MATERIAL FIELD - Cis_star
     #--------------------------------------------------------------------------------------------

Cis_star_MaterialsField = iron.Field()
Cis_star_EquationsSet.MaterialsCreateStart(Cis_star_MaterialsFieldUserNumber,Cis_star_MaterialsField)

Cis_star_MaterialsField.VariableLabelSet(iron.FieldVariableTypes.U,"Cis_star Materials Field")
Cis_star_MaterialsField.ComponentInterpolationSet(iron.FieldVariableTypes.U,1,iron.FieldInterpolationTypes.NODE_BASED)
Cis_star_MaterialsField.ComponentInterpolationSet(iron.FieldVariableTypes.U,2,iron.FieldInterpolationTypes.NODE_BASED)
Cis_star_MaterialsField.ComponentInterpolationSet(iron.FieldVariableTypes.U,3,iron.FieldInterpolationTypes.NODE_BASED)
Cis_star_EquationsSet.MaterialsCreateFinish()


#Diffusion Coefficient in X
Cis_star_MaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                     iron.FieldParameterSetTypes.VALUES,
                                                     1,Dx_Cis_star)

#Diffusion Coefficient in Y
Cis_star_MaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                     iron.FieldParameterSetTypes.VALUES,
                                                     2,Dy_Cis_star)

#Storage Coefficient
Cis_star_MaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                     iron.FieldParameterSetTypes.VALUES,
                                                     3,store_coeff)


     #--------------------------------------------------------------------------------------------
     # SOURCE FIELD - Cis_star
     #--------------------------------------------------------------------------------------------

#Setting up the source field for reaction diffusion equation set.
#For split problem subtype, source field is not used
Cis_star_SourceField = iron.Field()
Cis_star_EquationsSet.SourceCreateStart(Cis_star_SourceFieldUserNumber, Cis_star_SourceField)
Cis_star_SourceField.VariableLabelSet(iron.FieldVariableTypes.U,"Cis_star Source Field")
Cis_star_EquationsSet.SourceCreateFinish()

#Initialise to 0.0
Cis_star_SourceField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                  iron.FieldParameterSetTypes.VALUES,
                                                  1,0.0)


#Update Cis_star Source Field
Cis_star_SourceField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,
                                              iron.FieldParameterSetTypes.VALUES)
Cis_star_SourceField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,
                                               iron.FieldParameterSetTypes.VALUES)

#Update Cis_star Dependent Field
Cis_star_Field.ParameterSetUpdateStart(iron.FieldVariableTypes.U,
                                        iron.FieldParameterSetTypes.VALUES)
Cis_star_Field.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,
                                         iron.FieldParameterSetTypes.VALUES)

#-------------------------------------------------------------------------------------------------
# Cis EQUATIONS (THE FIRST PROTEIN)
#-------------------------------------------------------------------------------------------------

     #--------------------------------------------------------------------------------------------
     # EQUATION SETS - Cis
     #--------------------------------------------------------------------------------------------

# Create standard Laplace equations set
Cis_EquationsSetField         = iron.Field()
Cis_EquationsSet              = iron.EquationsSet()
Cis_EquationsSetSpecification = [iron.EquationsSetClasses.CLASSICAL_FIELD,
                                       iron.EquationsSetTypes.REACTION_DIFFUSION_EQUATION,
                                       iron.EquationsSetSubtypes.CELLML_REAC_SPLIT_REAC_DIFF]
Cis_EquationsSet.CreateStart(Cis_EquationsSetUserNumber,region,
                                   geometricField,Cis_EquationsSetSpecification,
                                   Cis_EquationsSetFieldUserNumber,Cis_EquationsSetField)
Cis_EquationsSet.CreateFinish()



     #--------------------------------------------------------------------------------------------
     # DEPENDENT FIELD - Cis
     #--------------------------------------------------------------------------------------------

#Start and label Cis Field Dependent Field
Cis_Field = iron.Field()
Cis_EquationsSet.DependentCreateStart(Cis_FieldUserNumber, Cis_Field)
Cis_Field.VariableLabelSet(iron.FieldVariableTypes.U,"Cis Field")
Cis_Field.VariableLabelSet(iron.FieldVariableTypes.DELUDELN,"Cis Field DELUDELN")

#Finish Cis Field Dependent Field
Cis_EquationsSet.DependentCreateFinish()

#Initialise Cis Field Dependent field
Cis_Field.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,init_Cis)


     #--------------------------------------------------------------------------------------------
     # MATERIAL FIELD - Cis
     #--------------------------------------------------------------------------------------------

Cis_MaterialsField = iron.Field()
Cis_EquationsSet.MaterialsCreateStart(Cis_MaterialsFieldUserNumber,Cis_MaterialsField)

Cis_MaterialsField.VariableLabelSet(iron.FieldVariableTypes.U,"Cis Materials Field")
Cis_MaterialsField.ComponentInterpolationSet(iron.FieldVariableTypes.U,1,iron.FieldInterpolationTypes.NODE_BASED)
Cis_MaterialsField.ComponentInterpolationSet(iron.FieldVariableTypes.U,2,iron.FieldInterpolationTypes.NODE_BASED)
Cis_MaterialsField.ComponentInterpolationSet(iron.FieldVariableTypes.U,3,iron.FieldInterpolationTypes.NODE_BASED)
Cis_EquationsSet.MaterialsCreateFinish()


#Diffusion Coefficient in X
Cis_MaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                     iron.FieldParameterSetTypes.VALUES,
                                                     1,Dx_Cis)

#Diffusion Coefficient in Y
Cis_MaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                     iron.FieldParameterSetTypes.VALUES,
                                                     2,Dy_Cis)

#Storage Coefficient
Cis_MaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                     iron.FieldParameterSetTypes.VALUES,
                                                     3,store_coeff)


     #--------------------------------------------------------------------------------------------
     # SOURCE FIELD - Cis
     #--------------------------------------------------------------------------------------------

#Setting up the source field for reaction diffusion equation set.
#For split problem subtype, source field is not used
Cis_SourceField = iron.Field()
Cis_EquationsSet.SourceCreateStart(Cis_SourceFieldUserNumber, Cis_SourceField)
Cis_SourceField.VariableLabelSet(iron.FieldVariableTypes.U,"Cis Source Field")
Cis_EquationsSet.SourceCreateFinish()

#Initialise to 0.0
Cis_SourceField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                  iron.FieldParameterSetTypes.VALUES,
                                                  1,0.0)


#Update Cis Source Field
Cis_SourceField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,
                                              iron.FieldParameterSetTypes.VALUES)
Cis_SourceField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,
                                               iron.FieldParameterSetTypes.VALUES)

#Update Cis Dependent Field
Cis_Field.ParameterSetUpdateStart(iron.FieldVariableTypes.U,
                                        iron.FieldParameterSetTypes.VALUES)
Cis_Field.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,
                                         iron.FieldParameterSetTypes.VALUES)


#-------------------------------------------------------------------------------------------------
# CONSTANT FIELDS
#-------------------------------------------------------------------------------------------------

     #--------------------------------------------------------------------------------------------
     # k1 FIELD
     #--------------------------------------------------------------------------------------------

k1_Field = iron.Field()
k1_Field.CreateStart(k1_FieldUserNumber,region)

k1_Field.TypeSet(iron.FieldTypes.GENERAL)
k1_Field.MeshDecompositionSet(decomposition)
k1_Field.GeometricFieldSet(geometricField)
k1_Field.NumberOfVariablesSet(1)
k1_Field.VariableTypesSet([iron.FieldVariableTypes.U])
k1_Field.DataTypeSet(iron.FieldVariableTypes.U,iron.FieldDataTypes.DP)
k1_Field.DimensionSet(iron.FieldVariableTypes.U,iron.FieldDimensionTypes.SCALAR)
k1_Field.NumberOfComponentsSet(iron.FieldVariableTypes.U,1)
k1_Field.VariableLabelSet(iron.FieldVariableTypes.U,"k1 Field")

geometricMeshComponent = geometricField.ComponentMeshComponentGet(iron.FieldVariableTypes.U,1)
k1_Field.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,geometricMeshComponent)
k1_Field.ComponentInterpolationSet(iron.FieldVariableTypes.U,1,
                                      iron.FieldInterpolationTypes.NODE_BASED)

k1_Field.CreateFinish()

k1_Field.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                        iron.FieldParameterSetTypes.VALUES,
                                        1,k1)

	 #--------------------------------------------------------------------------------------------
     # k_xp FIELD
     #--------------------------------------------------------------------------------------------

k_xp_Field = iron.Field()
k_xp_Field.CreateStart(k_xp_FieldUserNumber,region)

k_xp_Field.TypeSet(iron.FieldTypes.GENERAL)
k_xp_Field.MeshDecompositionSet(decomposition)
k_xp_Field.GeometricFieldSet(geometricField)
k_xp_Field.NumberOfVariablesSet(1)
k_xp_Field.VariableTypesSet([iron.FieldVariableTypes.U])
k_xp_Field.DataTypeSet(iron.FieldVariableTypes.U,iron.FieldDataTypes.DP)
k_xp_Field.DimensionSet(iron.FieldVariableTypes.U,iron.FieldDimensionTypes.SCALAR)
k_xp_Field.NumberOfComponentsSet(iron.FieldVariableTypes.U,1)
k_xp_Field.VariableLabelSet(iron.FieldVariableTypes.U,"k_xp Field")

geometricMeshComponent = geometricField.ComponentMeshComponentGet(iron.FieldVariableTypes.U,1)
k_xp_Field.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,geometricMeshComponent)
k_xp_Field.ComponentInterpolationSet(iron.FieldVariableTypes.U,1,
                                      iron.FieldInterpolationTypes.NODE_BASED)

k_xp_Field.CreateFinish()

k_xp_Field.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                        iron.FieldParameterSetTypes.VALUES,
                                        1,k_xp)

	 #--------------------------------------------------------------------------------------------
     # k_x FIELD
     #--------------------------------------------------------------------------------------------

k_x_Field = iron.Field()
k_x_Field.CreateStart(k_x_FieldUserNumber,region)

k_x_Field.TypeSet(iron.FieldTypes.GENERAL)
k_x_Field.MeshDecompositionSet(decomposition)
k_x_Field.GeometricFieldSet(geometricField)
k_x_Field.NumberOfVariablesSet(1)
k_x_Field.VariableTypesSet([iron.FieldVariableTypes.U])
k_x_Field.DataTypeSet(iron.FieldVariableTypes.U,iron.FieldDataTypes.DP)
k_x_Field.DimensionSet(iron.FieldVariableTypes.U,iron.FieldDimensionTypes.SCALAR)
k_x_Field.NumberOfComponentsSet(iron.FieldVariableTypes.U,1)
k_x_Field.VariableLabelSet(iron.FieldVariableTypes.U,"k_x Field")

geometricMeshComponent = geometricField.ComponentMeshComponentGet(iron.FieldVariableTypes.U,1)
k_x_Field.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,geometricMeshComponent)
k_x_Field.ComponentInterpolationSet(iron.FieldVariableTypes.U,1,
                                      iron.FieldInterpolationTypes.NODE_BASED)

k_x_Field.CreateFinish()

k_x_Field.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                        iron.FieldParameterSetTypes.VALUES,
                                        1,k_x)

	 #--------------------------------------------------------------------------------------------
     # k_sp FIELD
     #--------------------------------------------------------------------------------------------

k_sp_Field = iron.Field()
k_sp_Field.CreateStart(k_sp_FieldUserNumber,region)

k_sp_Field.TypeSet(iron.FieldTypes.GENERAL)
k_sp_Field.MeshDecompositionSet(decomposition)
k_sp_Field.GeometricFieldSet(geometricField)
k_sp_Field.NumberOfVariablesSet(1)
k_sp_Field.VariableTypesSet([iron.FieldVariableTypes.U])
k_sp_Field.DataTypeSet(iron.FieldVariableTypes.U,iron.FieldDataTypes.DP)
k_sp_Field.DimensionSet(iron.FieldVariableTypes.U,iron.FieldDimensionTypes.SCALAR)
k_sp_Field.NumberOfComponentsSet(iron.FieldVariableTypes.U,1)
k_sp_Field.VariableLabelSet(iron.FieldVariableTypes.U,"k_sp Field")

geometricMeshComponent = geometricField.ComponentMeshComponentGet(iron.FieldVariableTypes.U,1)
k_sp_Field.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,geometricMeshComponent)
k_sp_Field.ComponentInterpolationSet(iron.FieldVariableTypes.U,1,
                                      iron.FieldInterpolationTypes.NODE_BASED)

k_sp_Field.CreateFinish()

k_sp_Field.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                        iron.FieldParameterSetTypes.VALUES,
                                        1,k_sp)


	 #--------------------------------------------------------------------------------------------
     # k_s FIELD
     #--------------------------------------------------------------------------------------------

k_s_Field = iron.Field()
k_s_Field.CreateStart(k_s_FieldUserNumber,region)

k_s_Field.TypeSet(iron.FieldTypes.GENERAL)
k_s_Field.MeshDecompositionSet(decomposition)
k_s_Field.GeometricFieldSet(geometricField)
k_s_Field.NumberOfVariablesSet(1)
k_s_Field.VariableTypesSet([iron.FieldVariableTypes.U])
k_s_Field.DataTypeSet(iron.FieldVariableTypes.U,iron.FieldDataTypes.DP)
k_s_Field.DimensionSet(iron.FieldVariableTypes.U,iron.FieldDimensionTypes.SCALAR)
k_s_Field.NumberOfComponentsSet(iron.FieldVariableTypes.U,1)
k_s_Field.VariableLabelSet(iron.FieldVariableTypes.U,"k_s Field")

geometricMeshComponent = geometricField.ComponentMeshComponentGet(iron.FieldVariableTypes.U,1)
k_s_Field.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,geometricMeshComponent)
k_s_Field.ComponentInterpolationSet(iron.FieldVariableTypes.U,1,
                                      iron.FieldInterpolationTypes.NODE_BASED)

k_s_Field.CreateFinish()

k_s_Field.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                        iron.FieldParameterSetTypes.VALUES,
                                        1,k_s)

	 #--------------------------------------------------------------------------------------------
     # k_xs FIELD
     #--------------------------------------------------------------------------------------------

k_xs_Field = iron.Field()
k_xs_Field.CreateStart(k_xs_FieldUserNumber,region)

k_xs_Field.TypeSet(iron.FieldTypes.GENERAL)
k_xs_Field.MeshDecompositionSet(decomposition)
k_xs_Field.GeometricFieldSet(geometricField)
k_xs_Field.NumberOfVariablesSet(1)
k_xs_Field.VariableTypesSet([iron.FieldVariableTypes.U])
k_xs_Field.DataTypeSet(iron.FieldVariableTypes.U,iron.FieldDataTypes.DP)
k_xs_Field.DimensionSet(iron.FieldVariableTypes.U,iron.FieldDimensionTypes.SCALAR)
k_xs_Field.NumberOfComponentsSet(iron.FieldVariableTypes.U,1)
k_xs_Field.VariableLabelSet(iron.FieldVariableTypes.U,"k_xs Field")

geometricMeshComponent = geometricField.ComponentMeshComponentGet(iron.FieldVariableTypes.U,1)
k_xs_Field.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,geometricMeshComponent)
k_xs_Field.ComponentInterpolationSet(iron.FieldVariableTypes.U,1,
                                      iron.FieldInterpolationTypes.NODE_BASED)

k_xs_Field.CreateFinish()

k_xs_Field.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                        iron.FieldParameterSetTypes.VALUES,
                                        1,k_xs)

	 #--------------------------------------------------------------------------------------------
     # k_sx FIELD
     #--------------------------------------------------------------------------------------------

k_sx_Field = iron.Field()
k_sx_Field.CreateStart(k_sx_FieldUserNumber,region)

k_sx_Field.TypeSet(iron.FieldTypes.GENERAL)
k_sx_Field.MeshDecompositionSet(decomposition)
k_sx_Field.GeometricFieldSet(geometricField)
k_sx_Field.NumberOfVariablesSet(1)
k_sx_Field.VariableTypesSet([iron.FieldVariableTypes.U])
k_sx_Field.DataTypeSet(iron.FieldVariableTypes.U,iron.FieldDataTypes.DP)
k_sx_Field.DimensionSet(iron.FieldVariableTypes.U,iron.FieldDimensionTypes.SCALAR)
k_sx_Field.NumberOfComponentsSet(iron.FieldVariableTypes.U,1)
k_sx_Field.VariableLabelSet(iron.FieldVariableTypes.U,"k_sx Field")

geometricMeshComponent = geometricField.ComponentMeshComponentGet(iron.FieldVariableTypes.U,1)
k_sx_Field.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,geometricMeshComponent)
k_sx_Field.ComponentInterpolationSet(iron.FieldVariableTypes.U,1,
                                      iron.FieldInterpolationTypes.NODE_BASED)

k_sx_Field.CreateFinish()

k_sx_Field.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                        iron.FieldParameterSetTypes.VALUES,
                                        1,k_sx)

     #--------------------------------------------------------------------------------------------
     # k2 FIELD
     #--------------------------------------------------------------------------------------------

k2_Field = iron.Field()
k2_Field.CreateStart(k2_FieldUserNumber,region)

k2_Field.TypeSet(iron.FieldTypes.GENERAL)
k2_Field.MeshDecompositionSet(decomposition)
k2_Field.GeometricFieldSet(geometricField)
k2_Field.NumberOfVariablesSet(1)
k2_Field.VariableTypesSet([iron.FieldVariableTypes.U])
k2_Field.DataTypeSet(iron.FieldVariableTypes.U,iron.FieldDataTypes.DP)
k2_Field.DimensionSet(iron.FieldVariableTypes.U,iron.FieldDimensionTypes.SCALAR)
k2_Field.NumberOfComponentsSet(iron.FieldVariableTypes.U,1)
k2_Field.VariableLabelSet(iron.FieldVariableTypes.U,"k2 Field")

geometricMeshComponent = geometricField.ComponentMeshComponentGet(iron.FieldVariableTypes.U,1)
k2_Field.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,geometricMeshComponent)
k2_Field.ComponentInterpolationSet(iron.FieldVariableTypes.U,1,
                                      iron.FieldInterpolationTypes.NODE_BASED)

k2_Field.CreateFinish()

k2_Field.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                        iron.FieldParameterSetTypes.VALUES,
                                        1,k2)

     #--------------------------------------------------------------------------------------------
     # k_cisp FIELD
     #--------------------------------------------------------------------------------------------

k_cisp_Field = iron.Field()
k_cisp_Field.CreateStart(k_cisp_FieldUserNumber,region)

k_cisp_Field.TypeSet(iron.FieldTypes.GENERAL)
k_cisp_Field.MeshDecompositionSet(decomposition)
k_cisp_Field.GeometricFieldSet(geometricField)
k_cisp_Field.NumberOfVariablesSet(1)
k_cisp_Field.VariableTypesSet([iron.FieldVariableTypes.U])
k_cisp_Field.DataTypeSet(iron.FieldVariableTypes.U,iron.FieldDataTypes.DP)
k_cisp_Field.DimensionSet(iron.FieldVariableTypes.U,iron.FieldDimensionTypes.SCALAR)
k_cisp_Field.NumberOfComponentsSet(iron.FieldVariableTypes.U,1)
k_cisp_Field.VariableLabelSet(iron.FieldVariableTypes.U,"k_cisp Field")

geometricMeshComponent = geometricField.ComponentMeshComponentGet(iron.FieldVariableTypes.U,1)
k_cisp_Field.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,geometricMeshComponent)
k_cisp_Field.ComponentInterpolationSet(iron.FieldVariableTypes.U,1,
                                      iron.FieldInterpolationTypes.NODE_BASED)

k_cisp_Field.CreateFinish()

k_cisp_Field.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                        iron.FieldParameterSetTypes.VALUES,
                                        1,k_cisp)

     #--------------------------------------------------------------------------------------------
     # k_cis FIELD
     #--------------------------------------------------------------------------------------------

k_cis_Field = iron.Field()
k_cis_Field.CreateStart(k_cis_FieldUserNumber,region)

k_cis_Field.TypeSet(iron.FieldTypes.GENERAL)
k_cis_Field.MeshDecompositionSet(decomposition)
k_cis_Field.GeometricFieldSet(geometricField)
k_cis_Field.NumberOfVariablesSet(1)
k_cis_Field.VariableTypesSet([iron.FieldVariableTypes.U])
k_cis_Field.DataTypeSet(iron.FieldVariableTypes.U,iron.FieldDataTypes.DP)
k_cis_Field.DimensionSet(iron.FieldVariableTypes.U,iron.FieldDimensionTypes.SCALAR)
k_cis_Field.NumberOfComponentsSet(iron.FieldVariableTypes.U,1)
k_cis_Field.VariableLabelSet(iron.FieldVariableTypes.U,"k_cis Field")

geometricMeshComponent = geometricField.ComponentMeshComponentGet(iron.FieldVariableTypes.U,1)
k_cis_Field.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,geometricMeshComponent)
k_cis_Field.ComponentInterpolationSet(iron.FieldVariableTypes.U,1,
                                      iron.FieldInterpolationTypes.NODE_BASED)

k_cis_Field.CreateFinish()

k_cis_Field.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                        iron.FieldParameterSetTypes.VALUES,
                                        1,k_cis)

     #--------------------------------------------------------------------------------------------
     # k_sa FIELD
     #--------------------------------------------------------------------------------------------

k_sa_Field = iron.Field()
k_sa_Field.CreateStart(k_sa_FieldUserNumber,region)

k_sa_Field.TypeSet(iron.FieldTypes.GENERAL)
k_sa_Field.MeshDecompositionSet(decomposition)
k_sa_Field.GeometricFieldSet(geometricField)
k_sa_Field.NumberOfVariablesSet(1)
k_sa_Field.VariableTypesSet([iron.FieldVariableTypes.U])
k_sa_Field.DataTypeSet(iron.FieldVariableTypes.U,iron.FieldDataTypes.DP)
k_sa_Field.DimensionSet(iron.FieldVariableTypes.U,iron.FieldDimensionTypes.SCALAR)
k_sa_Field.NumberOfComponentsSet(iron.FieldVariableTypes.U,1)
k_sa_Field.VariableLabelSet(iron.FieldVariableTypes.U,"k_sa Field")

geometricMeshComponent = geometricField.ComponentMeshComponentGet(iron.FieldVariableTypes.U,1)
k_sa_Field.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,geometricMeshComponent)
k_sa_Field.ComponentInterpolationSet(iron.FieldVariableTypes.U,1,
                                      iron.FieldInterpolationTypes.NODE_BASED)

k_sa_Field.CreateFinish()

k_sa_Field.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                        iron.FieldParameterSetTypes.VALUES,
                                        1,k_sa)

     #--------------------------------------------------------------------------------------------
     # k_cisp_m FIELD
     #--------------------------------------------------------------------------------------------

k_cisp_m_Field = iron.Field()
k_cisp_m_Field.CreateStart(k_cisp_m_FieldUserNumber,region)

k_cisp_m_Field.TypeSet(iron.FieldTypes.GENERAL)
k_cisp_m_Field.MeshDecompositionSet(decomposition)
k_cisp_m_Field.GeometricFieldSet(geometricField)
k_cisp_m_Field.NumberOfVariablesSet(1)
k_cisp_m_Field.VariableTypesSet([iron.FieldVariableTypes.U])
k_cisp_m_Field.DataTypeSet(iron.FieldVariableTypes.U,iron.FieldDataTypes.DP)
k_cisp_m_Field.DimensionSet(iron.FieldVariableTypes.U,iron.FieldDimensionTypes.SCALAR)
k_cisp_m_Field.NumberOfComponentsSet(iron.FieldVariableTypes.U,1)
k_cisp_m_Field.VariableLabelSet(iron.FieldVariableTypes.U,"k_cisp_m Field")

geometricMeshComponent = geometricField.ComponentMeshComponentGet(iron.FieldVariableTypes.U,1)
k_cisp_m_Field.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,geometricMeshComponent)
k_cisp_m_Field.ComponentInterpolationSet(iron.FieldVariableTypes.U,1,
                                      iron.FieldInterpolationTypes.NODE_BASED)

k_cisp_m_Field.CreateFinish()

k_cisp_m_Field.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                        iron.FieldParameterSetTypes.VALUES,
                                        1,k_cisp_m)
                                        
     #--------------------------------------------------------------------------------------------
     # k_cispstar_m FIELD
     #--------------------------------------------------------------------------------------------

k_cispstar_m_Field = iron.Field()
k_cispstar_m_Field.CreateStart(k_cispstar_m_FieldUserNumber,region)

k_cispstar_m_Field.TypeSet(iron.FieldTypes.GENERAL)
k_cispstar_m_Field.MeshDecompositionSet(decomposition)
k_cispstar_m_Field.GeometricFieldSet(geometricField)
k_cispstar_m_Field.NumberOfVariablesSet(1)
k_cispstar_m_Field.VariableTypesSet([iron.FieldVariableTypes.U])
k_cispstar_m_Field.DataTypeSet(iron.FieldVariableTypes.U,iron.FieldDataTypes.DP)
k_cispstar_m_Field.DimensionSet(iron.FieldVariableTypes.U,iron.FieldDimensionTypes.SCALAR)
k_cispstar_m_Field.NumberOfComponentsSet(iron.FieldVariableTypes.U,1)
k_cispstar_m_Field.VariableLabelSet(iron.FieldVariableTypes.U,"k_cispstar_m Field")

geometricMeshComponent = geometricField.ComponentMeshComponentGet(iron.FieldVariableTypes.U,1)
k_cispstar_m_Field.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,geometricMeshComponent)
k_cispstar_m_Field.ComponentInterpolationSet(iron.FieldVariableTypes.U,1,
                                      iron.FieldInterpolationTypes.NODE_BASED)

k_cispstar_m_Field.CreateFinish()

k_cispstar_m_Field.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                        iron.FieldParameterSetTypes.VALUES,
                                        1,k_cispstar_m)
     #--------------------------------------------------------------------------------------------
     # alpha FIELD
     #--------------------------------------------------------------------------------------------

alpha_Field = iron.Field()
alpha_Field.CreateStart(alpha_FieldUserNumber,region)

alpha_Field.TypeSet(iron.FieldTypes.GENERAL)
alpha_Field.MeshDecompositionSet(decomposition)
alpha_Field.GeometricFieldSet(geometricField)
alpha_Field.NumberOfVariablesSet(1)
alpha_Field.VariableTypesSet([iron.FieldVariableTypes.U])
alpha_Field.DataTypeSet(iron.FieldVariableTypes.U,iron.FieldDataTypes.DP)
alpha_Field.DimensionSet(iron.FieldVariableTypes.U,iron.FieldDimensionTypes.SCALAR)
alpha_Field.NumberOfComponentsSet(iron.FieldVariableTypes.U,1)
alpha_Field.VariableLabelSet(iron.FieldVariableTypes.U,"alpha Field")

geometricMeshComponent = geometricField.ComponentMeshComponentGet(iron.FieldVariableTypes.U,1)
alpha_Field.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,geometricMeshComponent)
alpha_Field.ComponentInterpolationSet(iron.FieldVariableTypes.U,1,
                                      iron.FieldInterpolationTypes.NODE_BASED)

alpha_Field.CreateFinish()

alpha_Field.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                        iron.FieldParameterSetTypes.VALUES,
                                        1,alpha)

     #--------------------------------------------------------------------------------------------
     # beta FIELD
     #--------------------------------------------------------------------------------------------

beta_Field = iron.Field()
beta_Field.CreateStart(beta_FieldUserNumber,region)

beta_Field.TypeSet(iron.FieldTypes.GENERAL)
beta_Field.MeshDecompositionSet(decomposition)
beta_Field.GeometricFieldSet(geometricField)
beta_Field.NumberOfVariablesSet(1)
beta_Field.VariableTypesSet([iron.FieldVariableTypes.U])
beta_Field.DataTypeSet(iron.FieldVariableTypes.U,iron.FieldDataTypes.DP)
beta_Field.DimensionSet(iron.FieldVariableTypes.U,iron.FieldDimensionTypes.SCALAR)
beta_Field.NumberOfComponentsSet(iron.FieldVariableTypes.U,1)
beta_Field.VariableLabelSet(iron.FieldVariableTypes.U,"beta Field")

geometricMeshComponent = geometricField.ComponentMeshComponentGet(iron.FieldVariableTypes.U,1)
beta_Field.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,geometricMeshComponent)
beta_Field.ComponentInterpolationSet(iron.FieldVariableTypes.U,1,
                                      iron.FieldInterpolationTypes.NODE_BASED)

beta_Field.CreateFinish()

beta_Field.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                        iron.FieldParameterSetTypes.VALUES,
                                        1,beta)
     #--------------------------------------------------------------------------------------------
     # jon FIELD
     #--------------------------------------------------------------------------------------------

jon_Field = iron.Field()
jon_Field.CreateStart(jon_FieldUserNumber,region)

jon_Field.TypeSet(iron.FieldTypes.GENERAL)
jon_Field.MeshDecompositionSet(decomposition)
jon_Field.GeometricFieldSet(geometricField)
jon_Field.NumberOfVariablesSet(1)
jon_Field.VariableTypesSet([iron.FieldVariableTypes.U])
jon_Field.DataTypeSet(iron.FieldVariableTypes.U,iron.FieldDataTypes.DP)
jon_Field.DimensionSet(iron.FieldVariableTypes.U,iron.FieldDimensionTypes.SCALAR)
jon_Field.NumberOfComponentsSet(iron.FieldVariableTypes.U,1)
jon_Field.VariableLabelSet(iron.FieldVariableTypes.U,"jon Field")

geometricMeshComponent = geometricField.ComponentMeshComponentGet(iron.FieldVariableTypes.U,1)
jon_Field.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,geometricMeshComponent)
jon_Field.ComponentInterpolationSet(iron.FieldVariableTypes.U,1,
                                      iron.FieldInterpolationTypes.NODE_BASED)

jon_Field.CreateFinish()

jon_Field.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                        iron.FieldParameterSetTypes.VALUES,
                                        1,jon)
     #--------------------------------------------------------------------------------------------
     # joff FIELD
     #--------------------------------------------------------------------------------------------

joff_Field = iron.Field()
joff_Field.CreateStart(joff_FieldUserNumber,region)

joff_Field.TypeSet(iron.FieldTypes.GENERAL)
joff_Field.MeshDecompositionSet(decomposition)
joff_Field.GeometricFieldSet(geometricField)
joff_Field.NumberOfVariablesSet(1)
joff_Field.VariableTypesSet([iron.FieldVariableTypes.U])
joff_Field.DataTypeSet(iron.FieldVariableTypes.U,iron.FieldDataTypes.DP)
joff_Field.DimensionSet(iron.FieldVariableTypes.U,iron.FieldDimensionTypes.SCALAR)
joff_Field.NumberOfComponentsSet(iron.FieldVariableTypes.U,1)
joff_Field.VariableLabelSet(iron.FieldVariableTypes.U,"joff Field")

geometricMeshComponent = geometricField.ComponentMeshComponentGet(iron.FieldVariableTypes.U,1)
joff_Field.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,geometricMeshComponent)
joff_Field.ComponentInterpolationSet(iron.FieldVariableTypes.U,1,
                                      iron.FieldInterpolationTypes.NODE_BASED)

joff_Field.CreateFinish()

joff_Field.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                        iron.FieldParameterSetTypes.VALUES,
                                        1,joff)

#-------------------------------------------------------------------------------------------------
# CELLML
#-------------------------------------------------------------------------------------------------

     #--------------------------------------------------------------------------------------------
     # CELLML FIELD
     #--------------------------------------------------------------------------------------------

#Initialise cellml
cellml = iron.CellML()
cellml.CreateStart(cellmlUserNumber,region)


#Importing the cellml model
constantModelIndex = cellml.ModelImport(EcadDynamicsModel)

#Parameters that will not change --> parameters field
cellml.VariableSetAsKnown(constantModelIndex,"membrane/k1")
cellml.VariableSetAsKnown(constantModelIndex,"membrane/k_xp")
cellml.VariableSetAsKnown(constantModelIndex,"membrane/k_x")
cellml.VariableSetAsKnown(constantModelIndex,"membrane/k_sp")
cellml.VariableSetAsKnown(constantModelIndex,"membrane/k_s")
cellml.VariableSetAsKnown(constantModelIndex,"membrane/k_xs")
cellml.VariableSetAsKnown(constantModelIndex,"membrane/k_sx")

cellml.VariableSetAsKnown(constantModelIndex,"membrane/k2")
cellml.VariableSetAsKnown(constantModelIndex,"membrane/k_cisp")
cellml.VariableSetAsKnown(constantModelIndex,"membrane/k_cis")
cellml.VariableSetAsKnown(constantModelIndex,"membrane/k_sa")
cellml.VariableSetAsKnown(constantModelIndex,"membrane/k_cisp_m")
cellml.VariableSetAsKnown(constantModelIndex,"membrane/k_cispstar_m")
cellml.VariableSetAsKnown(constantModelIndex,"membrane/alpha")
cellml.VariableSetAsKnown(constantModelIndex,"membrane/beta")
cellml.VariableSetAsKnown(constantModelIndex,"membrane/jon")
cellml.VariableSetAsKnown(constantModelIndex,"membrane/joff")
#Finish cellml
cellml.CreateFinish()


     #--------------------------------------------------------------------------------------------
     # CELLML OPENCMISS FIELD MAPS
     #--------------------------------------------------------------------------------------------

#Initialise Field Maps
cellml.FieldMapsCreateStart()

#Create the Field Maps


#------------------------#
# Parameters #
#------------------------#

#k1
cellml.CreateFieldToCellMLMap(k1_Field,iron.FieldVariableTypes.U,1,
                              iron.FieldParameterSetTypes.VALUES,constantModelIndex,
                              "membrane/k1",iron.FieldParameterSetTypes.VALUES)
cellml.CreateCellMLToFieldMap(constantModelIndex,"membrane/k1",
                              iron.FieldParameterSetTypes.VALUES,k1_Field,
                              iron.FieldVariableTypes.U,1,iron.FieldParameterSetTypes.VALUES)

#k_xp
cellml.CreateFieldToCellMLMap(k_xp_Field,iron.FieldVariableTypes.U,1,
                              iron.FieldParameterSetTypes.VALUES,constantModelIndex,
                              "membrane/k_xp",iron.FieldParameterSetTypes.VALUES)
cellml.CreateCellMLToFieldMap(constantModelIndex,"membrane/k_xp",
                              iron.FieldParameterSetTypes.VALUES,k_xp_Field,
                              iron.FieldVariableTypes.U,1,iron.FieldParameterSetTypes.VALUES)

#k_x
cellml.CreateFieldToCellMLMap(k_x_Field,iron.FieldVariableTypes.U,1,
                              iron.FieldParameterSetTypes.VALUES,constantModelIndex,
                              "membrane/k_x",iron.FieldParameterSetTypes.VALUES)
cellml.CreateCellMLToFieldMap(constantModelIndex,"membrane/k_x",
                              iron.FieldParameterSetTypes.VALUES,k_x_Field,
                              iron.FieldVariableTypes.U,1,iron.FieldParameterSetTypes.VALUES)

#k_sp
cellml.CreateFieldToCellMLMap(k_sp_Field,iron.FieldVariableTypes.U,1,
                              iron.FieldParameterSetTypes.VALUES,constantModelIndex,
                              "membrane/k_sp",iron.FieldParameterSetTypes.VALUES)
cellml.CreateCellMLToFieldMap(constantModelIndex,"membrane/k_sp",
                              iron.FieldParameterSetTypes.VALUES,k_sp_Field,
                              iron.FieldVariableTypes.U,1,iron.FieldParameterSetTypes.VALUES)

#k_s
cellml.CreateFieldToCellMLMap(k_s_Field,iron.FieldVariableTypes.U,1,
                              iron.FieldParameterSetTypes.VALUES,constantModelIndex,
                              "membrane/k_s",iron.FieldParameterSetTypes.VALUES)
cellml.CreateCellMLToFieldMap(constantModelIndex,"membrane/k_s",
                              iron.FieldParameterSetTypes.VALUES,k_s_Field,
                              iron.FieldVariableTypes.U,1,iron.FieldParameterSetTypes.VALUES)

#k_xs
cellml.CreateFieldToCellMLMap(k_xs_Field,iron.FieldVariableTypes.U,1,
                              iron.FieldParameterSetTypes.VALUES,constantModelIndex,
                              "membrane/k_xs",iron.FieldParameterSetTypes.VALUES)
cellml.CreateCellMLToFieldMap(constantModelIndex,"membrane/k_xs",
                              iron.FieldParameterSetTypes.VALUES,k_xs_Field,
                              iron.FieldVariableTypes.U,1,iron.FieldParameterSetTypes.VALUES)

#k_sx
cellml.CreateFieldToCellMLMap(k_sx_Field,iron.FieldVariableTypes.U,1,
                              iron.FieldParameterSetTypes.VALUES,constantModelIndex,
                              "membrane/k_sx",iron.FieldParameterSetTypes.VALUES)
cellml.CreateCellMLToFieldMap(constantModelIndex,"membrane/k_sx",
                              iron.FieldParameterSetTypes.VALUES,k_sx_Field,
                              iron.FieldVariableTypes.U,1,iron.FieldParameterSetTypes.VALUES)

#k2
cellml.CreateFieldToCellMLMap(k2_Field,iron.FieldVariableTypes.U,1,
                              iron.FieldParameterSetTypes.VALUES,constantModelIndex,
                              "membrane/k2",iron.FieldParameterSetTypes.VALUES)
cellml.CreateCellMLToFieldMap(constantModelIndex,"membrane/k2",
                              iron.FieldParameterSetTypes.VALUES,k2_Field,
                              iron.FieldVariableTypes.U,1,iron.FieldParameterSetTypes.VALUES)

#k_cisp
cellml.CreateFieldToCellMLMap(k_cisp_Field,iron.FieldVariableTypes.U,1,
                              iron.FieldParameterSetTypes.VALUES,constantModelIndex,
                              "membrane/k_cisp",iron.FieldParameterSetTypes.VALUES)
cellml.CreateCellMLToFieldMap(constantModelIndex,"membrane/k_cisp",
                              iron.FieldParameterSetTypes.VALUES,k_cisp_Field,
                              iron.FieldVariableTypes.U,1,iron.FieldParameterSetTypes.VALUES)

#k_cis
cellml.CreateFieldToCellMLMap(k_cis_Field,iron.FieldVariableTypes.U,1,
                              iron.FieldParameterSetTypes.VALUES,constantModelIndex,
                              "membrane/k_cis",iron.FieldParameterSetTypes.VALUES)
cellml.CreateCellMLToFieldMap(constantModelIndex,"membrane/k_cis",
                              iron.FieldParameterSetTypes.VALUES,k_cis_Field,
                              iron.FieldVariableTypes.U,1,iron.FieldParameterSetTypes.VALUES)
#k_sa
cellml.CreateFieldToCellMLMap(k_sa_Field,iron.FieldVariableTypes.U,1,
                              iron.FieldParameterSetTypes.VALUES,constantModelIndex,
                              "membrane/k_sa",iron.FieldParameterSetTypes.VALUES)
cellml.CreateCellMLToFieldMap(constantModelIndex,"membrane/k_sa",
                              iron.FieldParameterSetTypes.VALUES,k_sa_Field,
                              iron.FieldVariableTypes.U,1,iron.FieldParameterSetTypes.VALUES)

#k_cisp_m
cellml.CreateFieldToCellMLMap(k_cisp_m_Field,iron.FieldVariableTypes.U,1,
                              iron.FieldParameterSetTypes.VALUES,constantModelIndex,
                              "membrane/k_cisp_m",iron.FieldParameterSetTypes.VALUES)
cellml.CreateCellMLToFieldMap(constantModelIndex,"membrane/k_cisp_m",
                              iron.FieldParameterSetTypes.VALUES,k_cisp_m_Field,
                              iron.FieldVariableTypes.U,1,iron.FieldParameterSetTypes.VALUES)
#k_cispstar_m
cellml.CreateFieldToCellMLMap(k_cispstar_m_Field,iron.FieldVariableTypes.U,1,
                              iron.FieldParameterSetTypes.VALUES,constantModelIndex,
                              "membrane/k_cispstar_m",iron.FieldParameterSetTypes.VALUES)
cellml.CreateCellMLToFieldMap(constantModelIndex,"membrane/k_cispstar_m",
                              iron.FieldParameterSetTypes.VALUES,k_cispstar_m_Field,
                              iron.FieldVariableTypes.U,1,iron.FieldParameterSetTypes.VALUES)

#alpha
cellml.CreateFieldToCellMLMap(alpha_Field,iron.FieldVariableTypes.U,1,
                              iron.FieldParameterSetTypes.VALUES,constantModelIndex,
                              "membrane/alpha",iron.FieldParameterSetTypes.VALUES)
cellml.CreateCellMLToFieldMap(constantModelIndex,"membrane/alpha",
                              iron.FieldParameterSetTypes.VALUES,alpha_Field,
                              iron.FieldVariableTypes.U,1,iron.FieldParameterSetTypes.VALUES)

#beta
cellml.CreateFieldToCellMLMap(beta_Field,iron.FieldVariableTypes.U,1,
                              iron.FieldParameterSetTypes.VALUES,constantModelIndex,
                              "membrane/beta",iron.FieldParameterSetTypes.VALUES)
cellml.CreateCellMLToFieldMap(constantModelIndex,"membrane/beta",
                              iron.FieldParameterSetTypes.VALUES,beta_Field,
                              iron.FieldVariableTypes.U,1,iron.FieldParameterSetTypes.VALUES)

#jon
cellml.CreateFieldToCellMLMap(jon_Field,iron.FieldVariableTypes.U,1,
                              iron.FieldParameterSetTypes.VALUES,constantModelIndex,
                              "membrane/jon",iron.FieldParameterSetTypes.VALUES)
cellml.CreateCellMLToFieldMap(constantModelIndex,"membrane/jon",
                              iron.FieldParameterSetTypes.VALUES,jon_Field,
                              iron.FieldVariableTypes.U,1,iron.FieldParameterSetTypes.VALUES)

#joff
cellml.CreateFieldToCellMLMap(joff_Field,iron.FieldVariableTypes.U,1,
                              iron.FieldParameterSetTypes.VALUES,constantModelIndex,
                              "membrane/joff",iron.FieldParameterSetTypes.VALUES)
cellml.CreateCellMLToFieldMap(constantModelIndex,"membrane/joff",
                              iron.FieldParameterSetTypes.VALUES,joff_Field,
                              iron.FieldVariableTypes.U,1,iron.FieldParameterSetTypes.VALUES)

#---------------------#
# Dependent Variables #
#---------------------#

#M1
cellml.CreateFieldToCellMLMap(M1_Field,iron.FieldVariableTypes.U,1,
                              iron.FieldParameterSetTypes.VALUES,constantModelIndex,
                              "membrane/M1",iron.FieldParameterSetTypes.VALUES)
cellml.CreateCellMLToFieldMap(constantModelIndex,"membrane/M1",
                              iron.FieldParameterSetTypes.VALUES,M1_Field,
                              iron.FieldVariableTypes.U,1,iron.FieldParameterSetTypes.VALUES)

#M2
cellml.CreateFieldToCellMLMap(M2_Field,iron.FieldVariableTypes.U,1,
                              iron.FieldParameterSetTypes.VALUES,constantModelIndex,
                              "membrane/M2",iron.FieldParameterSetTypes.VALUES)
cellml.CreateCellMLToFieldMap(constantModelIndex,"membrane/M2",
                              iron.FieldParameterSetTypes.VALUES,M2_Field,
                              iron.FieldVariableTypes.U,1,iron.FieldParameterSetTypes.VALUES)
                              
#MM
cellml.CreateFieldToCellMLMap(MM_Field,iron.FieldVariableTypes.U,1,
                              iron.FieldParameterSetTypes.VALUES,constantModelIndex,
                              "membrane/MM",iron.FieldParameterSetTypes.VALUES)
cellml.CreateCellMLToFieldMap(constantModelIndex,"membrane/MM",
                              iron.FieldParameterSetTypes.VALUES,MM_Field,
                              iron.FieldVariableTypes.U,1,iron.FieldParameterSetTypes.VALUES)                             

#X
cellml.CreateFieldToCellMLMap(X_Field,iron.FieldVariableTypes.U,1,
                              iron.FieldParameterSetTypes.VALUES,constantModelIndex,
                              "membrane/X",iron.FieldParameterSetTypes.VALUES)
cellml.CreateCellMLToFieldMap(constantModelIndex,"membrane/X",
                              iron.FieldParameterSetTypes.VALUES,X_Field,
                              iron.FieldVariableTypes.U,1,iron.FieldParameterSetTypes.VALUES)

#S
cellml.CreateFieldToCellMLMap(S_Field,iron.FieldVariableTypes.U,1,
                              iron.FieldParameterSetTypes.VALUES,constantModelIndex,
                              "membrane/S",iron.FieldParameterSetTypes.VALUES)
cellml.CreateCellMLToFieldMap(constantModelIndex,"membrane/S",
                              iron.FieldParameterSetTypes.VALUES,S_Field,
                              iron.FieldVariableTypes.U,1,iron.FieldParameterSetTypes.VALUES)

#SS
cellml.CreateFieldToCellMLMap(SS_Field,iron.FieldVariableTypes.U,1,
                              iron.FieldParameterSetTypes.VALUES,constantModelIndex,
                              "membrane/SS",iron.FieldParameterSetTypes.VALUES)
cellml.CreateCellMLToFieldMap(constantModelIndex,"membrane/SS",
                              iron.FieldParameterSetTypes.VALUES,SS_Field,
                              iron.FieldVariableTypes.U,1,iron.FieldParameterSetTypes.VALUES)

#Cis_star
cellml.CreateFieldToCellMLMap(Cis_star_Field,iron.FieldVariableTypes.U,1,
                              iron.FieldParameterSetTypes.VALUES,constantModelIndex,
                              "membrane/Cis_star",iron.FieldParameterSetTypes.VALUES)
cellml.CreateCellMLToFieldMap(constantModelIndex,"membrane/Cis_star",
                              iron.FieldParameterSetTypes.VALUES,Cis_star_Field,
                              iron.FieldVariableTypes.U,1,iron.FieldParameterSetTypes.VALUES)

#Cis
cellml.CreateFieldToCellMLMap(Cis_Field,iron.FieldVariableTypes.U,1,
                              iron.FieldParameterSetTypes.VALUES,constantModelIndex,
                              "membrane/Cis",iron.FieldParameterSetTypes.VALUES)
cellml.CreateCellMLToFieldMap(constantModelIndex,"membrane/Cis",
                              iron.FieldParameterSetTypes.VALUES,Cis_Field,
                              iron.FieldVariableTypes.U,1,iron.FieldParameterSetTypes.VALUES)


#Finish Field Maps
cellml.FieldMapsCreateFinish()


     #--------------------------------------------------------------------------------------------
     # MODEL FIELD - CELLML
     #--------------------------------------------------------------------------------------------

#Initialise the Model Field
cellmlModelsField = iron.Field()
cellml.ModelsFieldCreateStart(cellmlModelsFieldUserNumber,cellmlModelsField)
cellml.ModelsFieldCreateFinish()


     #--------------------------------------------------------------------------------------------
     # STATE FIELD - CELLML
     #--------------------------------------------------------------------------------------------

#Initialise the State Field
cellmlStateField = iron.Field()
cellml.StateFieldCreateStart(cellmlStateFieldUserNumber,cellmlStateField)
cellml.StateFieldCreateFinish()

     #-------------------------------------------------------------------------------------------------
     # PARAMETERS FIELD - CELLML
     #-------------------------------------------------------------------------------------------------

#Initialise the Parameters Field
cellmlParametersField = iron.Field()
cellml.ParametersFieldCreateStart(cellmlParametersFieldUserNumber,cellmlParametersField)
cellml.ParametersFieldCreateFinish()


#-------------------------------------------------------------------------------------------------
# REACTION DIFFUSION EQUATIONS
#-------------------------------------------------------------------------------------------------

	#-------------------------------------------------------------------------------------------------
	# M1
	#-------------------------------------------------------------------------------------------------

#M1 Equations Set
M1_Equations = iron.Equations()
M1_EquationsSet.EquationsCreateStart(M1_Equations)
M1_Equations.sparsityType = iron.EquationsSparsityTypes.SPARSE

#Equations Outputs
M1_Equations.outputType = iron.EquationsOutputTypes.NONE

#Finish
M1_EquationsSet.EquationsCreateFinish()

	#-------------------------------------------------------------------------------------------------
	# M2
	#-------------------------------------------------------------------------------------------------

#M2 Equations Set
M2_Equations = iron.Equations()
M2_EquationsSet.EquationsCreateStart(M2_Equations)
M2_Equations.sparsityType = iron.EquationsSparsityTypes.SPARSE

#Equations Outputs
M2_Equations.outputType = iron.EquationsOutputTypes.NONE

#Finish
M2_EquationsSet.EquationsCreateFinish()

	#-------------------------------------------------------------------------------------------------
	# MM
	#-------------------------------------------------------------------------------------------------

#MM Equations Set
MM_Equations = iron.Equations()
MM_EquationsSet.EquationsCreateStart(MM_Equations)
MM_Equations.sparsityType = iron.EquationsSparsityTypes.SPARSE

#Equations Outputs
MM_Equations.outputType = iron.EquationsOutputTypes.NONE

#Finish
MM_EquationsSet.EquationsCreateFinish()

	#-------------------------------------------------------------------------------------------------
	# X
	#-------------------------------------------------------------------------------------------------

#X Equations Set
X_Equations = iron.Equations()
X_EquationsSet.EquationsCreateStart(X_Equations)
X_Equations.sparsityType = iron.EquationsSparsityTypes.SPARSE

#Equations Outputs
X_Equations.outputType = iron.EquationsOutputTypes.NONE

#Finish
X_EquationsSet.EquationsCreateFinish()

	#-------------------------------------------------------------------------------------------------
	# S
	#-------------------------------------------------------------------------------------------------

#S Equations Set
S_Equations = iron.Equations()
S_EquationsSet.EquationsCreateStart(S_Equations)
S_Equations.sparsityType = iron.EquationsSparsityTypes.SPARSE

#Equations Outputs
S_Equations.outputType = iron.EquationsOutputTypes.NONE

#Finish
S_EquationsSet.EquationsCreateFinish()

	#-------------------------------------------------------------------------------------------------
	# SS
	#-------------------------------------------------------------------------------------------------

#SS Equations Set
SS_Equations = iron.Equations()
SS_EquationsSet.EquationsCreateStart(SS_Equations)
SS_Equations.sparsityType = iron.EquationsSparsityTypes.SPARSE

#Equations Outputs
SS_Equations.outputType = iron.EquationsOutputTypes.NONE

#Finish
SS_EquationsSet.EquationsCreateFinish()

	#-------------------------------------------------------------------------------------------------
	# Cis_star
	#-------------------------------------------------------------------------------------------------

#Cis_star Equations Set
Cis_star_Equations = iron.Equations()
Cis_star_EquationsSet.EquationsCreateStart(Cis_star_Equations)
Cis_star_Equations.sparsityType = iron.EquationsSparsityTypes.SPARSE

#Equations Outputs
Cis_star_Equations.outputType = iron.EquationsOutputTypes.NONE

#Finish
Cis_star_EquationsSet.EquationsCreateFinish()

	#-------------------------------------------------------------------------------------------------
	# Cis
	#-------------------------------------------------------------------------------------------------

#Cis Equations Set
Cis_Equations = iron.Equations()
Cis_EquationsSet.EquationsCreateStart(Cis_Equations)
Cis_Equations.sparsityType = iron.EquationsSparsityTypes.SPARSE

#Equations Outputs
Cis_Equations.outputType = iron.EquationsOutputTypes.NONE

#Finish
Cis_EquationsSet.EquationsCreateFinish()




# -------------------------------------------------------------------------------------------------
# Define the problem
# -------------------------------------------------------------------------------------------------

Mech_Problem = iron.Problem()
ProblemSpecification = [iron.ProblemClasses.ELASTICITY,
                        iron.ProblemTypes.FINITE_ELASTICITY,
                        iron.ProblemSubtypes.NONE]
Mech_Problem.CreateStart(Mech_problemUserNumber, ProblemSpecification)
Mech_Problem.CreateFinish()

# -------------------------------------------------------------------------------------------------
# Problem control loop
# -------------------------------------------------------------------------------------------------
Mech_Problem.ControlLoopCreateStart()
Mech_ControlLoop = iron.ControlLoop()
Mech_Problem.ControlLoopGet([iron.ControlLoopIdentifiers.NODE], Mech_ControlLoop)
Mech_ControlLoop.MaximumIterationsSet(load_steps)
Mech_ControlLoop.LoadOutputSet(1)
Mech_Problem.ControlLoopCreateFinish()

# -------------------------------------------------------------------------------------------------
# Problem solvers
# -------------------------------------------------------------------------------------------------
Mech_NonLinearSolver = iron.Solver()
Mech_LinearSolver = iron.Solver()
Mech_Problem.SolversCreateStart()
Mech_Problem.SolverGet([iron.ControlLoopIdentifiers.NODE], 1, Mech_NonLinearSolver)
Mech_NonLinearSolver.outputType = iron.SolverOutputTypes.NONE

Mech_NonLinearSolver.NewtonJacobianCalculationTypeSet(iron.JacobianCalculationTypes.FD)

#Mech_NonLinearSolver.NewtonAbsoluteToleranceSet(1e-14)
#Mech_NonLinearSolver.NewtonSolutionToleranceSet(1e-14)
#Mech_NonLinearSolver.NewtonRelativeToleranceSet(1e-14)
Mech_NonLinearSolver.NewtonLinearSolverGet(Mech_LinearSolver)
Mech_LinearSolver.linearType = iron.LinearSolverTypes.DIRECT
# linearSolver.libraryType = iron.SolverLibraries.LAPACK
Mech_Problem.SolversCreateFinish()

# -------------------------------------------------------------------------------------------------
# Solver and solver equations
# -------------------------------------------------------------------------------------------------

Mech_Solver = iron.Solver()
Mech_SolverEquations = iron.SolverEquations()
Mech_Problem.SolverEquationsCreateStart()
Mech_Problem.SolverGet([iron.ControlLoopIdentifiers.NODE], 1, Mech_Solver)

Mech_Solver.SolverEquationsGet(Mech_SolverEquations)
Mech_SolverEquations.sparsityType = iron.SolverEquationsSparsityTypes.SPARSE

EquationsSetIndex = Mech_SolverEquations.EquationsSetAdd(Mech_EquationsSet)
Mech_Problem.SolverEquationsCreateFinish()

# -------------------------------------------------------------------------------------------------
# Boundary conditions
# -------------------------------------------------------------------------------------------------
Mech_BoundaryConditions = iron.BoundaryConditions()
Mech_SolverEquations.BoundaryConditionsCreateStart(Mech_BoundaryConditions)
for i in range(num_nodes):
#    Mech_BoundaryConditions.AddNode(Mech_DependentField, iron.FieldVariableTypes.U, 1, 1, NodeNums[i][0], 1, iron.BoundaryConditionsTypes.FIXED, 0.0)
#    Mech_BoundaryConditions.AddNode(Mech_DependentField, iron.FieldVariableTypes.U, 1, 1, NodeNums[i][0], 2, iron.BoundaryConditionsTypes.FIXED, 0.0)
    if NodeNums[i][0] == 1:
        Mech_BoundaryConditions.AddNode(Mech_DependentField, iron.FieldVariableTypes.DELUDELN, 1, 1, NodeNums[i][0], 1, iron.BoundaryConditionsTypes.FIXED, 0.0)
        Mech_BoundaryConditions.AddNode(Mech_DependentField, iron.FieldVariableTypes.DELUDELN, 1, 1, NodeNums[i][0], 2, iron.BoundaryConditionsTypes.FIXED, 0.0)

Mech_SolverEquations.BoundaryConditionsCreateFinish()

#_________________________________________________________________________________________________
# REACTION DIFFUSION PROBLEM AND SOLVER
#_________________________________________________________________________________________________

#-------------------------------------------------------------------------------------------------
# PROBLEM - REACTION DIFFUSION
#-------------------------------------------------------------------------------------------------
    
# Create Problem
rd_problem = iron.Problem()
rd_problemSpecification = [iron.ProblemClasses.CLASSICAL_FIELD,
                           iron.ProblemTypes.REACTION_DIFFUSION_EQUATION,
                           iron.ProblemSubtypes.CELLML_REAC_INTEG_REAC_DIFF_STRANG_SPLIT]
rd_problem.CreateStart(rd_problemUserNumber, rd_problemSpecification)
rd_problem.CreateFinish()


# Create control loops
rd_problem.ControlLoopCreateStart()
rd_controlLoop = iron.ControlLoop()

rd_problem.ControlLoopGet([iron.ControlLoopIdentifiers.NODE],rd_controlLoop)
# start time, stop time, time increment
rd_controlLoop.TimesSet(startT,startT+Tstep,Tstep)


#Control Loop Outputs
#rd_controlLoop.TimeOutputSet(outputfreq)
#rd_controlLoop.LoadOutputSet(1)
#rd_controlLoop.OutputTypeSet(iron.ControlLoopOutputTypes.PROGRESS)
#rd_controlLoop.OutputTypeSet(iron.ControlLoopOutputTypes.TIMING)

rd_problem.ControlLoopCreateFinish()


#-------------------------------------------------------------------------------------------------
# SOLVER - REACTION DIFFUSION
#-------------------------------------------------------------------------------------------------

#
#    1st Solver --> DAE
#         |
#         v
#    2nd Solver --> Dynamic 
#         |
#         v      
#    3rd Solver --> DAE
#

#Create rd_problem solver for Strang splitting
rd_problem.SolversCreateStart()


#Create first solver --> DAE Solver
rd_solver = iron.Solver()
rd_problem.SolverGet([iron.ControlLoopIdentifiers.NODE],1,rd_solver)
#rd_solver.DAESolverTypeSet(iron.DAESolverTypes.EULER)
rd_solver.DAETimeStepSet(ODE_TIME_STEP)
rd_solver.OutputTypeSet(iron.SolverOutputTypes.NONE)

#Create second solver --> Dynamic solver for parabolic equation
rd_solver = iron.Solver()
rd_linearsolver = iron.Solver()
rd_problem.SolverGet([iron.ControlLoopIdentifiers.NODE],2,rd_solver)

#Set theta - backward vs forward time step parameter
rd_solver.DynamicThetaSet([1.0])

#Set output type
rd_solver.OutputTypeSet(iron.SolverOutputTypes.NONE)

#Obtain dynamic linear solver from the solver
rd_solver.DynamicLinearSolverGet(rd_linearsolver)

#Set Library
#rd_solver.LibraryTypeSet(iron.SolverLibraries.LAPACK)
#rd_solver.LibraryTypeSet(iron.SolverLibraries.CMISS)
#rd_solver.LibraryTypeSet(iron.SolverLibraries.MUMPS)
#rd_solver.LinearTypeSet(iron.LinearSolverTypes.DIRECT)

#rd_solver.LinearDirectTypeSet(iron.DirectLinearSolverTypes.LU)

rd_linearsolver.LinearIterativeMaximumIterationsSet(100000)
#rd_linearsolver.linearIterativeAbsoluteTolerance = 1.0E-12
#rd_linearsolver.linearIterativeRelativeTolerance = 1.0E-12


#Create third solver --> Another DAE Solver
rd_solver = iron.Solver()
rd_problem.SolverGet([iron.ControlLoopIdentifiers.NODE],3,rd_solver)
#rd_solver.DAESolverTypeSet(iron.DAESolverTypes.EULER)
rd_solver.DAETimeStepSet(ODE_TIME_STEP)
rd_solver.OutputTypeSet(iron.SolverOutputTypes.NONE)

#Finish the rd_problem
rd_problem.SolversCreateFinish()



#-------------------------------------------------------------------------------------------------
# SOLVER CELLML EQUATIONS - REACTION DIFFUSION
#-------------------------------------------------------------------------------------------------

#Start rd_solver Cellml Equations
rd_problem.CellMLEquationsCreateStart()

#Create first solver
#cellml equations
rd_solver = iron.Solver()
rd_problem.SolverGet([iron.ControlLoopIdentifiers.NODE],1,rd_solver)
cellmlEquations = iron.CellMLEquations()
rd_solver.CellMLEquationsGet(cellmlEquations)
#Add into Cellml Environment
cellmlIndex = cellmlEquations.CellMLAdd(cellml)

#Create third solver
#cellml equations
rd_solver = iron.Solver()
rd_problem.SolverGet([iron.ControlLoopIdentifiers.NODE],3,rd_solver)
cellmlEquations = iron.CellMLEquations()
rd_solver.CellMLEquationsGet(cellmlEquations)
#Add into Cellml Environment
cellmlIndex = cellmlEquations.CellMLAdd(cellml)

#Finish the solver Cellml Equations
rd_problem.CellMLEquationsCreateFinish()



#-------------------------------------------------------------------------------------------------
# SOLVER EQUATIONS - REACTION DIFFUSION
#-------------------------------------------------------------------------------------------------

#Start solver equations
rd_problem.SolverEquationsCreateStart()

#Create second solver
#Solver equations
rd_solver = iron.Solver()
rd_solverEquations = iron.SolverEquations()
rd_problem.SolverGet([iron.ControlLoopIdentifiers.NODE],2,rd_solver)
rd_solver.SolverEquationsGet(rd_solverEquations)


rd_solverEquations.sparsityType = iron.SolverEquationsSparsityTypes.SPARSE

equationsSetIndex = rd_solverEquations.EquationsSetAdd(M1_EquationsSet)
equationsSetIndex = rd_solverEquations.EquationsSetAdd(M2_EquationsSet)
equationsSetIndex = rd_solverEquations.EquationsSetAdd(MM_EquationsSet)
equationsSetIndex = rd_solverEquations.EquationsSetAdd(X_EquationsSet)
equationsSetIndex = rd_solverEquations.EquationsSetAdd(S_EquationsSet)
equationsSetIndex = rd_solverEquations.EquationsSetAdd(SS_EquationsSet)
equationsSetIndex = rd_solverEquations.EquationsSetAdd(Cis_star_EquationsSet)
equationsSetIndex = rd_solverEquations.EquationsSetAdd(Cis_EquationsSet)

rd_problem.SolverEquationsCreateFinish()



#-------------------------------------------------------------------------------------------------
# BOUNDARY CONDITIONS - REACTION DIFFUSION
#-------------------------------------------------------------------------------------------------

# Set up Boundary Conditions & Perturbation Initial Conditions
rd_boundaryConditions = iron.BoundaryConditions()
rd_solverEquations.BoundaryConditionsCreateStart(rd_boundaryConditions)

rd_condition = iron.BoundaryConditionsTypes.FIXED
rd_value = 0.0

for i in range(num_nodes):
    node = NodeNums[i][0]
    nodeDomain = decomposition.NodeDomainGet(node,1)
    if nodeDomain==computationalNodeNumber:
        #boundary conditions
        if NodeNums[i][1] == 1:
            rd_boundaryConditions.SetNode(M1_Field,iron.FieldVariableTypes.DELUDELN,
                           1,1,node,1,rd_condition,rd_value)
            rd_boundaryConditions.SetNode(M2_Field,iron.FieldVariableTypes.DELUDELN,
                           1,1,node,1,rd_condition,rd_value)
            rd_boundaryConditions.SetNode(MM_Field,iron.FieldVariableTypes.DELUDELN,
                           1,1,node,1,rd_condition,rd_value)
            rd_boundaryConditions.SetNode(X_Field,iron.FieldVariableTypes.DELUDELN,
                           1,1,node,1,rd_condition,rd_value)
            rd_boundaryConditions.SetNode(S_Field,iron.FieldVariableTypes.DELUDELN,
                           1,1,node,1,rd_condition,rd_value)
            rd_boundaryConditions.SetNode(SS_Field,iron.FieldVariableTypes.DELUDELN,
                           1,1,node,1,rd_condition,rd_value)
            rd_boundaryConditions.SetNode(Cis_star_Field,iron.FieldVariableTypes.DELUDELN,
                           1,1,node,1,rd_condition,rd_value)
            rd_boundaryConditions.SetNode(Cis_Field,iron.FieldVariableTypes.DELUDELN,
                           1,1,node,1,rd_condition,rd_value)
rd_solverEquations.BoundaryConditionsCreateFinish()


#_________________________________________________________________________________________________
# TIME LOOP - COUPLING MECHANICS AND REACTION DIFFUSION
#_________________________________________________________________________________________________
f = open("./appliedforce_3.txt","w+")
f.write("# of Nodes: %d # of Boundary Nodes: %d # of Time steps: %d\r\n" % (num_nodes,Num_BoundaryNodes,endT))

# Update k1 with a random number generator
for i in range(num_nodes):
    #Obtain Node Number
    node = NodeNums[i][0]
    k1_Field.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,node,1,random.uniform(k1-0.2,k1+0.2))

#total number of loops
n_loops = int((endT - startT)/Tstep)

#counter to keep track of the current start time
startT_current = startT

# ring_generation_t
ring_generation_t = 600.0

breaked = False 
#Looping through all the time steps
for t in range (0,n_loops):
    
    
    
    
    
    
    
    
    threshold_all = []
    node_index = []
    for i in range(num_nodes):
        #Obtain Node Number
        node = NodeNums[i][0]
        
        #Parallel Programming
        nodeDomain = decomposition.NodeDomainGet(node,1)
        if nodeDomain==computationalNodeNumber:
            radius = math.sqrt(NodeCoords[i][0]**2+NodeCoords[i][1]**2)
            if radius>4.0:

                concentration_X = X_Field.ParameterSetGetNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,node,1)
                concentration_S = S_Field.ParameterSetGetNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,node,1)
                concentration_Cis_star = Cis_star_Field.ParameterSetGetNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,node,1)
                concentration_Cis = Cis_Field.ParameterSetGetNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,node,1)
#                
#                print concentration_S 
#                print concentration_X
                concentration_dimerandclusters = concentration_X+concentration_S+concentration_Cis_star+concentration_Cis
                threshold_singlepoint = (concentration_dimerandclusters/num_nodes)*f_crit
                threshold_all.append(threshold_singlepoint)
                node_index.append(node)
            	f.write("Time: %d Concentration_Cis_star: %f Concentration_Cis: %f Threshold: %f Node Index: %d\r\n" % (t,concentration_Cis_star,concentration_Cis,threshold_singlepoint,node))
    
    force_threshold = sum(threshold_all)
    f.write("Force needed to separate the cell-cell contact: %f\r\n" %force_threshold)	
    threshold_new = sorted(threshold_all)
    if t > aimed_t or t == aimed_t:
        if break_mode == 1:     
            force_mag = threshold_new[-1] + 0.001
        elif break_mode == 2:
            force_mag = threshold_new[len(threshold_all)/2] + 0.001
        elif break_mode == 3:
            force_mag = threshold_new[0] + 0.001
        elif break_mode == 4:
            force_mag = 1.0
    else:
        force_mag = 1.0
	
    for i in range(len(threshold_all)):
    #Parallel Programming
        nodeDomain = decomposition.NodeDomainGet(node,1)
        if nodeDomain==computationalNodeNumber:
            
            concentration_X = X_Field.ParameterSetGetNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,node_index[i],1)
            concentration_S = S_Field.ParameterSetGetNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,node_index[i],1)
            concentration_Cis_star = Cis_star_Field.ParameterSetGetNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,node_index[i],1)
            concentration_Cis = Cis_Field.ParameterSetGetNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,node_index[i],1)
            
            if t > ring_generation_t:
                k_cisp_Field.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,node_index[i],1,0.0)

#    for i in range(num_nodes):
#        node = NodeNums[i][0]          
#        if breaked:
#            f.write('Break = True')
#            
#            concentration_X = X_Field.ParameterSetGetNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,node,1)
#            concentration_S = S_Field.ParameterSetGetNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,node,1)
#            concentration_Cis_star = Cis_star_Field.ParameterSetGetNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,node,1)
#            concentration_Cis = Cis_Field.ParameterSetGetNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,node,1)
#            
#            k1_Field.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,node,1,0.0)
#            k_xs_Field.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,node,1,0.0)
#            k_s_Field.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,node,1,k_s)
#            k_x_Field.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,node,1,k_x)
#            k_sx_Field.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,node,1,0.0)
#            k_xp_Field.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,node,1,concentration_X)
#            k_sp_Field.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,node,1,concentration_S)
#            
#            k2_Field.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,node,1,0.0)
#            k_cis_Field.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,node,1,0.0)
#            k_cisp_Field.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,node,1,0.0)
#            k_cisp_m_Field.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,node,1,concentration_Cis_star)
#            k_cispstar_m_Field.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,node,1,concentration_Cis)

#                    
#-------------------------------------------------------------------------------------------------
# SOLVE - TIME LOOP
#-------------------------------------------------------------------------------------------------        

    #Reaction Diffusion Solver
    rd_problem.Solve()

    #Mechanics Solver
    Mech_Problem.Solve()
    
    #-------------------------------------------------------------------------------------------------
    # CONTROL LOOP UPDATE - TIME LOOP
    #-------------------------------------------------------------------------------------------------
    startT_current = startT_current + Tstep
    
    # start time, stop time, time increment 
    rd_controlLoop.TimesSet(startT_current,startT_current+Tstep,Tstep)

    #-------------------------------------------------------------------------------------------------
    # OUTPUT - TIME LOOP
    #-------------------------------------------------------------------------------------------------
    # Export results
    fields = iron.Fields()
    fields.CreateRegion(region)
	
    # Create new path
    newpath = r"./output_circle_3"
    if not os.path.exists(newpath):
        os.makedirs(newpath)

    #Numbering the files from 0 to n_loops
    num_file = "{0:0" + str(len(str(n_loops))) + "}"
	
    #Output files based on output frequency
    if t % outputfreq == 0:

        #Naming the file properly
        name = "output_circle_3/AdherensJunctions_" + str(t)

        #Exporting the files
        fields.NodesExport(name,"FORTRAN")
        fields.ElementsExport(name,"FORTRAN")
        print "Writing", name, "of", str(n_loops)



    #Only need one element file. Update when elements are compatable with CMGUI
#    if t == 0:
#        fields.ElementsExport("ActinWaves","FORTRAN")
        
    fields.Finalise()


# Finalise OpenCMISS-Iron
iron.Finalise()
f.close()












## -------------------------------------------------------------------------------------------------
## Solver the problem
## -------------------------------------------------------------------------------------------------
#problem.Solve()

## -------------------------------------------------------------------------------------------------
## Export
## -------------------------------------------------------------------------------------------------
#fields = iron.Fields()
#fields.CreateRegion(region)
#name = "cube_2.0/SimpleMigration"
#fields.NodesExport(name, "FORTRAN")
#fields.ElementsExport(name, "FORTRAN")
#fields.Finalise()
#f.close()
