#!/usr/bin/env python

# > Main script
# Add Python bindings directory to PATH
# Intialise OpenCMISS library calling

import re
import numpy as np
# load up the opencmiss library of functions into an object called iron.
from opencmiss.iron import iron
# Intialise OpenCMISS

# Set problem parameters

numberOfXi = 3
numberOfLoadIncrements = 1

coordinateSystemUserNumber = 1
regionUserNumber = 2
basisUserNumber = 3
pressureBasisUserNumber = 4
generatedMeshUserNumber = 5
meshUserNumber = 6
decompositionUserNumber = 7
geometricFieldUserNumber = 8
FibreFieldUserNumber = 9
MaterialsFieldUserNumber = 10
DependentFieldUserNumber = 11
ActiveTensionFieldUserNumber = 12
EquationsSetUserNumber = 13
EquationsSetFieldUserNumber = 14
FieldUserNumber = 15
ProblemUserNumber = 16

# Set all diganostic levels on for testing
# iron.DiagnosticsSetOn(iron.DiagnosticTypes.ALL,[1,2,3,4,5],"Diagnostics",["DOMAIN_MAPPINGS_LOCAL_FROM_GLOBAL_CALCULATE"])

# Get the computational nodes information
numberOfComputationalNodes = iron.ComputationalNumberOfNodesGet()
computationalNodeNumber = iron.ComputationalNodeNumberGet()

# -------------------------------------------------------------------------------------------------
# COORDINATE SYSTEM
# -------------------------------------------------------------------------------------------------

# Three Dimensional Coordinate System
coordinateSystem = iron.CoordinateSystem()
coordinateSystem.CreateStart(coordinateSystemUserNumber)
coordinateSystem.dimension = 3
coordinateSystem.CreateFinish()

# -------------------------------------------------------------------------------------------------
# REGION
# -------------------------------------------------------------------------------------------------

# Start Region
region = iron.Region()
region.CreateStart(regionUserNumber, iron.WorldRegion)
# region.label = "cube"
region.coordinateSystem = coordinateSystem
region.CreateFinish()

# -------------------------------------------------------------------------------------------------
# BASIS
# -------------------------------------------------------------------------------------------------

# Simplex Basis
basis = iron.Basis()
basis.CreateStart(basisUserNumber)
basis.type = iron.BasisTypes.SIMPLEX
basis.numberOfXi = numberOfXi
basis.interpolationXi = [iron.BasisInterpolationSpecifications.QUADRATIC_SIMPLEX] * numberOfXi
# what is quadrature order?
basis.QuadratureOrderSet(2)
# basis.quadratureNumberOfGaussXi = [2]*2
basis.CreateFinish()

# -------------------------------------------------------------------------------------------------
# MESH
# -------------------------------------------------------------------------------------------------

# ___________________________________________________________________________________
#                              Initializing Node File

# Inputing node file
node_file = open('sphere.3.node', 'r')

# Reading the values of the first line
num_nodes, num_coords, num_attributes, boundary_markers = node_file.readline().split()

# Converting the numbers to integers
num_nodes = int(num_nodes)
num_coords = int(num_coords)
num_attributes = int(num_attributes)
boundary_markers = int(boundary_markers)

# Creating variables to store node number & boundary marker
# Creating a list consisting of zeros
NodeNums = [[0 for m in range(1)] for n in range(num_nodes)]

# Creating variable to store x, y and z coordinates
NodeCoords = [[0 for m in range(num_coords)] for n in range(num_nodes)]

# Reading data from nodefiler
for i in range(num_nodes):
    NodeNums[i][0], NodeCoords[i][0], NodeCoords[i][1], NodeCoords[i][2] = node_file.readline().split()

    # Converting from string to int/float
    NodeNums[i][0] = int(NodeNums[i][0])
    NodeCoords[i][0] = float(NodeCoords[i][0])
    NodeCoords[i][1] = float(NodeCoords[i][1])
    NodeCoords[i][2] = float(NodeCoords[i][2])

# Closing the file
node_file.close()

# ___________________________________________________________________________________
#                              Initializing Element File

# Inputing element file
elem_file = open('sphere.3.ele', 'r')

# Reading the values of the first line
num_elements, nodes_per_ele, ele_attributes = elem_file.readline().split()

# Converting the numbers to integers
num_elements = int(num_elements)
nodes_per_ele = int(nodes_per_ele)
ele_attributes = int(ele_attributes)

# Creating variable to store element map
ElemMap = [[0 for x in range(10)] for y in range(num_elements)]
ElemMap_corr = [[0 for x in range(10)] for y in range(num_elements)]

Elemindex = [[0 for m in range(1)] for n in range(num_elements)]
# Reading data from elemfile
for i in range(num_elements):


    Elemindex[i][0], ElemMap[i][0], ElemMap[i][1], ElemMap[i][2], ElemMap[i][3], ElemMap[i][4], ElemMap[i][5], ElemMap[i][6], ElemMap[i][7], ElemMap[i][8], ElemMap[i][9] = elem_file.readline().split()


    # Converting from string to int
    Elemindex[i][0] = int(Elemindex[i][0])
    ElemMap_corr[i][0] = int(ElemMap[i][0])
    ElemMap_corr[i][1] = int(ElemMap[i][3])
    ElemMap_corr[i][2] = int(ElemMap[i][1])
    ElemMap_corr[i][3] = int(ElemMap[i][2])
    ElemMap_corr[i][4] = int(ElemMap[i][5])
    ElemMap_corr[i][5] = int(ElemMap[i][6])
    ElemMap_corr[i][6] = int(ElemMap[i][9])
    ElemMap_corr[i][7] = int(ElemMap[i][8])
    ElemMap_corr[i][8] = int(ElemMap[i][7])
    ElemMap_corr[i][9] = int(ElemMap[i][4])
    
# Closing the file
elem_file.close()

# Initialise Nodes
nodes = iron.Nodes()
nodes.CreateStart(region, num_nodes)
nodes.CreateFinish()

# Initialise Mesh
mesh = iron.Mesh()
mesh.CreateStart(meshUserNumber, region, num_coords)
mesh.NumberOfElementsSet(num_elements)
mesh.NumberOfComponentsSet(1)

# Initialise Elements
meshElements = iron.MeshElements()
meshElements.CreateStart(mesh, 1, basis)
for i in range(num_elements):
    element = Elemindex[i][0]
    meshElements.NodesSet(element, [ElemMap_corr[i][0], ElemMap_corr[i][1], ElemMap_corr[i][2], ElemMap_corr[i][3], ElemMap_corr[i][4], ElemMap_corr[i][5], ElemMap_corr[i][6], ElemMap_corr[i][7], ElemMap_corr[i][8], ElemMap_corr[i][9]])

meshElements.CreateFinish()

# Finilise Mesh
mesh.CreateFinish()

# -------------------------------------------------------------------------------------------------
# MESH DECOMPOSITION
# -------------------------------------------------------------------------------------------------

# Parallelization
decomposition = iron.Decomposition()
decomposition.CreateStart(decompositionUserNumber, mesh)
decomposition.type = iron.DecompositionTypes.CALCULATED
decomposition.numberOfDomains = numberOfComputationalNodes
decomposition.CreateFinish()

# -------------------------------------------------------------------------------------------------
# GEOMETRIC FIELD
# -------------------------------------------------------------------------------------------------

# Geometric Field
geometricField = iron.Field()
geometricField.CreateStart(geometricFieldUserNumber, region)
geometricField.MeshDecompositionSet(decomposition)
geometricField.TypeSet(iron.FieldTypes.GEOMETRIC)

geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 1, 1)
geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 2, 1)
geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 3, 1)
geometricField.CreateFinish()

# Update Geometric Field from customized mesh
for i in range(num_nodes):
    node = NodeNums[i][0]
    nodeDomain = decomposition.NodeDomainGet(node, 1)
    if nodeDomain == computationalNodeNumber:
        nodex = NodeCoords[i][0]
        nodey = NodeCoords[i][1]
        nodez = NodeCoords[i][2]
        geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,
                                                iron.FieldParameterSetTypes.VALUES,
                                                1, 1, node, 1, nodex)
        geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,
                                                iron.FieldParameterSetTypes.VALUES,
                                                1, 1, node, 2, nodey)
        geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,
                                                iron.FieldParameterSetTypes.VALUES,
                                                1, 1, node, 3, nodez)

## Update Geometric Field
#geometricField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,
#                                       iron.FieldParameterSetTypes.VALUES)
#geometricField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,
#                                        iron.FieldParameterSetTypes.VALUES)

## -------------------------------------------------------------------------------------------------
## FIBRE FIELD
## -------------------------------------------------------------------------------------------------

## Setting up the fibre field for reaction diffusion equation set.
## For split problem subtype, fibre field is not used
#FibreField = iron.Field()
#FibreField.CreateStart(FibreFieldUserNumber, region)
#FibreField.TypeSet(iron.FieldTypes.FIBRE)
#FibreField.MeshDecompositionSet(decomposition)

#FibreField.GeometricFieldSet(geometricField)

#FibreField.NumberOfVariablesSet(1)

#FibreField.NumberOfComponentsSet(iron.FieldVariableTypes.U, 3)
#FibreField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 1, 1)
#FibreField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 2, 1)
#FibreField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 3, 1)
#FibreField.VariableLabelSet(iron.FieldVariableTypes.U, "Fibre Field")
#FibreField.CreateFinish()

## -------------------------------------------------------------------------------------------------
## Equations and fields
## -------------------------------------------------------------------------------------------------

## -------------------------------------------------------------------------------------------------
## EQUATIONS SET
## What is the relationship between euqations set and dependent field
## Equations set field put all field together, such as dependent field, material field
## -------------------------------------------------------------------------------------------------

## Equation Set
#EquationsSet = iron.EquationsSet()
#EquationsSetField = iron.Field()
#EquationsSetSpecification = [iron.EquationsSetClasses.ELASTICITY,
#                             iron.EquationsSetTypes.FINITE_ELASTICITY,
#                             iron.EquationsSetSubtypes.GUCCIONE_ACTIVECONTRACTION]


#EquationsSet.CreateStart(EquationsSetUserNumber, region,
#                         FibreField, EquationsSetSpecification,
#                         EquationsSetFieldUserNumber, EquationsSetField)

#EquationsSet.CreateFinish()

## -------------------------------------------------------------------------------------------------
## DEPENDENT FIELD
## -------------------------------------------------------------------------------------------------

## Start and label deformation Field Dependent Field
#dependentField = iron.Field()
#dependentField.CreateStart(DependentFieldUserNumber, region)
#dependentField.TypeSet(iron.FieldTypes.GEOMETRIC_GENERAL)
#dependentField.MeshDecompositionSet(decomposition)
#dependentField.GeometricFieldSet(geometricField)
#dependentField.DependentTypeSet(iron.FieldDependentTypes.DEPENDENT)

#dependentField.NumberOfVariablesSet(2)
#dependentField.NumberOfComponentsSet(iron.FieldVariableTypes.U, 4)
#dependentField.NumberOfComponentsSet(iron.FieldVariableTypes.DELUDELN, 4)
#dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 1, 1)
#dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 2, 1)
#dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 3, 1)
#dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 4, 1)
#dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.DELUDELN, 1, 1)
#dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.DELUDELN, 2, 1)
#dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.DELUDELN, 3, 1)
#dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.DELUDELN, 4, 1)

#dependentField.VariableLabelSet(iron.FieldVariableTypes.U, "deformationField")
#dependentField.VariableLabelSet(iron.FieldVariableTypes.DELUDELN, "reaction force Field")

## Finish deformation Field Dependent Field
#dependentField.CreateFinish()

##dependentField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
##                                           iron.FieldParameterSetTypes.VALUES,
##                                           1, 0.0)
##dependentField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
##                                           iron.FieldParameterSetTypes.VALUES,
##                                           2, 0.0)
##dependentField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
##                                           iron.FieldParameterSetTypes.VALUES,
##                                           3, 0.0)
##dependentField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
##                                           iron.FieldParameterSetTypes.VALUES,
##                                           4, 0.0)

## --------------------------------------------------------------------------------------------
## MATERIAL FIELD - DEFORMATION
## --------------------------------------------------------------------------------------------

#MaterialsField = iron.Field()
#MaterialsField.CreateStart(MaterialsFieldUserNumber, region)
#MaterialsField.TypeSet(iron.FieldTypes.MATERIAL)
#MaterialsField.MeshDecompositionSet(decomposition)
#MaterialsField.GeometricFieldSet(geometricField)
#MaterialsField.NumberOfVariablesSet(1)
#MaterialsField.NumberOfComponentsSet(iron.FieldVariableTypes.U, 4)
#MaterialsField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 1, 1)
#MaterialsField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 2, 1)
#MaterialsField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 3, 1)
#MaterialsField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 4, 1)

#MaterialsField.VariableLabelSet(iron.FieldVariableTypes.U, "deformation Materials Field")
#MaterialsField.CreateFinish()

## Set up the constants in the equation
## Mooney-Rivlin constant C1
#MaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
#                                           iron.FieldParameterSetTypes.VALUES,
#                                           1, 1.0)
## Mooney-Rivlin constant C2
#MaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
#                                           iron.FieldParameterSetTypes.VALUES,
#                                           2, 1.0)
## Mooney-Rivlin constant c3
#MaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
#                                           iron.FieldParameterSetTypes.VALUES,
#                                           3, 1.0)
## Mooney-Rivlin constant c4
#MaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
#                                           iron.FieldParameterSetTypes.VALUES,
#                                           4, 1.0)

## --------------------------------------------------------------------------------------------
## INDEPENDENT FIELD - DEFORMATION
## --------------------------------------------------------------------------------------------
#                                           
#ActiveTensionField = iron.Field()
#EquationsSet.IndependentCreateStart(ActiveTensionFieldUserNumber, ActiveTensionField)          
#ActiveTensionField.VariableLabelSet(iron.FieldVariableTypes.U,"Active contraction")
##ActiveTensionField.MeshDecompositionSet(decomposition)
#ActiveTensionField.CreateFinish()  

#ActiveTensionField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,0.0)
#ActiveTensionField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,2,0.0)         
#ActiveTensionField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,3,0.0)           
#                    
#EquationsSet.DependentCreateStart(DependentFieldUserNumber, dependentField)
#EquationsSet.DependentCreateFinish()

#EquationsSet.MaterialsCreateStart(MaterialsFieldUserNumber, MaterialsField)
#EquationsSet.MaterialsCreateFinish()

## -------------------------------------------------------------------------------------------------
## Define the equation, problem, control loop, solver
## -------------------------------------------------------------------------------------------------

## -------------------------------------------------------------------------------------------------
## Create equations
## -------------------------------------------------------------------------------------------------
#equations = iron.Equations()
#EquationsSet.EquationsCreateStart(equations)
#equations.sparsityType = iron.EquationsSparsityTypes.SPARSE
#equations.outputType = iron.EquationsOutputTypes.NONE
#EquationsSet.EquationsCreateFinish()

## what does this part do?
#iron.Field.ParametersToFieldParametersComponentCopy(
#    geometricField, iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1,
#    dependentField, iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1)

#iron.Field.ParametersToFieldParametersComponentCopy(
#    geometricField, iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 2,
#    dependentField, iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 2)

#iron.Field.ParametersToFieldParametersComponentCopy(
#    geometricField, iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 3,
#    dependentField, iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 3)

## -------------------------------------------------------------------------------------------------
## Define the problem
## -------------------------------------------------------------------------------------------------

#problem = iron.Problem()
#problemSpecification = [iron.ProblemClasses.ELASTICITY,
#                        iron.ProblemTypes.FINITE_ELASTICITY,
#                        iron.ProblemSubtypes.NONE]
#problem.CreateStart(ProblemUserNumber, problemSpecification)
#problem.CreateFinish()

## -------------------------------------------------------------------------------------------------
## Problem control loop
## -------------------------------------------------------------------------------------------------
#problem.ControlLoopCreateStart()
#controlLoop = iron.ControlLoop()
#problem.ControlLoopGet([iron.ControlLoopIdentifiers.NODE], controlLoop)
#controlLoop.MaximumIterationsSet(numberOfLoadIncrements)
#controlLoop.LoadOutputSet(1)
#problem.ControlLoopCreateFinish()

## -------------------------------------------------------------------------------------------------
## Problem solvers
## -------------------------------------------------------------------------------------------------
#nonLinearSolver = iron.Solver()
#linearSolver = iron.Solver()
#problem.SolversCreateStart()
#problem.SolverGet([iron.ControlLoopIdentifiers.NODE], 1, nonLinearSolver)
#nonLinearSolver.outputType = iron.SolverOutputTypes.PROGRESS

#nonLinearSolver.NewtonJacobianCalculationTypeSet(iron.JacobianCalculationTypes.FD)

##nonLinearSolver.NewtonAbsoluteToleranceSet(1e-14)
##nonLinearSolver.NewtonSolutionToleranceSet(1e-14)
##nonLinearSolver.NewtonRelativeToleranceSet(1e-14)
#nonLinearSolver.NewtonLinearSolverGet(linearSolver)
#linearSolver.linearType = iron.LinearSolverTypes.DIRECT
## linearSolver.libraryType = iron.SolverLibraries.LAPACK
#problem.SolversCreateFinish()

## -------------------------------------------------------------------------------------------------
## Solver and solver equations
## -------------------------------------------------------------------------------------------------

#solver = iron.Solver()
#solverEquations = iron.SolverEquations()
#problem.SolverEquationsCreateStart()
#problem.SolverGet([iron.ControlLoopIdentifiers.NODE], 1, solver)

#solver.SolverEquationsGet(solverEquations)
#solverEquations.sparsityType = iron.SolverEquationsSparsityTypes.SPARSE

#EquationsSetIndex = solverEquations.EquationsSetAdd(EquationsSet)
#problem.SolverEquationsCreateFinish()

##f = open("cube_2.0/force.txt","w+")

## -------------------------------------------------------------------------------------------------
## Boundary conditions
## -------------------------------------------------------------------------------------------------
#boundaryConditions = iron.BoundaryConditions()
#solverEquations.BoundaryConditionsCreateStart(boundaryConditions)
#boundTol = 1.0E-6
#for i in range(num_nodes):
#    if NodeCoords[i][0]**2+NodeCoords[i][1]**2+NodeCoords[i][2]**2 > 9:
#        ActiveTensionField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,0.0)
#        ActiveTensionField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,2,0.0)
#        ActiveTensionField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,3,0.0)
#		
#    elif NodeCoords[i][2] < 10*np.cos(13/16*np.pi)+boundTol:
#    	
#        boundaryConditions.AddNode(dependentField, iron.FieldVariableTypes.U,
#                                   1, 1, NodeNums[i][0], 3, iron.BoundaryConditionsTypes.FIXED, 0.0)

#solverEquations.BoundaryConditionsCreateFinish()
## -------------------------------------------------------------------------------------------------
## Solver the problem
## -------------------------------------------------------------------------------------------------
#problem.Solve()

## -------------------------------------------------------------------------------------------------
## Export
## -------------------------------------------------------------------------------------------------
#fields = iron.Fields()
#fields.CreateRegion(region)
#name = "./sphere"
#fields.NodesExport(name, "FORTRAN")
#fields.ElementsExport(name, "FORTRAN")
#fields.Finalise()
