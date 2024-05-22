#################################################################
###########################OPTIONS###############################

set nDim 2

set GRBL 1.125
set BLSpacing 1e-6

set NMaxLayers 100
set IsotropicHeight 1.0
set NFullLayers 1
set BLHeight 5e-3

set TRexBCName "Wall"

set BLTreatment 1

set Angle 5

set path "./Meshes/Mesh_"

set iter 4
set iterm1 [expr {$iter - 1}]
# set iterm1 "0_ModifiedLETE"

set path2Mesh $path$iter
set extension ".pw"
set ReferenceMeshName $path$iterm1$extension

set Boundaries2Ref [list "AIRFOIL"]
set Boundaries2Preserve [list "FARFIELD"]

set existingDomainsOrBlocks [list "dom-1"]
set isNSOrEuler [list "NS"]
set MaxEdge [list 4.0]
set NSBoundaries [list "AIRFOIL"]

set OutputFormat "SU2"

#################################################################
#################################################################
