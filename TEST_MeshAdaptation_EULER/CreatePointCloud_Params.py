filenames = ["flow.vtu"]
fileFormatIn= "CGNS"
fileNameIn= "flow2Refine.cgns"
fileNameOut= "PointCloud.dat"

nDim= 2
actualNDim= 2


# General parameters for point cloud and spacing
decay= 0.8
maxSpacing= 4.0
minSpacing= 1e-4

# Check if the complexity ratio paradigm has to be used
UseComplexityRatio= True
ComplexityRatio= 2.0

# Check if metric smoothing has to be performed. Default: simple mean over neighbors
Smoothing= True
UseLogExpSmoothing= True

# Check if all of the mesh points have to be written or only the ones whose spacing
# has been modified
WriteAllPoints= False

# Boundary layer strategy for removing points which are between the surface and
# a user-defined distance
BLTreatment= False
wallDistanceMin= 1e-3

# Set if I want to perform more cycles of refinement by removing
# each time from the pool the points previously refined
# !!!!!!!! Only works with Pointwise_OnlyRef and Pointwise_WithCoarsening !!!!!!!!
refineBySteps= False
nStepRefinement= 2

# Boundaries for the region of adaptation
# Limits= "Box"
# xLimMin= -2
# xLimMax= 2
# yLimMin= -2
# yLimMax= 2.5
# zLimMin= -2
# zLimMax= 2.5
#
# BoundingBoxOptions= {'Limits': Limits,
#                      'Parameters': [xLimMin, xLimMax, yLimMin, yLimMax, zLimMin, zLimMax]}

Limits= "Sphere"
Center_X= 0.5
Center_Y= 0.0
Center_Z= 0.0
Radius= 20

BoundingBoxOptions= {'Limits': Limits,
                     'Parameters': [Center_X, Center_Y, Center_Z, Radius]}



options = {'nDim': nDim,
           'actualNDim': actualNDim,
           'fileNameIn': fileNameIn,
           'fileFormatIn': fileFormatIn,
           'fileNameOut': fileNameOut,
           'UseComplexityRatio': UseComplexityRatio,
           'ComplexityRatio': ComplexityRatio,
           'decay': decay,
           'maxSpacing': maxSpacing,
           'minSpacing': minSpacing,
           'BLTreatment': BLTreatment,
           'wallDistanceMin': wallDistanceMin,
           'refineBySteps': refineBySteps,
           'nStepRefinement': nStepRefinement,
           'WriteAllPoints': WriteAllPoints}


# Adaptation method: Pointwise_OnlyRef, Pointwise_WithCoarsening, Pointwise_Clusters, Re, Re_OnlyRef
AdaptationMethod= "Pointwise_OnlyRef"
# Error computation method: Pointwise, Pointwise_Mod, Re
ErrorComputationMethod= "Pointwise"
AdaptField= "Mach"

# Parameters for adaptation from Pointwise methodology
p= 1

NRef = 2
RefScaling = [1.2, 1.4, 1.6]
RefCoeffs = [0.5, 1.0]
NCoarse = 1
CoarseScaling = [1.2, 1.4]
CoarseCoeffs = [0.5]

# Only for AdaptationMethod = "Pointwise_WithCoarsening"
CoarseningFactor= 1.0
Eps_Coarse= 1.1

AdaptMethodOptions= {'AdaptationMethod': AdaptationMethod,
                     'ErrorComputationMethod': ErrorComputationMethod,
                     'p': p,
                     'AdaptField': AdaptField,
                     'NRef': NRef,
                     'RefScaling': RefScaling,
                     'RefCoeffs': RefCoeffs,
                     'NCoarse': NCoarse,
                     'CoarseScaling': CoarseScaling,
                     'CoarseCoeffs': CoarseCoeffs,
                     'CoarseningFactor': CoarseningFactor,
                     'Eps_Coarse': Eps_Coarse,
                     'Smoothing': Smoothing,
                     'UseLogExpSmoothing': UseLogExpSmoothing}
