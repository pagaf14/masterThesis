from paraview.simple import *

AdaptField= "Mach"
PressureInf= 64318
PressureAdimFactor= 14421
nDim= 2
gamma= 1.4
ComputeGradients= False
BLTreatment= True


fileNameIn= "flow2Refine.vtu"
fileFormat= "VTU"

# fileFormatOut= "CSV"

fileFormatOut= "CGNS"
fileNameOut= "flow2Refine.cgns"

options = {'fileNameIn': fileNameIn,
           'fileFormat': fileFormat,
           'fileNameOut': fileNameOut,
           'fileFormatOut': fileFormatOut,
           'nDim': nDim,
           'AdaptField': AdaptField,
           'PressureInf': PressureInf,
           'PressureAdimFactor': PressureAdimFactor,
           'gamma': PressureAdimFactor,
           'ComputeGradients': ComputeGradients,
           'BLTreatment': BLTreatment}
