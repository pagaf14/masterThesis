from paraview.simple import *

AdaptField= "Pressure"
PressureInf= 64318
PressureAdimFactor= 14421
nDim= 2
gamma= 1.4
ComputeGradients= False
BLTreatment= False


fileNameIn= "surface_flow2Refine.vtu"
fileFormat= "VTU"

# fileFormatOut= "CSV"

fileFormatOut= "CSV"
fileNameOut= "PointFields.csv"

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
