#### import the simple module from the paraview
from paraview.simple import *

#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'XML Unstructured Grid Reader'
flow = XMLUnstructuredGridReader(registrationName='flow', FileName=['./flow2InterpolateFrom.vtu'])

# Properties modified on flow_00000vtu
flow.PointArrayStatus = []
flow.TimeArray = 'None'

UpdatePipeline(time=0.0, proxy=flow)

# create a new 'Tetrahedralize'
tetrahedralize1 = Tetrahedralize(registrationName='Tetrahedralize1', Input=flow)

UpdatePipeline(time=0.0, proxy=tetrahedralize1)

# save data
SaveData('./flow2InterpolateFrom_Triangulated.vtu', proxy=tetrahedralize1, ChooseArraysToWrite=0)

# destroy tetrahedralize1
Delete(tetrahedralize1)
del tetrahedralize1

# destroy flow_00000vtu
Delete(flow)
del flow
