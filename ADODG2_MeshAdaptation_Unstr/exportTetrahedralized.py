#### import the simple module from the paraview
from paraview.simple import *

#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

print("Reading flow2InterpolateFrom.vtu file...")

# create a new 'XML Unstructured Grid Reader'
flow = XMLUnstructuredGridReader(registrationName='flow', FileName=['./flow2InterpolateFrom.vtu'])

# Properties modified on flow_00000vtu
flow.PointArrayStatus = []
flow.TimeArray = 'None'

renderView1 = GetActiveViewOrCreate('RenderView')

UpdatePipeline(time=0.0, proxy=flow)

print("Tetrahedralize the mesh...")

# create a new 'Tetrahedralize'
tetrahedralize1 = Tetrahedralize(registrationName='Tetrahedralize1', Input=flow)

UpdatePipeline(time=0.0, proxy=tetrahedralize1)

# hide data in view
Hide(flow, renderView1)

# update the view to ensure updated data information
renderView1.Update()

# get layout
layout1 = GetLayout()

# split cell
layout1.SplitHorizontal(0, 0.5)

# set active view
SetActiveView(None)

# Create a new 'SpreadSheet View'
spreadSheetView1 = CreateView('SpreadSheetView')
spreadSheetView1.ColumnToSort = ''
spreadSheetView1.BlockSize = 1024

# show data in view
tetrahedralize1Display_1 = Show(tetrahedralize1, spreadSheetView1, 'SpreadSheetRepresentation')

# assign view to a particular cell in the layout
AssignViewToLayout(view=spreadSheetView1, layout=layout1, hint=2)

# Properties modified on spreadSheetView1
spreadSheetView1.FieldAssociation = 'Cell Data'

# Properties modified on spreadSheetView1
spreadSheetView1.GenerateCellConnectivity = 1

print("Export the connectivity of the tetrahedralized mesh...")
# export view
ExportView('PointIndices.csv', view=spreadSheetView1)

print("Deleting everything...")

Delete(spreadSheetView1)
del spreadSheetView1

Delete(tetrahedralize1)
del tetrahedralize1

Delete(flow)
del flow

print("Finished creation of tetrahedralized mesh!")
