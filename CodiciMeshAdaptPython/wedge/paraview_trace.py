# trace generated using paraview version 5.11.2
#import paraview
#paraview.compatibility.major = 5
#paraview.compatibility.minor = 11

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'XML Unstructured Grid Reader'
flowvtu = XMLUnstructuredGridReader(registrationName='flow.vtu', FileName=['/Users/matteopaganelli/Desktop/masterThesis/CodiciMeshAdaptPython/wedge/flow.vtu'])
flowvtu.PointArrayStatus = ['Density', 'Momentum', 'Energy', 'Pressure', 'Temperature', 'Mach', 'Pressure_Coefficient', 'Velocity']

# Properties modified on flowvtu
flowvtu.TimeArray = 'None'

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')

# show data in view
flowvtuDisplay = Show(flowvtu, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
flowvtuDisplay.Representation = 'Surface'
flowvtuDisplay.ColorArrayName = [None, '']
flowvtuDisplay.SelectTCoordArray = 'None'
flowvtuDisplay.SelectNormalArray = 'None'
flowvtuDisplay.SelectTangentArray = 'None'
flowvtuDisplay.OSPRayScaleArray = 'Density'
flowvtuDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
flowvtuDisplay.SelectOrientationVectors = 'Momentum'
flowvtuDisplay.ScaleFactor = 0.15000000000000002
flowvtuDisplay.SelectScaleArray = 'Density'
flowvtuDisplay.GlyphType = 'Arrow'
flowvtuDisplay.GlyphTableIndexArray = 'Density'
flowvtuDisplay.GaussianRadius = 0.0075
flowvtuDisplay.SetScaleArray = ['POINTS', 'Density']
flowvtuDisplay.ScaleTransferFunction = 'PiecewiseFunction'
flowvtuDisplay.OpacityArray = ['POINTS', 'Density']
flowvtuDisplay.OpacityTransferFunction = 'PiecewiseFunction'
flowvtuDisplay.DataAxesGrid = 'GridAxesRepresentation'
flowvtuDisplay.PolarAxes = 'PolarAxesRepresentation'
flowvtuDisplay.ScalarOpacityUnitDistance = 0.04446864695711832
flowvtuDisplay.OpacityArrayName = ['POINTS', 'Density']
flowvtuDisplay.SelectInputVectors = ['POINTS', 'Momentum']
flowvtuDisplay.WriteLog = ''

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
flowvtuDisplay.OSPRayScaleFunction.Points = [0.7770341038703918, 0.0, 0.5, 0.0, 1.7447890043258667, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
flowvtuDisplay.ScaleTransferFunction.Points = [0.9592652320861816, 0.0, 0.5, 0.0, 1.6455154418945312, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
flowvtuDisplay.OpacityTransferFunction.Points = [0.9592652320861816, 0.0, 0.5, 0.0, 1.6455154418945312, 1.0, 0.5, 0.0]

# reset view to fit data
renderView1.ResetCamera(False)

# get the material library
materialLibrary1 = GetMaterialLibrary()

# update the view to ensure updated data information
renderView1.Update()

# get layout
layout1 = GetLayout()

# split cell
layout1.SplitVertical(0, 0.5)

# set active view
SetActiveView(None)

# Create a new 'SpreadSheet View'
spreadSheetView1 = CreateView('SpreadSheetView')
spreadSheetView1.ColumnToSort = ''
spreadSheetView1.BlockSize = 1024

# show data in view
flowvtuDisplay_1 = Show(flowvtu, spreadSheetView1, 'SpreadSheetRepresentation')

# assign view to a particular cell in the layout
AssignViewToLayout(view=spreadSheetView1, layout=layout1, hint=2)

# Properties modified on flowvtuDisplay_1
flowvtuDisplay_1.Assembly = ''

# export view
ExportView('/Users/matteopaganelli/Desktop/masterThesis/CodiciMeshAdaptPython/wedge/data.csv', view=spreadSheetView1)

#================================================================
# addendum: following script captures some of the application
# state to faithfully reproduce the visualization during playback
#================================================================

#--------------------------------
# saving layout sizes for layouts

# layout/tab size in pixels
layout1.SetSize(2274, 1187)

#-----------------------------------
# saving camera placements for views

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.75, 0.5, 3.482695094980158]
renderView1.CameraFocalPoint = [0.75, 0.5, 0.0]
renderView1.CameraParallelScale = 0.9013878188659973

#--------------------------------------------
# uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).