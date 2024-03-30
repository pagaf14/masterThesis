# trace generated using paraview version 5.11.1
import paraview
paraview.compatibility.major = 5
paraview.compatibility.minor = 11

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'XML Unstructured Grid Reader'
flowvtu = XMLUnstructuredGridReader(registrationName='flow.vtu', FileName=['/Users/matteopaganelli/Desktop/masterThesis/CodiciMeshAdaptPython/wedge_2/flow.vtu'])
flowvtu.PointArrayStatus = ['Density', 'Momentum', 'Energy', 'Gradient_Density_X', 'Gradient_Density_Y', 'Gradient_Pressure_X', 'Gradient_Pressure_Y', 'Gradient_Mach_X', 'Gradient_Mach_Y', 'Hessian_Density_XX', 'Hessian_Density_XY', 'Hessian_Density_YX', 'Hessian_Density_YY', 'Hessian_Pressure_XX', 'Hessian_Pressure_XY', 'Hessian_Pressure_YX', 'Hessian_Pressure_YY', 'Hessian_Mach_XX', 'Hessian_Mach_XY', 'Hessian_Mach_YX', 'Hessian_Mach_YY', 'Eigenvalues_Density_0', 'Eigenvalues_Density_1', 'Eigenvalues_Pressure_0', 'Eigenvalues_Pressure_1', 'Eigenvalues_Mach_0', 'Eigenvalues_Mach_1']

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
flowvtuDisplay.ScalarOpacityUnitDistance = 0.056446048012770124
flowvtuDisplay.OpacityArrayName = ['POINTS', 'Density']
flowvtuDisplay.add_attribute("SelectInputVectors",['POINTS', 'Momentum'])
flowvtuDisplay.add_attribute("WriteLog" , '')

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
flowvtuDisplay.OSPRayScaleFunction.Points = [0.4575377106666565, 0.0, 0.5, 0.0, 2.800299882888794, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
flowvtuDisplay.ScaleTransferFunction.Points = [0.9665359258651733, 0.0, 0.5, 0.0, 1.615423560142517, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
flowvtuDisplay.OpacityTransferFunction.Points = [0.9665359258651733, 0.0, 0.5, 0.0, 1.615423560142517, 1.0, 0.5, 0.0]

# reset view to fit data
renderView1.ResetCamera(False)

#changing interaction mode based on data extents
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.75, 0.5, 5.025]
renderView1.CameraFocalPoint = [0.75, 0.5, 0.0]

# get the material library
materialLibrary1 = GetMaterialLibrary()

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
flowvtuDisplay_1 = Show(flowvtu, spreadSheetView1, 'SpreadSheetRepresentation')

# assign view to a particular cell in the layout
AssignViewToLayout(view=spreadSheetView1, layout=layout1, hint=2)

# Properties modified on flowvtuDisplay_1
flowvtuDisplay_1.Assembly = ''

# export view
ExportView('/Users/matteopaganelli/Desktop/masterThesis/CodiciMeshAdaptPython/wedge_2/prova.csv', view=spreadSheetView1, RealNumberNotation='Scientific')

#================================================================
# addendum: following script captures some of the application
# state to faithfully reproduce the visualization during playback
#================================================================

#--------------------------------
# saving layout sizes for layouts

# layout/tab size in pixels
layout1.SetSize(1196, 760)

#-----------------------------------
# saving camera placements for views

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.75, 0.5, 5.025]
renderView1.CameraFocalPoint = [0.75, 0.5, 0.0]
renderView1.CameraParallelScale = 0.9013878188659973

#--------------------------------------------
# uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).