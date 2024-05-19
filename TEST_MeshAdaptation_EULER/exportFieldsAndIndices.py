from exportFieldsAndIndices_Params import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()



print("Opening file...")

# create a new 'XML Unstructured Grid Reader'
if options['fileFormat'] == "VTU":
    flow = XMLUnstructuredGridReader(registrationName='flow', FileName=[options['fileNameIn']])
    flow.TimeArray = 'None'
elif options['fileFormat'] == "CGNS":
    flow = CGNSSeriesReader(registrationName='flow', FileNames=[options['fileNameIn']])
else:
    print("ERROR! Format ", options['fileFormat'], " is unknown!")
    exit()

# Properties modified on flow_00000vtu
# flow.PointArrayStatus = ['Gradient_Density_X', 'Gradient_Density_Y', 'Gradient_Pressure_X', 'Gradient_Pressure_Y', 'Gradient_Mach_X', 'Gradient_Mach_Y', 'Gradient_Density_Z', 'Gradient_Pressure_Z', 'Gradient_Mach_Z', 'Eigenvalues_Density_0', 'Eigenvalues_Density_1', 'Eigenvalues_Pressure_0', 'Eigenvalues_Pressure_1', 'Eigenvalues_Mach_0', 'Eigenvalues_Mach_1', 'Hessian_Density_XX', 'Hessian_Density_XY', 'Hessian_Density_XZ', 'Hessian_Density_YX', 'Hessian_Density_YY', 'Hessian_Density_YZ', 'Hessian_Density_ZX', 'Hessian_Density_ZY', 'Hessian_Density_ZZ', 'Eigenvalues_Density_2', 'Hessian_Pressure_XX', 'Hessian_Pressure_XY', 'Hessian_Pressure_XZ', 'Hessian_Pressure_YX', 'Hessian_Pressure_YY', 'Hessian_Pressure_YZ', 'Hessian_Pressure_ZX', 'Hessian_Pressure_ZY', 'Hessian_Pressure_ZZ', 'Eigenvalues_Pressure_2', 'Hessian_Mach_XX', 'Hessian_Mach_XY', 'Hessian_Mach_XZ', 'Hessian_Mach_YX', 'Hessian_Mach_YY', 'Hessian_Mach_YZ', 'Hessian_Mach_ZX', 'Hessian_Mach_ZY', 'Hessian_Mach_ZZ', 'Eigenvalues_Mach_2']


renderView1 = GetActiveViewOrCreate('RenderView')

fields2Load = ['Gradient_'+options['AdaptField']+'_X', 'Gradient_'+options['AdaptField']+'_Y', 'Gradient_'+options['AdaptField']+'_Z']

if options['BLTreatment']:
    fields2Load.append('Wall_Distance')

flow.PointArrayStatus = fields2Load
UpdatePipeline(time=0.0, proxy=flow)

if options['ComputeGradients']:
    fields2Load = [options['AdaptField']]
    if options['BLTreatment']:
        fields2Load.append('Wall_Distance')

    flow.PointArrayStatus = fields2Load
    gradient1 = Gradient(registrationName='Gradient1', Input=flow)
    gradient1.ScalarArray = ['POINTS', options['AdaptField']]
    calculator1 = Calculator(registrationName='Calculator1', Input=gradient1)
    calculator1.ResultArrayName = 'Gradient_'+options['AdaptField']+'_X'
    calculator1.Function = 'Gradient_X'
    calculator2 = Calculator(registrationName='Calculator2', Input=calculator1)
    calculator2.ResultArrayName = 'Gradient_'+options['AdaptField']+'_Y'
    calculator2.Function = 'Gradient_Y'
    ToUse = calculator2
    if options['nDim'] == 3:
        calculator3 = Calculator(registrationName='Calculator3', Input=calculator2)
        calculator3.ResultArrayName = 'Gradient_'+options['AdaptField']+'_Z'
        calculator3.Function = 'Gradient_Z'
        ToUse = calculator3
    flow = ToUse


if options['AdaptField'] == "Pressure":
    calculator1 = Calculator(registrationName='Calculator1', Input=flow)
    calculator1.ResultArrayName = 'Gradient_Pressure_X'
    calculator1.Function = 'Gradient_Pressure_X/'+str(options['PressureAdimFactor'])
    calculator2 = Calculator(registrationName='Calculator2', Input=calculator1)
    calculator2.ResultArrayName = 'Gradient_Pressure_Y'
    calculator2.Function = 'Gradient_Pressure_Y/'+str(options['PressureAdimFactor'])
    ToUse = calculator2
    if options['nDim'] == 3:
        calculator3 = Calculator(registrationName='Calculator3', Input=calculator2)
        calculator3.ResultArrayName = 'Gradient_Pressure_Z'
        calculator3.Function = 'Gradient_Pressure_Z/'+str(options['PressureAdimFactor'])
        calculator2 = calculator3
        ToUse = calculator3
    flow = ToUse

if options['AdaptField'] == "TotalPressure":

    # For now I do not think this work with ComputeGradients set to True

    fields2Load = ['Gradient_Mach_X', 'Gradient_Mach_Y', 'Gradient_Mach_Z', 'Mach', 'Gradient_Pressure_X', 'Gradient_Pressure_Y', 'Gradient_Pressure_Z', 'Pressure']
    if options['BLTreatment']:
        fields2Load.append('Wall_Distance')

    flow.PointArrayStatus = fields2Load
    # Properties modified on calculator7
    UpdatePipeline(time=0.0, proxy=flow)
    calculator1 = Calculator(registrationName='Calculator1', Input=flow)
    calculator1.ResultArrayName = 'InsideParenthesis'
    calculator1.Function = '1+'+str((options['gamma']-1)/2)+'*Mach^2'
    UpdatePipeline(time=1.0, proxy=calculator1)
    calculator2 = Calculator(registrationName='Calculator2', Input=calculator1)
    calculator2.ResultArrayName = 'GradP'
    if options['nDim'] == 2:
        calculator2.Function = '(Gradient_Pressure_X*iHat+Gradient_Pressure_Y*jHat)/'+str(options['PressureAdimFactor'])
    else:
        calculator2.Function = '(Gradient_Pressure_X*iHat+Gradient_Pressure_Y*jHat+Gradient_Pressure_Z*kHat)/'+str(options['PressureAdimFactor'])

    UpdatePipeline(time=1.0, proxy=calculator2)
    calculator3 = Calculator(registrationName='Calculator3', Input=calculator2)
    calculator3.ResultArrayName = 'GradM'
    if options['nDim'] == 2:
        calculator3.Function = 'Gradient_Mach_X*iHat+Gradient_Mach_Y*jHat'
    else:
        calculator3.Function = 'Gradient_Mach_X*iHat+Gradient_Mach_Y*jHat+Gradient_Mach_Z*kHat'

    UpdatePipeline(time=1.0, proxy=calculator3)
    calculator4 = Calculator(registrationName='Calculator4', Input=calculator3)
    calculator4.ResultArrayName = 'FirstDeriv'
    calculator4.Function = 'GradP*(InsideParenthesis)^('+str(options['gamma'])+'/'+str(options['gamma']-1)+')'
    UpdatePipeline(time=1.0, proxy=calculator4)
    calculator5 = Calculator(registrationName='Calculator5', Input=calculator4)
    calculator5.ResultArrayName = 'SecondDeriv'
    calculator5.Function = '((Pressure-'+str(options['PressureInf'])+')/'+str(options['PressureAdimFactor'])+')*('+str(options['gamma'])+'/'+str(options['gamma']-1)+')*((InsideParenthesis)^(1/'+str(options['gamma']-1)+'))*'+str(options['gamma']-1)+'*Mach*GradM'
    UpdatePipeline(time=1.0, proxy=calculator5)
    calculator6 = Calculator(registrationName='Calculator6', Input=calculator5)
    calculator6.ResultArrayName = 'Gradient_TotalPressure'
    calculator6.Function = 'FirstDeriv+SecondDeriv'
    UpdatePipeline(time=1.0, proxy=calculator6)
    SaveState('State.pvsm')

    flow = calculator6


if options['fileFormat'] == "CGNS":
    mergeBlocks1 = MergeBlocks(registrationName='MergeBlocks1', Input=flow)
    flow = mergeBlocks1


if options['fileFormatOut'] == "CSV":
    # update the view to ensure updated data information
    renderView1.Update()

    # get layout
    layout1 = GetLayout()

    # split cell
    layout1.SplitHorizontal(0, 0.5)

    # set active view
    SetActiveView(None)

    print("Creating and saving point Fields spreadhseet...")

    # Create a new 'SpreadSheet View'
    spreadSheetView1 = CreateView('SpreadSheetView')
    spreadSheetView1.ColumnToSort = ''
    spreadSheetView1.BlockSize = 1024

    # show data in view
    flow1Display_1 = Show(flow, spreadSheetView1, 'SpreadSheetRepresentation')

    # assign view to a particular cell in the layout
    AssignViewToLayout(view=spreadSheetView1, layout=layout1, hint=2)

    ExportView(options['fileNameOut'], view=spreadSheetView1, RealNumberPrecision=12)

    print("Creating and saving point IDs spreadhseet...")

    # Properties modified on spreadSheetView1
    spreadSheetView1.FieldAssociation = 'Cell Data'

    # Properties modified on spreadSheetView1
    spreadSheetView1.GenerateCellConnectivity = 1

    # export view
    ExportView('PointIndices.csv', view=spreadSheetView1)

    # destroy spreadSheetView1
    Delete(spreadSheetView1)
    del spreadSheetView1

elif options['fileFormatOut'] == "CGNS":

    SaveData(options['fileNameOut'], proxy=flow)

else:
    print("ERROR! Format ", options['fileFormatOut'], " is unknown!")
    exit()

if AdaptField == "TotalPressure":

    Delete(calculator6)
    del calculator6
    Delete(calculator5)
    del calculator5
    Delete(calculator4)
    del calculator4
    Delete(calculator3)
    del calculator3
    Delete(calculator2)
    del calculator2
    Delete(calculator1)
    del calculator1


# destroy flow_00000vtu
Delete(flow)
del flow
