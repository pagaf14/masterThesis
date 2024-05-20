from CreatePointCloud_Fun import *
from CreatePointCloud_Params import options
from CreatePointCloud_Params import AdaptMethodOptions
from CreatePointCloud_Params import BoundingBoxOptions

AdaptField = AdaptMethodOptions['AdaptField']
AdaptationMethod = AdaptMethodOptions['AdaptationMethod']
ErrorComputationMethod = AdaptMethodOptions['ErrorComputationMethod']

startAll = time.time()


print("Reading point coordinates...")
start = time.time()

if options['fileFormatIn'] == "CSV":

    # -------------------------------------------------
    with open(options['fileNameIn']) as f:
        first_line = f.readline()
        first_line = np.array(first_line.split(','))
        first_line[-1] = first_line[-1].split('\n')[0]
        index = np.where(first_line == 'Points_0')[0][0]
        x = np.loadtxt(options['fileNameIn'], delimiter=',', skiprows=1, dtype=np.float64)[:,index]
        index = np.where(first_line == 'Points_1')[0][0]
        y = np.loadtxt(options['fileNameIn'], delimiter=',', skiprows=1, dtype=np.float64)[:,index]
        z = np.zeros((len(x),))
        if options['actualNDim'] == 3:
            index = np.where(first_line == 'Points_2')[0][0]
            z = np.loadtxt(options['fileNameIn'], delimiter=',', skiprows=1, dtype=np.float64)[:,index]

        grad_x = []
        grad_y = []
        grad_z = np.zeros((len(x),))
        print("Extracting Gradient of ", AdaptField, " in the x-direction...")
        index = np.where(np.logical_or(first_line == "Gradient_"+AdaptField+"_X", first_line == "Gradient_"+AdaptField+"_0"))[0][0]
        grad_x = np.loadtxt(options['fileNameIn'], delimiter=',', skiprows=1, dtype=np.float64)[:,index]
        print("Extracting Gradient of ", AdaptField, " in the y-direction...")
        index = np.where(np.logical_or(first_line == "Gradient_"+AdaptField+"_Y", first_line == "Gradient_"+AdaptField+"_1"))[0][0]
        grad_y = np.loadtxt(options['fileNameIn'], delimiter=',', skiprows=1, dtype=np.float64)[:,index]

        if options['actualNDim'] == 3:
            print("Extracting Gradient of ", AdaptField, " in the z-direction...")
            index = np.where(np.logical_or(first_line == "Gradient_"+AdaptField+"_Z", first_line == "Gradient_"+AdaptField+"_2"))[0][0]
            grad_z = np.loadtxt(options['fileNameIn'], delimiter=',', skiprows=1, dtype=np.float64)[:,index]

        grad = np.transpose(np.stack((grad_x, grad_y, grad_z)))

        wall_distance = np.zeros((len(grad_y),))
        if options['BLTreatment']:
            print("Reading wall distance field...")
            start = time.time()

            index = np.where(first_line == "Wall_Distance")[0][0]
            wall_distance = np.loadtxt(options['fileNameIn'], delimiter=',', skiprows=1, dtype=np.float64)[:,index]

            print("Done! Elapsed time = ", time.time()-start, "s")

    # exit()
    # -------------------------------------------------
    # x = np.array(PointFieldCoords['Points_0'])

    coords = np.transpose(np.stack((x, y, z)))

    end = time.time()
    print("Done! Elapsed time", end-start,"s")
    # print(coords)

    # Now get all of the cell centroids
    print("Reading point indices...")
    start = time.time()
    PointIDs = np.loadtxt('PointIndices.csv', delimiter=',', skiprows=1, dtype=int)[:,2:]
    CellTypes = np.loadtxt('PointIndices.csv', delimiter=',', skiprows=1, dtype=int)[:,1]
    print("Done! Elapsed time = ", time.time()-start, "s")

elif options['fileFormatIn'] == "CGNS":

    f = h5py.File(options['fileNameIn'],'r')
    field2Read = 'Zone 1'
    # if options['nDim'] == 3:
    #     field2Read = 'blk-1'

    print(list(f['Base'][field2Read].keys()))

    print("Reading point indices...")
    start = time.time()

    listOfKeys = list(f['Base'][field2Read].keys())
    Trias = []
    Quads = []
    Hexas = []
    Tetras = []
    Wedges = []
    Pyras = []

    for i in range(len(listOfKeys)):
        if listOfKeys[i].split('_')[0] == 'Elem':
            field = np.array(f['Base'][field2Read][listOfKeys[i]]['ElementConnectivity'][' data'], dtype=int)-1
            elemName = listOfKeys[i].split('_')[1]
            if elemName == 'Triangles':
                Trias = field.reshape((int(len(field)/3),3))
            elif elemName == 'Quads':
                Quads = field.reshape((int(len(field)/4),4))
            elif elemName == 'Tetras':
                Tetras = field.reshape((int(len(field)/4),4))
            elif elemName == 'Hexas':
                Hexas = field.reshape((int(len(field)/8),8))
            elif elemName == 'Wedges':
                Wedges = field.reshape((int(len(field)/6),6))
            elif elemName == 'Pyramids':
                Pyras = field.reshape((int(len(field)/5),5))



    PointIDs = []
    CellTypes = []

    if options['nDim'] == 2:

        if len(Trias) > 0:
            PointIDs = Trias
            CellTypes = 5*np.ones((len(Trias),), dtype=int)

        if len(Quads) > 0:
            cells2Add = 9*np.ones((len(Quads),), dtype=int)
            if len(Trias) > 0:
                PointIDs = np.append(PointIDs, np.ones((len(Trias),1), dtype=int), axis=1)
                PointIDs = np.append(PointIDs, Quads, axis=0)
                CellTypes = np.append(CellTypes, cells2Add, axis=0)
            else:
                PointIDs = Quads
                CellTypes = cells2Add
    else:
        nPoints = 4
        if len(Pyras) > 0:
            nPoints = 5
        if len(Wedges) > 0:
            nPoints = 6
        if len(Hexas) > 0:
            nPoints = 8

        if len(Tetras) > 0:
            PointIDs = np.append(Tetras, np.ones((len(Tetras), nPoints-4), dtype=int),axis=1)
            CellTypes = 10*np.ones((len(Tetras),), dtype=int)

        if len(Pyras) > 0:
            vec2Add = np.append(Pyras, np.ones((len(Pyras), nPoints-5), dtype=int),axis=1)
            cell2Add = 14*np.ones((len(Pyras),), dtype=int)
            if len(PointIDs) > 0:
                PointIDs = np.append(PointIDs, vec2Add, axis=0)
                CellTypes = np.append(CellTypes, cell2Add, axis=0)
            else:
                PointIDs = vec2Add
                CellTypes = cell2Add

        if len(Wedges) > 0:
            vec2Add = np.append(Wedges, np.ones((len(Wedges), nPoints-6), dtype=int),axis=1)
            cell2Add = 13*np.ones((len(Wedges),), dtype=int)
            if len(PointIDs) > 0:
                PointIDs = np.append(PointIDs, vec2Add, axis=0)
                CellTypes = np.append(CellTypes, cell2Add, axis=0)
            else:
                PointIDs = vec2Add
                CellTypes = cell2Add

        if len(Hexas) > 0:
            vec2Add = Hexas
            cell2Add = 12*np.ones((len(Hexas),), dtype=int)
            if len(PointIDs) > 0:
                PointIDs = np.append(PointIDs, vec2Add, axis=0)
                CellTypes = np.append(CellTypes, cell2Add, axis=0)
            else:
                PointIDs = vec2Add
                CellTypes = cell2Add

    print("Done! Elapsed time = ", time.time()-start, "s")

    print("Reading gradient of adaptation field...")
    start = time.time()
    x = np.array(f['Base'][field2Read]['GridCoordinates']['CoordinateX'][' data'])
    y = np.array(f['Base'][field2Read]['GridCoordinates']['CoordinateY'][' data'])
    z = np.zeros((len(y),))
    if options['actualNDim'] == 3:
        z = np.array(f['Base'][field2Read]['GridCoordinates']['CoordinateZ'][' data'])


    coords = np.transpose(np.stack((x, y, z)))

    print("Done! Elapsed time = ", time.time()-start, "s")


    print("Reading gradient of adaptation field...")
    start = time.time()

    fields = np.array(list(f['Base'][field2Read]['PointData'].keys()))
    index = np.where(np.logical_or(fields=='Gradient_'+AdaptField+'_X', fields=='Gradient_'+AdaptField+'_0'))[0]
    grad_x = np.array(f['Base'][field2Read]['PointData'][fields[index][0]][' data'])
    index = np.where(np.logical_or(fields=='Gradient_'+AdaptField+'_Y', fields=='Gradient_'+AdaptField+'_1'))[0]
    grad_y = np.array(f['Base'][field2Read]['PointData'][fields[index][0]][' data'])
    grad_z = np.zeros((len(grad_y),))
    if options['actualNDim'] == 3:
        index = np.where(np.logical_or(fields=='Gradient_'+AdaptField+'_Z', fields=='Gradient_'+AdaptField+'_2'))[0]
        grad_z = np.array(f['Base'][field2Read]['PointData'][fields[index][0]][' data'])


    grad = np.transpose(np.stack((grad_x, grad_y, grad_z)))

    print("Done! Elapsed time = ", time.time()-start, "s")

    wall_distance = np.zeros((len(grad_y),))
    if options['BLTreatment']:
        print("Reading wall distance field...")
        start = time.time()

        wall_distance = np.array(f['Base'][field2Read]['PointData']['Wall_Distance'][' data'])

        print("Done! Elapsed time = ", time.time()-start, "s")

    # exit()

else:
    print("ERROR! Format ", options['fileFormatIn'], " is unknown!")
    exit(2)

# print(cell_connectivity_matrix)

print("Compute cell volumes/areas of all of the cells...")
AllCellVolumesOrAreas = computeCellVolumesOrAreas(PointIDs, CellTypes, coords, options)

# Compute neighbors for each point
print("Extracting Neighbors points...")
Neighbors, NNeighbors, whereNeighbors, Cells, NCells = getNeighborsPoints(PointIDs, CellTypes, len(coords))

# npNeighbors = np.array(Neighbors, dtype="object")


avg_S = np.zeros((len(coords),))
tot_S = np.zeros((len(coords),))
avg_Length = np.zeros((len(coords),))
tot_dist = np.zeros((len(coords),))
# xVec = [0 for j in range(len(coords))]
# yVec = [0 for j in range(len(coords))]
EdgeLengthMaxS = np.zeros((len(coords),))
isInsideTheBox = False*np.ones((len(coords),))

timeCheckInside = 0.0
timeComputeInitialMetric = 0.0

start = time.time()

isInsideTheBox = []

print("Computing if a point is inside the prescribed bounding box...")

if BoundingBoxOptions['Limits'] == "Box":
    xLimMin = BoundingBoxOptions['Parameters'][0]
    yLimMin = BoundingBoxOptions['Parameters'][2]
    zLimMin = BoundingBoxOptions['Parameters'][4]

    xLimMax = BoundingBoxOptions['Parameters'][1]
    yLimMax = BoundingBoxOptions['Parameters'][3]
    zLimMax = BoundingBoxOptions['Parameters'][5]

    limitsMin = np.array([xLimMin, yLimMin, zLimMin])
    limitsMax = np.array([xLimMax, yLimMax, zLimMax])
    isInsideTheBox = np.all(np.logical_and(coords <= np.tile(np.transpose(limitsMax), (len(coords), 1)), coords >= np.tile(np.transpose(limitsMin), (len(coords), 1))), axis=1)

elif BoundingBoxOptions['Limits'] == "Sphere":

    Center_X = BoundingBoxOptions['Parameters'][0]
    Center_Y = BoundingBoxOptions['Parameters'][1]
    Center_Z = BoundingBoxOptions['Parameters'][2]
    Radius = BoundingBoxOptions['Parameters'][3]

    CenterCoords = np.tile(np.transpose(np.array([Center_X, Center_Y, Center_Z])), (len(coords), 1))
    isInsideTheBox = np.sum( (coords-CenterCoords)*(coords-CenterCoords), axis = 1) <= Radius*Radius*np.ones((len(coords),))
else:
    print("Bounding limit",BoundingBoxOptions['Limits'], "is not implemented!")
    exit(2)

indices = np.array(range(len(coords)))
ActualIndices2Analyze = indices[isInsideTheBox]
ActualIndices2NotAnalyze = indices[np.logical_not(isInsideTheBox)]

end = time.time()

timeCheckInside+= end-start

print("Done! ", len(ActualIndices2Analyze), " found inside the box, ", len(ActualIndices2NotAnalyze), " found outside. Elapsed time = ", timeCheckInside, " s")

if options['BLTreatment']:
    print("Checking which points are farther from the wall than the minimum distance set...")

    isOutsideBL = wall_distance > options['wallDistanceMin']
    ActualIndices2Analyze = indices[np.logical_and(isInsideTheBox, isOutsideBL)]
    ActualIndices2NotAnalyze = indices[np.logical_not(np.logical_and(isInsideTheBox, isOutsideBL))]

    timeCheckInside+= end-start

    print("Done! ", len(ActualIndices2Analyze), " found outside the BL and inside the box, ", len(ActualIndices2NotAnalyze), " found inside the BL. Elapsed time = ", timeCheckInside, " s")


print("Computing initial metric...")
start = time.time()

S_thresh = 0.0
S_Std = 0.0

# For every point compute distance from its neighbors
coordsNeigh = coords[Neighbors, :]
diffCoords = coordsNeigh - np.repeat(coords[:, np.newaxis, :], len(Neighbors[0]), axis=1)
dist = np.linalg.norm(diffCoords, axis=2)

if ErrorComputationMethod == "Pointwise" or ErrorComputationMethod == "Pointwise_Mod":

    h = diffCoords[ActualIndices2Analyze]/np.repeat(dist[ActualIndices2Analyze, :, np.newaxis], 3, axis=2)

    delta_grad = np.repeat(grad[ActualIndices2Analyze, np.newaxis, :], len(Neighbors[0]), axis=1)
    if ErrorComputationMethod == "Pointwise":
        delta_grad -= grad[Neighbors[ActualIndices2Analyze], :]

    Sij = (dist[ActualIndices2Analyze]**AdaptMethodOptions['p']) * np.absolute(np.sum(h*delta_grad, axis=2))
    Sij = np.nan_to_num(Sij, nan=-1)
    indexMaximum = np.argmax(Sij, axis=1)
    tot_S[ActualIndices2Analyze] = Sij[np.arange(len(Sij)), indexMaximum]
    EdgeLengthMaxS[ActualIndices2Analyze] = dist[np.arange(len(dist[ActualIndices2Analyze])), indexMaximum]

    avg_S = tot_S.copy()

elif ErrorComputationMethod == "Re":

    avg_S[ActualIndices2Analyze] = np.linalg.norm(grad[ActualIndices2Analyze, :], axis=1) * np.mean(dist[ActualIndices2Analyze, :], axis=1)
    S_Std = np.std([avg_S[ActualIndices2Analyze]])

else:
    print("Error computation method", ErrorComputationMethod, "is not among the recognized options!")
    exit(2)

tot_dist = np.sum(dist, axis=1)
avg_Length = tot_dist/NNeighbors

end = time.time()

print("Done! Elapsed time ", end-start, "s")


S_thresh = np.mean([avg_S[ActualIndices2Analyze]])
print("np.mean(avg_S) = ", np.mean(avg_S[ActualIndices2Analyze]))
print("np.std(avg_S) = ", np.std(avg_S[ActualIndices2Analyze]))


print("Computing Cells volumes...")

start= time.time()
CellAreasOrVolumes = np.divide(np.sum(AllCellVolumesOrAreas[Cells], axis=1), NCells)
end = time.time()
print("Done! Elapsed time = ", end-start, " s")

# Compute current complexity
CurrentComplexity = np.sum(CellAreasOrVolumes[ActualIndices2Analyze] / (avg_Length[ActualIndices2Analyze]**options['nDim']), axis=0)

print("Current complexity of region to adapt = ", CurrentComplexity )

newSpacings = avg_Length.copy()

toWrite = False * np.ones((len(coords),))
pointsStats, notRefinedIndices, newSpacings = computeSpacing(AdaptationMethod, newSpacings, coords, ActualIndices2Analyze, avg_Length, avg_S, S_thresh, S_Std, toWrite, options, AdaptMethodOptions)

nPoints2Coarse = pointsStats[0]
nPoints2Refine = pointsStats[1]

if options['refineBySteps'] and not AdaptMethodOptions['AdaptationMethod'] == "Re" :
    print("nPoints2Coarse = ", nPoints2Coarse, " nPoints2Refine = ", nPoints2Refine)

    print("Performing step by step refinement...")
    for iStep in range(options['nStepRefinement']):
        print("Step n", iStep)
        # Recompute the threshold with the remaining points
        S_thresh = np.mean([avg_S[notRefinedIndices]])
        print("New threshold:", S_thresh)
        if AdaptMethodOptions['AdaptationMethod'] == "Re_OnlyRef":
            S_Std = np.std([avg_S[notRefinedIndices]])
            print("New std:", S_Std)

        pointsStats, notRefinedIndices, newSpacings = computeSpacing(AdaptationMethod, newSpacings, coords, notRefinedIndices, avg_Length, avg_S, S_thresh, S_Std, toWrite, options, AdaptMethodOptions, iStep)

        nPoints2Coarse += pointsStats[0]
        nPoints2Refine += pointsStats[1]

        print("nPoints2Coarse = ", pointsStats[0], " nPoints2Refine = ", pointsStats[1])

AdaptedComplexity = np.sum(CellAreasOrVolumes[ActualIndices2Analyze] / (newSpacings[ActualIndices2Analyze]**options['nDim']), axis=0)

print("Total number of points to coarse = ", nPoints2Coarse, ", Total number of points to refine = ", nPoints2Refine)

if AdaptMethodOptions['Smoothing']:

    print("Smoothing the metric distribution...")

    start = time.time()

    spacingToMean = newSpacings[Neighbors]
    if AdaptMethodOptions['UseLogExpSmoothing']:
        spacingToMean = np.log(spacingToMean)

    SmoothedSpacings = np.sum(spacingToMean, axis=1, where=whereNeighbors)/(NNeighbors+1)

    if AdaptMethodOptions['UseLogExpSmoothing']:
        SmoothedSpacings = np.exp(SmoothedSpacings)

    AdaptedComplexity = np.sum(CellAreasOrVolumes/ (SmoothedSpacings**options['nDim']), axis=0)
    newSpacings = SmoothedSpacings.copy()

    end = time.time()

    print("Done! Elapsed time ", end-start, "s")

print("Final AdaptedComplexity = ", AdaptedComplexity )


if options['UseComplexityRatio']:

    print("A complexity ratio of", options['ComplexityRatio'], "has been requested")
    print("Current complexity ratio = ", AdaptedComplexity/CurrentComplexity)
    print("Modify metric to enforce desired complexity ratio")

    tol = 1e-3
    error = tol + 1e-6
    nMax = 100
    iter = 0

    Factor = 0.0
    actualDesiredComplexity = options['ComplexityRatio'] - options['ComplexityRatio']*Factor

    while error > tol and iter < nMax:

        newSpacings[ActualIndices2Analyze] *= ((AdaptedComplexity/CurrentComplexity) / actualDesiredComplexity)**(1.0/options['nDim'])

        indices = newSpacings[ActualIndices2Analyze] < options["minSpacing"]
        newSpacings[ActualIndices2Analyze[indices]] = options["minSpacing"]
        nPointsHittingMin = max(len(np.where(indices)[0])-1, 0)
        indices = newSpacings[ActualIndices2Analyze] > options["maxSpacing"]
        newSpacings[ActualIndices2Analyze[indices]] = options["maxSpacing"]
        nPointsHittingMax = max(len(np.where(indices)[0])-1, 0)

        AdaptedComplexity = np.sum(CellAreasOrVolumes[ActualIndices2Analyze] / (newSpacings[ActualIndices2Analyze]**options['nDim']), axis=0)

        iter += 1
        error = abs(AdaptedComplexity/CurrentComplexity + Factor * AdaptedComplexity/CurrentComplexity - options['ComplexityRatio'])/options['ComplexityRatio']

    print("Reached an error =", error, "after", iter, "iteration")
    print("Final complexity ratio =", AdaptedComplexity/CurrentComplexity + Factor * AdaptedComplexity/CurrentComplexity)

else:

    # Now check if any points exceed the maximum or minimum values
    indices = newSpacings[ActualIndices2Analyze] < options["minSpacing"]
    newSpacings[ActualIndices2Analyze[indices]] = options["minSpacing"]
    newSpacings[ActualIndices2Analyze[indices]] = options["minSpacing"]
    nPointsHittingMin = max(len(np.where(indices)[0])-1, 0)
    indices = newSpacings[ActualIndices2Analyze] > options["maxSpacing"]
    newSpacings[ActualIndices2Analyze[indices]] = options["maxSpacing"]
    nPointsHittingMax = max(len(np.where(indices)[0])-1, 0)

print("N of points hitting max spacing = ", nPointsHittingMax, "  N of points hitting min spacing = ", nPointsHittingMin )


fid = open(options['fileNameOut'], 'w')
if options['WriteAllPoints']:
    toWrite[ActualIndices2Analyze] = True

for i in range(len(avg_S)):
    if toWrite[i]:
        coordsHere = coords[i]
        fid.write(str(coordsHere[0]) + " " + str(coordsHere[1]) + " " + str(coordsHere[2]) + " " + str(newSpacings[i]) + " " + str(options['decay']) + "\n")

fid.close()

endAll = time.time()

print("Total elapsed time = ", endAll-startAll, "s")
