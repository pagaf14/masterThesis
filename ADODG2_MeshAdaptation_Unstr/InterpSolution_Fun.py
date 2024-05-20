from InterpSolution_Params import options

import numpy as np
try:
    import polars as pl
    print("The Polars package is installed")
    usePolar= True
except ImportError as error:
    print("The Polars package is not installed")
    usePolar= False
try:
    from scipy.spatial import KDTree
    from scipy.interpolate import NearestNDInterpolator
    print("The scipy package is installed")
except ImportError as error:
    print(error)
    exit()

if options['MeshType'] == "CGNS":
    print("Trying to import h5py for reading CGNS mesh")
    try:
        import h5py
        print("The h5py package is installed")
    except ImportError as error:
        print(error)
        exit()

import os
import time


def readMesh_CGNS(filename, options):

    f = h5py.File(options['MeshName'],'r')
    field2Read = 'dom-1'
    if options['nDim'] == 3:
        field2Read = 'blk-1'

    # print(list(f['Base'][field2Read].keys()))
    # print(f['Base'][field2Read]['HexElements']['ElementConnectivity'][' data'][0:4])
    # print(f['Base'][field2Read]['GridCoordinates']['CoordinateY'][' data'][0])
    x = np.array(f['Base'][field2Read]['GridCoordinates']['CoordinateX'][' data'])
    y = np.array(f['Base'][field2Read]['GridCoordinates']['CoordinateY'][' data'])
    z = np.zeros((len(x),))
    if options['nDim'] == 3:
        z = np.array(f['Base'][field2Read]['GridCoordinates']['CoordinateZ'][' data'])
    Points2InterpTo = np.transpose(np.stack((x, y, z)))

    return Points2InterpTo

def readMesh_SU2(filename, options):

    fid = open(filename, 'r')
    points = []

    while True:

        # Get next line from file
        line = fid.readline()
        line = line.split('=')
        if line[0] == "NPOIN" :
            points = [[]]*int(line[1])

            for i in range(int(line[1])):
                line = fid.readline()
                line = line.split()
                x = float(line[0])
                y = float(line[1])
                z = 0.0
                if options['nDim'] == 3:
                    z = float(line[2])

                coordinates = [x, y, z]

                points[i] = coordinates

            break


    fid.close()

    return np.array(points)




def SameSide(v1, v2, v3, v4, coord, Tol):

    normal = np.cross(v2-v1, v3-v1, axis=2)
    dotV4 = np.sum(normal*(v4-v1), axis=2)
    dotP = np.sum(normal*(coord-v1), axis=2)
    # Including a tolerance in dotP == 0.0 would require to normalize all of the vectors
    # involved. This costs more than the benefits.
    return np.logical_or(np.sign(dotV4) == np.sign(dotP), dotP == 0.0);

def PointInTetra(v1, v2, v3, v4, coord, Tol):

    coord = np.repeat(coord[:, np.newaxis, :], len(v4[0]), axis=1)
    return np.logical_and(np.logical_and(SameSide(v1, v2, v3, v4, coord, Tol), SameSide(v2, v3, v4, v1, coord, Tol)), np.logical_and(SameSide(v3, v4, v1, v2, coord, Tol), SameSide(v4, v1, v2, v3, coord, Tol)))


def puntoTriangoloST (P, V0, V1, V2):

    # P = np.tile(P, (len(V0), 1))
    P = np.repeat(P[:, np.newaxis, :], len(V0[0]), axis=1)

    u = V1-V0
    v = V2-V0
    w = P-V0

    uu = np.sum(u*u, axis=2)
    vv = np.sum(v*v, axis=2)
    uv = np.sum(u*v, axis=2)
    wu = np.sum(w*u, axis=2)
    wv = np.sum(w*v, axis=2)

    den = uv*uv - uu*vv;

    s = (uv*wv-vv*wu)/den
    t = (uv*wu-uu*wv)/den

    a = np.transpose(np.stack((1-s-t, s, t)), (1, 2, 0))

    return a



def weightsForTetrahedron(p, p1, p2, p3, p4):
    # Tetrahedron vertices

    P1mP4 = p1 - p4
    P2mP4 = p2 - p4
    P3mP4 = p3 - p4
    # Matrix A
    A_matrix = np.transpose(np.stack((p1 - p4, p2 - p4, p3 - p4)), (1, 2, 0))

    # Right-hand side vector b
    b_vector = p-p4

    # Solve for the barycentric coordinates
    weights = np.linalg.solve(A_matrix, b_vector)

    # Calculate the fourth barycentric coordinate
    weights = np.hstack((weights, 1 - np.sum(weights, axis=1,keepdims=True)))


    return weights



def findPointInsideTetra(PointIDs, coords, options, PointsTo, tree):

    nElem2Check = options['nNeigh']
    # Now search among all of these tria/tetra which one actually contains the point
    nNewPoints = len(PointsTo)
    weightsAll = np.zeros((nNewPoints, options['nDim']+1))
    rightElement = np.zeros((nNewPoints, ), dtype=int)
    isInside = 0
    timeForPoints = 0.0
    skipPoints = options['skipPoints']
    nMinComputations = options['nMinComputations']
    iPointStart = 0
    timeForFun = 0.0
    isFound = np.full((nNewPoints, ), False)
    Tol = options['triaToll']
    timeForTree = 0.0


    if iPointStart+skipPoints > nNewPoints:
        skipPoints = nNewPoints

    startAll = time.time()

    while iPointStart < nNewPoints:

        Points2Compute = iPointStart+np.arange(skipPoints, dtype=int)
        Points2ComputeLocal = np.arange(skipPoints, dtype=int)

        print("Computing neighbors tria/tetra from tree structure...")
        start = time.time()
        dd, Neighbors = tree.query(PointsTo[Points2Compute], k=options['nNeigh'], workers=options['nOfCores'])

        end = time.time()
        timeForTree += end-start


        print("Computing associated triangles...")

        print("Going from ", iPointStart, " to ", iPointStart+skipPoints, " points")
        iElem = 0
        skipElements = max(options['skipElements'], int(options['nMinComputations']/skipPoints))
        if iElem+skipElements > nElem2Check:
            skipElements = nElem2Check

        startIn = time.time()

        while len(Points2Compute)>0 and iElem < nElem2Check:

            print("Going from ", iElem, " to ", iElem+skipElements, " elements")
            print("Analyzing ", len(Points2Compute), " points")

            ElemID = Neighbors[Points2ComputeLocal,iElem:(iElem+skipElements)]
            start = time.time()
            PointsCoords = coords[PointIDs[ElemID],:]
            p1 = PointsCoords[:, :, 0, :]
            p2 = PointsCoords[:, :, 1, :]
            p3 = PointsCoords[:, :, 2, :]
            p4 = PointsCoords[:, :, 3, :]
            end = time.time()
            timeForPoints+= end-start
            # Cell is a Tetra
            start = time.time()
            isInsideTetra = PointInTetra(p1, p2, p3, p4, PointsTo[Points2Compute], Tol)
            end = time.time()
            timeForFun+= end-start
            first_all_positive_index = np.argmax(isInsideTetra, axis=1)
            goodIndices = np.any(isInsideTetra, axis=1)
            notGoodIndices = np.where(np.logical_not(goodIndices))[0]
            isFound[Points2Compute] = goodIndices
            rightElement[Points2Compute[goodIndices]] = ElemID[np.where(goodIndices), first_all_positive_index[goodIndices]]

            Points2ComputeLocal = Points2ComputeLocal[notGoodIndices]
            Points2Compute = Points2Compute[notGoodIndices]


            if iElem+skipElements > nElem2Check:
                skipElements = nElem2Check-iElem
            iElem += skipElements
            if(len(Points2ComputeLocal)>0):
                skipElements = max(options['skipElements'], int(options['nMinComputations']/len(Points2ComputeLocal)))
            if iElem+skipElements > nElem2Check:
                skipElements = nElem2Check-iElem


        if iPointStart+skipPoints > nNewPoints:
            skipPoints = nNewPoints-iPointStart
        iPointStart += skipPoints
        if iPointStart+skipPoints > nNewPoints:
            skipPoints = nNewPoints-iPointStart


        endIn = time.time()
        print("Cycle time = ", endIn-startIn, "s")
        print("Points Not Found ", len(Points2ComputeLocal))


    endAll = time.time()
    print("Elapsed time for computing neighbors = ", timeForTree, "s")
    print("Elapsed time for computing associated triangles = " , endAll - startAll, " s")

    return isFound, weightsAll, rightElement, timeForPoints, timeForFun


def findPointInsideTria(PointIDs, coords, options, PointsTo, tree):

    nElem2Check = options['nNeigh']
    # Now search among all of these tria/tetra which one actually contains the point
    nNewPoints = len(PointsTo)
    weightsAll = np.zeros((nNewPoints, options['nDim']+1))
    rightElement = np.zeros((nNewPoints, ), dtype=int)
    isInside = 0
    timeForPoints = 0.0
    skipPoints = options['skipPoints']
    iPointStart = 0
    timeForFun = 0.0
    isFound = np.full((nNewPoints, ), False)
    Tol = options['triaToll']
    timeForTree = 0.0


    if iPointStart+skipPoints > nNewPoints:
        skipPoints = nNewPoints

    startAll = time.time()

    while iPointStart < nNewPoints:

        Points2Compute = iPointStart+np.arange(skipPoints, dtype=int)
        Points2ComputeLocal = np.arange(skipPoints, dtype=int)

        print("Computing neighbors tria/tetra from tree structure...")
        start = time.time()
        dd, Neighbors = tree.query(PointsTo[Points2Compute], k=options['nNeigh'], workers=options['nOfCores'])

        end = time.time()
        timeForTree += end-start


        print("Computing associated triangles...")

        print("Going from ", iPointStart, " to ", iPointStart+skipPoints, " points")
        iElem = 0
        skipElements = max(options['skipElements'], int(options['nMinComputations']/skipPoints))
        if iElem+skipElements > nElem2Check:
            skipElements = nElem2Check

        startIn = time.time()

        while len(Points2Compute)>0 and iElem < nElem2Check:

            print("Going from ", iElem, " to ", iElem+skipElements, " elements")
            print("Analyzing ", len(Points2Compute), " points")

            ElemID = Neighbors[Points2ComputeLocal,iElem:(iElem+skipElements)]
            start = time.time()
            PointsCoords = coords[PointIDs[ElemID],:]
            p1 = PointsCoords[:, :, 0, :]
            p2 = PointsCoords[:, :, 1, :]
            p3 = PointsCoords[:, :, 2, :]
            end = time.time()
            timeForPoints+= end-start
            # It is a triangle
            # Check if the point is inside
            start = time.time()
            weights = puntoTriangoloST(PointsTo[Points2Compute], p1, p2, p3)
            end = time.time()
            timeForFun+= end-start

            rows_with_all_positive = (weights > -Tol).all(axis=2)
            first_all_positive_index = np.argmax(rows_with_all_positive == True, axis=1)
            goodIndices = np.any(rows_with_all_positive, axis=1)
            notGoodIndices = np.where(np.logical_not(goodIndices))[0]
            isFound[Points2Compute] = goodIndices
            rightElement[Points2Compute[goodIndices]] = ElemID[np.where(goodIndices), first_all_positive_index[goodIndices]]
            weightsAll[Points2Compute[goodIndices]] = weights[np.where(goodIndices), first_all_positive_index[goodIndices],:]
            Points2ComputeLocal = Points2ComputeLocal[notGoodIndices]
            Points2Compute = Points2Compute[notGoodIndices]



            if iElem+skipElements > nElem2Check:
                skipElements = nElem2Check-iElem
            iElem += skipElements
            if(len(Points2ComputeLocal)>0):
                skipElements = max(options['skipElements'], int(options['nMinComputations']/len(Points2ComputeLocal)))
            if iElem+skipElements > nElem2Check:
                skipElements = nElem2Check-iElem


        if iPointStart+skipPoints > nNewPoints:
            skipPoints = nNewPoints-iPointStart
        iPointStart += skipPoints
        if iPointStart+skipPoints > nNewPoints:
            skipPoints = nNewPoints-iPointStart


        endIn = time.time()
        print("Cycle time = ", endIn-startIn, "s")
        print("Points Not Found ", len(Points2ComputeLocal))


    endAll = time.time()
    print("Elapsed time for computing neighbors = ", timeForTree, "s")
    print("Elapsed time for computing associated triangles = " , endAll - startAll, " s")

    return isFound, weightsAll, rightElement, timeForPoints, timeForFun

def LinInterp_Mine(dataFrom, PointsTo, options):

    # Create array for interpolated restart and populate with coordinates

    startEverything = time.time()

    var2skip = options['nDim']
    rows, cols = (len(dataFrom), len(PointsTo))
    dataTo = np.zeros((len(dataFrom),len(PointsTo)))
    pointsCoords = np.transpose(PointsTo)
    dataTo[0,:] = pointsCoords[0,:]
    dataTo[1,:] = pointsCoords[1,:]
    if options['nDim'] == 3:
        dataTo[2,:] = pointsCoords[2,:]

    print("Reading point coordinates...")
    start = time.time()
    coords = np.zeros((3,len(dataFrom[0])), dtype=np.float64)
    coords[0,:] = dataFrom[0,:]
    coords[1,:] = dataFrom[1,:]
    if options['nDim'] == 3:
        coords[2,:] = dataFrom[2,:]
    coords = np.transpose(coords)
    end = time.time()
    print("Done! Elapsed time", end-start,"s")
    # print(coords)

    # Now get all of the cell centroids
    print("Computing cell centroids...")
    start = time.time()

    if usePolar:
        buffer = pl.read_csv('PointIndices.csv', has_header=True, separator=',')
        buffer = pl.DataFrame.to_numpy(buffer, writable=True)
        PointIDs = buffer[:,2:]
    else:
        PointIDs = np.loadtxt('PointIndices.csv', delimiter=',', skiprows=1, dtype=int)[:,2:]

    # I am sure that all of these are tetrahedras
    CellCentroids = np.mean(coords[PointIDs,:], axis=1)
    end = time.time()
    print("Elapsed time for computing centroids = " , end - start, " s")

    tree = KDTree(CellCentroids)

    eps = 1e-10

    # Now find the tria/tetra to which each points belong. If not able to find any, then save id to be used with neares neighbor
    notFound = []
    Found = []
    # weightsAll = [[]] * len(PointsTo)
    startAll = time.time()
    timeForPoints = 0.0
    timeForFun = 0.0
    timeForWeights = 0.0

    if options['nDim'] == 2:
        isFound, weightsAll, AssociatedTriaTetra, timeForPointsHere, timeForFunHere = findPointInsideTria(PointIDs, coords, options, PointsTo, tree)
    else:
        isFound, weightsAll, AssociatedTriaTetra, timeForPointsHere, timeForFunHere = findPointInsideTetra(PointIDs, coords, options, PointsTo, tree)

    timeForPoints += timeForPointsHere
    timeForFun += timeForFunHere

    endAll = time.time()

    print("Time for points = ", timeForPoints, "s time for fun = ", timeForFun, "s")

    # end = time.time()
    # print("Elapsed time for computing associated triangles when iElem = ", iElem, " = " , end - start, " s")

    # If it is not inside any of the tria/tetra then save it for NN interpolation
    if options['nDim'] == 3:
        start = time.time()
        PointsCoords = coords[PointIDs[AssociatedTriaTetra],:]
        p1 = PointsCoords[:, 0, :]
        p2 = PointsCoords[:, 1, :]
        p3 = PointsCoords[:, 2, :]
        p4 = PointsCoords[:, 3, :]
        weightsAll = weightsForTetrahedron(PointsTo, p1, p2, p3, p4)
        end = time.time()
        timeForWeights+= end-start
    Found = np.where(isFound)[0]
    notFound = np.where(np.logical_not(isFound))[0]


    if options['nDim'] == 3:
        print("Time for weights = ", timeForWeights, "s")

    print("Number of points not found in tria/tetra = ", len(notFound))
    # XTo = [0 for i in range(len(PointsTo))]
    # YTo = [0 for i in range(len(PointsTo))]
    # ZTo = [0 for i in range(len(PointsTo))]
    # for i in range(len(PointsTo)):
    #     XTo[i] = PointsTo[i][0]
    #     YTo[i] = PointsTo[i][1]
    #     ZTo[i] = PointsTo[i][2]

    # cm = plt.cm.get_cmap('RdYlBu')
    # sc = plt.scatter(XTo,YTo, c=isFound, vmin=0, vmax=1, s=2, cmap=cm)
    # plt.colorbar(sc)
    # plt.show()

    # weightsAll = np.array(weightsAll)

    print("Interpolating solution via barycentric coordinates...")
    start = time.time()

    skipPoints = int(1e6)
    iPointStart = 0
    nNewPoints = len(Found)

    if iPointStart+skipPoints > nNewPoints:
        skipPoints = nNewPoints

    while iPointStart < nNewPoints:
        print("Going from ", iPointStart, " to ", iPointStart+skipPoints, " points")
        elemID = AssociatedTriaTetra[Found[iPointStart:(iPointStart+skipPoints)]]
        # print(elemID.shape)
        weights = weightsAll[Found[iPointStart:(iPointStart+skipPoints)]]

        data = np.transpose(dataFrom[options['nDim']:, PointIDs[elemID]], (1,0,2))
        weights = np.repeat(weights[:, np.newaxis, :], len(data[0]), axis=1)

        dataTo[options['nDim']:,Found[iPointStart:(iPointStart+skipPoints)]] = np.transpose(np.sum(data*weights, axis=2), (1,0))

        if iPointStart+skipPoints > nNewPoints:
            skipPoints = nNewPoints-iPointStart
        iPointStart += skipPoints
        if iPointStart+skipPoints > nNewPoints:
            skipPoints = nNewPoints-iPointStart

    end = time.time()
    print("Elapsed time for Interpolating solution via barycentric coordinates = " , end - start, " s")

    #-------------------------------------------------------------------------------------------------
    # Interpolate points not found

    if len(notFound) > 0:


        print("Interpolating solution via Nearest Neighbor...")
        start = time.time()

        XFrom = dataFrom[0]
        YFrom = dataFrom[1]
        ZFrom = [0 for i in range(len(dataFrom[1]))]
        if options['nDim'] == 3:
            ZFrom = dataFrom[2]


        data2Interp = dataFrom[(options['nDim']):]
        interp = NearestNDInterpolator(list(zip(XFrom, YFrom, ZFrom)), data2Interp.transpose(), rescale=True)

        InterpData = np.transpose(interp(PointsTo[notFound,0], PointsTo[notFound,1], PointsTo[notFound,2]))


        dataTo[options['nDim']:,notFound] = InterpData


        end = time.time()
        print("Elapsed time for Interpolating solution via Nearest Neighbor = " , end - start, " s")


    endEverything = time.time()
    print("Total Elapsed time = ", endEverything - startEverything, " s")
    return dataTo


def NNInterp_Mine(dataFrom, PointsTo, options):

    # Create array for interpolated restart and populate with coordinates

    startEverything = time.time()

    var2skip = options['nDim']
    rows, cols = (len(dataFrom), len(PointsTo))
    dataTo = np.zeros((len(dataFrom),len(PointsTo)))
    pointsCoords = np.transpose(PointsTo)
    dataTo[0,:] = pointsCoords[0,:]
    dataTo[1,:] = pointsCoords[1,:]
    if options['nDim'] == 3:
        dataTo[2,:] = pointsCoords[2,:]



    print("Interpolating solution via Nearest Neighbor...")
    start = time.time()

    XFrom = dataFrom[0]
    YFrom = dataFrom[1]
    ZFrom = [0 for i in range(len(dataFrom[1]))]
    if options['nDim'] == 3:
        ZFrom = dataFrom[2]


    data2Interp = dataFrom[(options['nDim']):]
    interp = NearestNDInterpolator(list(zip(XFrom, YFrom, ZFrom)), data2Interp.transpose(), rescale=True)

    InterpData = np.transpose(interp(PointsTo[:,0], PointsTo[:,1], PointsTo[:,2]))

    dataTo[options['nDim']:,:] = InterpData


    end = time.time()
    print("Elapsed time for Interpolating solution via Nearest Neighbor = " , end - start, " s")


    endEverything = time.time()
    print("Total Elapsed time = ", endEverything - startEverything, " s")
    return dataTo


def adjustFields(restartData, restartFields, options):

    for i in range(len(options['Fields2Adjust'])):

        try:
            fieldIndex = np.where(restartFields == options['Fields2Adjust'][i])[0][0]
            lessThanZeroIndices = restartData[fieldIndex, :] < 0.0
            restartData[fieldIndex, lessThanZeroIndices] = 0.0
        except:
            print("Field",options['Fields2Adjust'][i], "is not present here!")

    return restartData
