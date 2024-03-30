import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import KDTree
from scipy.interpolate import RBFInterpolator
from scipy.interpolate import NearestNDInterpolator
from scipy.interpolate import LinearNDInterpolator
import os
import shutil
import vtk
import time
import math
import multiprocessing

def readMesh(filename, nDim):

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
                if nDim == 3:
                    z = float(line[2])

                coordinates = [x, y, z]

                points[i] = coordinates

            break


    fid.close()

    return points

def crossProduct(v1, v2):

    v3 = [0]*3

    v3[0] = ((v1[1] * v2[2]) - (v1[2] * v2[1]))
    v3[1] = ((v1[2] * v2[0]) - (v1[0] * v2[2]))
    v3[2] = ((v1[0] * v2[1]) - (v1[1] * v2[0]))

    return v3

def dotProduct(v1, v2):

    return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]

def subtract(v1, v2):

    v3 = [0]*3

    v3[0] = v2[0]-v1[0]
    v3[1] = v2[1]-v1[1]
    v3[2] = v2[2]-v1[2]

    return v3

def SameSide(v1, v2, v3, v4, coord):

    isSameSide = False
    normal = crossProduct(subtract(v1, v2), subtract(v1, v3))
    dotV4 = dotProduct(normal, subtract(v1, v4))
    dotP = dotProduct(normal, subtract(v1, coord))
    return np.sign(dotV4) == np.sign(dotP);


def SameSide_Vec(v1, v2, v3, v4, coord):

    normal = np.cross(v2-v1, v3-v1)
    dotV4 = np.sum(normal*(v4-v1), axis=1)
    dotP = np.sum(normal*(coord-v1), axis=1)
    return np.sign(dotV4) == np.sign(dotP);


def PointInTetra(v1, v2, v3, v4, coord):

    return SameSide(v1, v2, v3, v4, coord) and SameSide(v2, v3, v4, v1, coord) and SameSide(v3, v4, v1, v2, coord) and SameSide(v4, v1, v2, v3, coord)


def PointInTetra_Vec(v1, v2, v3, v4, coord):

    coord = np.tile(coord, (len(v1), 1))
    return np.logical_and(np.logical_and(SameSide_Vec(v1, v2, v3, v4, coord), SameSide_Vec(v2, v3, v4, v1, coord)), np.logical_and(SameSide_Vec(v3, v4, v1, v2, coord), SameSide_Vec(v4, v1, v2, v3, coord)))

def puntoTriangoloST (P, V0, V1, V2):

    u = [0]*3
    v = [0]*3

    u[0] = V1[0]-V0[0]
    u[1] = V1[1]-V0[1]
    u[2] = V1[2]-V0[2]

    v[0] = V2[0]-V0[0]
    v[1] = V2[1]-V0[1]
    v[2] = V2[2]-V0[2]

    n = [0]*3
    n[0] = ((u[1] * v[2]) - (u[2] * v[1]))
    n[1] = ((u[2] * v[0]) - (u[0] * v[2]))
    n[2] = ((u[0] * v[1]) - (u[1] * v[0]))

    # n = np.cross(u,v)
    # normN = np.linalg.norm(n)
    normN = math.sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2])
    n[0] /= normN
    n[1] /= normN
    n[2] /= normN

    I = [0]*3
    dotHere = n[0]*(P[0]-V0[0]) + n[1]*(P[1]-V0[1]) + n[2]*(P[2]-V0[2])
    I[0] = P[0]-dotHere*n[0]
    I[1] = P[1]-dotHere*n[1]
    I[2] = P[2]-dotHere*n[2]

    # Parametric coordinates (s, t) of the intersection point I

    w = [0]*3
    w[0] = I[0]-V0[0]
    w[1] = I[1]-V0[1]
    w[2] = I[2]-V0[2]

    uu = u[0]*u[0]  +  u[1]*u[1]  +  u[2]*u[2]
    vv = v[0]*v[0]  +  v[1]*v[1]  +  v[2]*v[2]
    uv = u[0]*v[0]  +  u[1]*v[1]  +  u[2]*v[2]
    wu = w[0]*u[0]  +  w[1]*u[1]  +  w[2]*u[2]
    wv = w[0]*v[0]  +  w[1]*v[1]  +  w[2]*v[2]

    den = uv*uv -uu*vv;

    s = (uv*wv-vv*wu)/den
    t = (uv*wu-uu*wv)/den

    a = [1-s-t, s, t]


    return a

def puntoTriangoloST_Vec (P, V0, V1, V2):

    P = np.tile(P, (len(V0), 1))

    u = V1-V0
    v = V2-V0
    w = P-V0

    uu = np.sum(u*u, axis=1)
    vv = np.sum(v*v, axis=1)
    uv = np.sum(u*v, axis=1)
    wu = np.sum(w*u, axis=1)
    wv = np.sum(w*v, axis=1)

    den = uv*uv -uu*vv;

    s = (uv*wv-vv*wu)/den
    t = (uv*wu-uu*wv)/den

    a = np.transpose(np.array([1-s-t, s, t]))

    return a

def weightsForTetrahedron(p, p1, p2, p3, p4):
    # Tetrahedron vertices


    # Matrix A
    A_matrix = np.array([
        [p1[0] - p4[0], p2[0] - p4[0], p3[0] - p4[0]],
        [p1[1] - p4[1], p2[1] - p4[1], p3[1] - p4[1]],
        [p1[2] - p4[2], p2[2] - p4[2], p3[2] - p4[2]]
    ])

    # Right-hand side vector b
    b_vector = p-p4

    # Solve for the barycentric coordinates
    weights = np.linalg.solve(A_matrix, b_vector)

    # Calculate the fourth barycentric coordinate
    weights = np.append(weights, 1 - np.sum(weights))

    return weights

def computeCentroids(vtkOut):
    PointIDs = np.zeros((vtkOut.GetNumberOfCells(), nDim+1), dtype=int)
    for i in range(vtkOut.GetNumberOfCells()):

        if i%100000 == 0:
            print("Cells analyzed: ", i)


        # In theory I do not need an if. Already, if nDim is 2 then I have trias, if nDim is 3 I have tetra
        for iVertex in range(nDim+1):
            IDHere = vtkOut.GetCell(i).GetPointIds().GetId(iVertex)
            PointIDs[i][iVertex] = IDHere

    return PointIDs

def computeCentroids_parallel(vtkOut, num_processes=None):
    num_processes = num_processes or multiprocessing.cpu_count()
    with multiprocessing.Pool(processes=num_processes) as pool:
        PointIDs = pool.map(computeCentroids, vtkOut)
    return PointIDs


def findPointInsideTriaTetra(PointIDs, coords, nElem2Check, iPoint, Neighbors, p):

    isInside = 0
    iElem = 0
    timeForPoints = 0.0
    timeForif = 0.0
    timeForFun = 0.0
    weights = []
    # print("Computing associated triangles...")
    # start = time.time()
    while not isInside and iElem < nElem2Check:

        ElemID = Neighbors[iElem]
        start = time.time()
        p1 = coords[PointIDs[ElemID][0]]
        p2 = coords[PointIDs[ElemID][1]]
        p3 = coords[PointIDs[ElemID][2]]
        end = time.time()
        timeForPoints+= end-start
        if nDim == 2:
            # It is a triangle
            # Check if the point is inside
            start = time.time()
            weights = puntoTriangoloST(p, p1, p2, p3)
            end = time.time()
            timeForFun+= end-start
            start = time.time()
            if not any(item < 0.0 for item in weights):
                isInside = 1
            end = time.time()
            timeForif+= end-start

        else:
            # Cell is a Tetra
            p4 = coords[PointIDs[ElemID][3]]
            start = time.time()
            isInside = PointInTetra(p1, p2, p3, p4, p)
            end = time.time()
            timeForFun+= end-start

        iElem += 1

    return iElem, isInside, weights, ElemID, timeForPoints, timeForif, timeForFun

def findPointInsideTriaTetra_Vec(PointIDs, coords, nElem2Check, iPoint, Neighbors, p):
    # Now search among all of these tria/tetra which one actually contains the point
    weights = []
    isInside = 0
    iElem = 0
    timeForPoints = 0.0
    timeForif = 0.0
    timeForFun = 0.0
    skipElements = int(nElem2Check/2)
    rightElement = 0
    # print("i = ", iPoint, " nElem2Check = ", nElem2Check, " skipElements = ", skipElements)
    # print("Computing associated triangles...")
    # start = time.time()
    while not isInside and iElem < nElem2Check:

        ElemID = Neighbors[iElem:(iElem+skipElements)]
        start = time.time()
        p1 = coords[np.ravel(PointIDs[ElemID])[0::(nDim+1)]]
        p2 = coords[np.ravel(PointIDs[ElemID])[1::(nDim+1)]]
        p3 = coords[np.ravel(PointIDs[ElemID])[2::(nDim+1)]]
        end = time.time()
        timeForPoints+= end-start
        if nDim == 2:
            # It is a triangle
            # Check if the point is inside
            start = time.time()
            weights = puntoTriangoloST_Vec(p, p1, p2, p3)
            end = time.time()
            timeForFun+= end-start
            start = time.time()
            rows_with_all_positive = (weights > 0.0).all(axis=1)
            first_all_positive_index = np.argmax(rows_with_all_positive)

            if rows_with_all_positive.any():
                isInside = 1
                rightElement = first_all_positive_index
                weights = weights[first_all_positive_index,:]

            end = time.time()
            timeForif+= end-start

        else:
            # Cell is a Tetra
            p4 = coords[np.ravel(PointIDs[ElemID])[3::(nDim+1)]]
            start = time.time()
            isInsideTetra = PointInTetra_Vec(p1, p2, p3, p4, p)
            end = time.time()
            # print("time for fun = ", end-start)
            timeForFun+= end-start
            first_all_positive_index = np.argmax(isInsideTetra)

            if isInsideTetra.any():
                isInside = 1
                rightElement = first_all_positive_index


        if iElem+skipElements > nElem2Check:
            skipElements = nElem2Check-iElem
        iElem += skipElements

    return iElem, isInside, weights, ElemID[rightElement], timeForPoints, timeForif, timeForFun

def LinInterp_Mine(dataFrom, PointsTo, nDim, nNeigh):

    # Create array for interpolated restart and populate with coordinates

    startEverything = time.time()

    # var2skip = nDim
    # rows, cols = (len(dataFrom), len(PointsTo))
    # dataTo = [[0 for i in range(cols)] for j in range(rows)]
    # for i in range(len(PointsTo)):
    #     for iDim in range(nDim):
    #         dataTo[iDim][i] = PointsTo[i][iDim]
    var2skip = nDim
    rows, cols = (len(dataFrom), len(PointsTo))
    dataTo = np.zeros((len(dataFrom),len(PointsTo)))
    pointsCoords = np.transpose(np.array(PointsTo))
    dataTo[0,:] = pointsCoords[0,:]
    dataTo[1,:] = pointsCoords[1,:]
    if nDim == 3:
        dataTo[2,:] = pointsCoords[2,:]



    # Read triangulated/tetralized mesh
    gridreader = vtk.vtkXMLUnstructuredGridReader()
    gridreader.SetFileName("flow2InterpolateFrom_Triangulated.vtu")
    gridreader.Update()
    vtkOut = gridreader.GetOutput()
    vtkData = vtkOut.GetPoints().GetData()
    coords = np.array([vtkData.GetTuple3(x)
                          for x in range(vtkData.GetNumberOfTuples())])


    # fext = open("Mach.sol", "w")
    # for i in range(vtkOut.GetNumberOfPoints()):
    #     fext.write(str(coords[i][0]) + " " + str(coords[i][1]) + " " + str(coords[i][2]) + " " + str(dataFrom[10][i]))
    #     fext.write("\n")
    # fext.close()

    # fext = open("PointsTo.sol", "w")
    # for i in range(len(PointsTo)):
    #     fext.write(str(PointsTo[i][0]) + " " + str(PointsTo[i][1]) + " 0.0")
    #     fext.write("\n")
    # fext.close()



    # Now get all of the cell centroids
    PointIDs = np.zeros((vtkOut.GetNumberOfCells(), nDim+1), dtype=int)
    print("Computing cell centroids...")
    start = time.time()
    for i in range(vtkOut.GetNumberOfCells()):

        if i%100000 == 0:
            print("Cells analyzed: ", i)


        # In theory I do not need an if. Already, if nDim is 2 then I have trias, if nDim is 3 I have tetra
        for iVertex in range(nDim+1):
            PointIDs[i][iVertex] = vtkOut.GetCell(i).GetPointIds().GetId(iVertex)


    CellCentroids = np.sum(coords[PointIDs,:], axis=1)/(nDim+1)
    end = time.time()
    print("Elapsed time for computing centroids = " , end - start, " s")


    print("Computing maximum distance...")
    start = time.time()

    CellCentroidsHere = np.repeat(CellCentroids[:, np.newaxis, :], nDim+1, axis=1)
    distances = CellCentroidsHere - coords[PointIDs, :]
    distancesValues = np.max(np.sqrt(np.sum(np.multiply(distances, distances), axis=2)), axis=1)
    maxDistance = np.amax(distancesValues)


    end = time.time()
    print("Elapsed time for computing maximum distance = " , end - start, " s")
    # Build a tree out of these centroids
    # tree = KDTree(np.array(CellCentroids))
    tree = KDTree(CellCentroids)

    eps = 1e-10

    # Now find the tria/tetra to which each points belong. If not able to find any, then save id to be used with neares neighbor
    notFound = []
    Found = []
    isFound = [False]*len(PointsTo)
    # weightsAll = [[]] * len(PointsTo)
    weightsAll = np.zeros((len(PointsTo), nDim+1))
    AssociatedTriaTetra = [0] * len(PointsTo)
    print("Computing associated triangles...")
    startAll = time.time()
    timeForPoints = 0.0
    timeForFun = 0.0
    timeForif = 0.0
    timeForTree = 0.0
    timeForWeights = 0.0
    shouldHaveDone = len(PointsTo)*nNeigh
    actuallyDone = 0
    actuallyActuallyDone = 0

    start = time.time()
    # dd, Neighbors = tree.query(p, k=nNeigh, distance_upper_bound=maxAllowedDistance)
    # ddTotal, NeighborsTotal = tree.query(PointsTo, k=nNeigh, workers=4, distance_upper_bound=maxDistance)
    ddTotal, NeighborsTotal = tree.query(PointsTo, k=nNeigh, workers=4)
    # print(len(distancesValues))
    # print(len(distancesValues[0]))
    distancesValuesMax = np.max(distancesValues[NeighborsTotal], axis=1)
    whichAreGood = ddTotal <= np.resize(distancesValuesMax, (len(distancesValuesMax), 1))
    ddTotal = np.where(whichAreGood, ddTotal, np.nan)
    NeighborsTotal = np.where(whichAreGood, NeighborsTotal, np.nan)
    # print(NeighborsTotal[0].astype(int))
    NeighborsTotal = NeighborsTotal.astype(int)

    end = time.time()
    timeForTree += end-start

    for i in range(len(PointsTo)):

        if i%100000 == 0:
            print("Points analyzed: ", i)

        p = PointsTo[i]
        # Create a boolean mask for NaN values
        dd = ddTotal[i]
        Neighbors = NeighborsTotal[i]
        Neighbors = Neighbors[~np.isnan(dd)]
        dd = dd[~np.isnan(dd)]
        # print(Neighbors)
        # nElem2Check = len(Neighbors)
        nElem2Check = len(Neighbors)
        actuallyDone+=nElem2Check

        # Now search among all of these tria/tetra which one actually contains the point
        if nElem2Check > 20:
            iElem, isInside, weights, rightElement, timeForPointsHere, timeForIfHere, timeForFunHere = findPointInsideTriaTetra_Vec(PointIDs, coords, nElem2Check, i, Neighbors, p)
        else:
            iElem, isInside, weights, rightElement, timeForPointsHere, timeForIfHere, timeForFunHere = findPointInsideTriaTetra(PointIDs, coords, nElem2Check, i, Neighbors, p)

        timeForPoints += timeForPointsHere
        timeForif += timeForIfHere
        timeForFun += timeForFunHere

        actuallyActuallyDone += iElem

        # end = time.time()
        # print("Elapsed time for computing associated triangles when iElem = ", iElem, " = " , end - start, " s")

        # If it is not inside any of the tria/tetra then save it for NN interpolation
        if isInside == 1:
            if nDim == 3:
                start = time.time()
                p1 = coords[PointIDs[rightElement][0]]
                p2 = coords[PointIDs[rightElement][1]]
                p3 = coords[PointIDs[rightElement][2]]
                p4 = coords[PointIDs[rightElement][3]]
                weights = weightsForTetrahedron(p, p1, p2, p3, p4)
                end = time.time()
                timeForWeights+= end-start
            AssociatedTriaTetra[i] = rightElement
            # print("Point ", i ," has weights = ", weights)
            weightsAll[i,:] = np.transpose(weights)
            Found.append(i)
            isFound[i] = True
        else:
            notFound.append(i)

    endAll = time.time()
    print("Elapsed time for computing associated triangles = " , endAll - startAll, " s")
    print("Actually done = ", actuallyDone, " really done = ", actuallyActuallyDone, " should have done = ", shouldHaveDone)
    print("Time for points = ", timeForPoints, "s time for fun = ", timeForFun, "s")
    if nDim == 3:
        print("Time for weights = ", timeForWeights, "s")
    print("time for if = ", timeForif, "s time for tree = ", timeForTree, "s")

    print("Number of points not found in tria/tetra = ", len(notFound))
    XTo = [0 for i in range(len(PointsTo))]
    YTo = [0 for i in range(len(PointsTo))]
    ZTo = [0 for i in range(len(PointsTo))]
    for i in range(len(PointsTo)):
        XTo[i] = PointsTo[i][0]
        YTo[i] = PointsTo[i][1]
        ZTo[i] = PointsTo[i][2]

    # cm = plt.cm.get_cmap('RdYlBu')
    # sc = plt.scatter(XTo,YTo, c=isFound, vmin=0, vmax=1, s=2, cmap=cm)
    # plt.colorbar(sc)
    # plt.show()

    # weightsAll = np.array(weightsAll)

    print("Interpolating solution via barycentric coordinates...")
    start = time.time()
    for i in Found:

        if i%100000 == 0:
            print("Points analyzed: ", i)
        # If it has a tria/tetra associated then
        elemID = AssociatedTriaTetra[i]
        weights = weightsAll[i]

        data = dataFrom[nDim:,PointIDs[elemID][0]]
        for k in range(len(weights)-1):
            data = np.column_stack((data, dataFrom[nDim:,PointIDs[elemID][k+1]]))

        dataTo[nDim:,i] = np.dot(data, weights)

        # for j in range(len(dataTo)-nDim):
        #     index = j+nDim
        #
        #     sumHere = 0.0
        #     for k in range(len(weights)):
        #         sumHere += weights[k] * dataFrom[index][PointIDs[elemID][k]]
        #
        #     dataTo[index][i] = sumHere

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
        if nDim == 3:
            ZFrom = dataFrom[2]

        XTo = [0 for i in range(len(notFound))]
        YTo = [0 for i in range(len(notFound))]
        ZTo = [0 for i in range(len(notFound))]
        for i in range(len(notFound)):
            XTo[i] = PointsTo[notFound[i]][0]
            YTo[i] = PointsTo[notFound[i]][1]
        if nDim == 3:
            for i in range(len(notFound)):
                ZTo[i] = PointsTo[notFound[i]][2]

        data2Interp = dataFrom[(nDim):]
        interp = NearestNDInterpolator(list(zip(XFrom, YFrom, ZFrom)), data2Interp.transpose(), rescale=True)

        InterpData = interp(XTo, YTo, ZTo)

        for i in range(len(notFound)):
            if i%1000 == 0:
                print("Points analyzed: ", i)

            dataTo[nDim:,notFound[i]] = np.transpose(InterpData[i,:])


        end = time.time()
        print("Elapsed time for Interpolating solution via Nearest Neighbor = " , end - start, " s")

        # number = len(InterpData[0])
        # for j in range(len(InterpData)-1):
        #     numberHere = len(InterpData[j+1])
        #     if numberHere != number:
        #         print("PD j = ", j, " number = ", number, " numberHere = ", numberHere)


    endEverything = time.time()
    print("Total Elapsed time = ", endEverything - startEverything, " s")
    return dataTo


nDim= 2
nNeigh= 200
InterpMethod= "Barycentric"


print("Reading restart file to interpolate from... ")
start = time.time()

# First read the restart file to interp from
filenameStart="restart_2InterpFrom.dat"
A = np.fromfile(filenameStart, dtype=np.int32, count=5)
restartFields = []
for i in range(A[1]):
    strbits = np.fromfile(filenameStart, dtype=np.int8, count=33, offset=20+33*i)
    a_string = ''.join([chr(item) for item in strbits])
    # print(i, " ", a_string)
    restartFields.append(a_string)


restartData = np.fromfile(filenameStart, dtype=np.float64, count=A[1]*A[2], offset=20+33*A[1])
restartData = restartData.reshape((A[1], A[2]), order='F')
# print(restartData)
end = time.time()
print("Done! Elapsed time ", end-start, "s")


print("Reading final mesh... ")
start = time.time()
# Then read the Mesh to import to. Maybe I just need the points?
Points2InterpTo = readMesh("Mesh_2InterpTo.su2", nDim)
end = time.time()
print("Done! Elapsed time ", end-start, "s")

print("Interpolating from mesh with ", len(restartData[0]), " points into the new mesh with ", len(Points2InterpTo), " points...")

# Interpolate the restart
restartDataTo = []
if InterpMethod == "Barycentric":
    restartDataTo = LinInterp_Mine(restartData, Points2InterpTo, nDim, nNeigh)

# Finally store the new restart.
# First read the restart file to interp from
filenames = []
filename="restart_Interp_head.dat"
filenames.append(filename)
A[1] = len(restartDataTo)
A[2] = len(Points2InterpTo)
A.tofile(filename)
for i in range(A[1]):
    filename="restart_Interp_"+str(i)+".dat"
    filenames.append(filename)
    strbits = np.fromfile(filenameStart, dtype=np.int8, count=33, offset=20+33*i)
    strbits.tofile(filename)

filename="restart_Interp_Data.dat"
filenames.append(filename)
# print("AAAAAAAAAAAAAAAAAAAAAAAAAAAA")
# print(restartData[5][10000])
# print(restartDataTo[5][10000])
# typehere = type(restartDataTo[0][0])
# for j in range(len(restartDataTo)):
#     for k in range(len(restartDataTo[j])):
#         if typehere != type(restartDataTo[j][k]):
#             print("PD restartDataTo[", j, "][", k, "] = ", restartDataTo[j][k])

# restartDataTo = np.array(restartDataTo)
# print(restartDataTo[5][10000])
restartDataTo = restartDataTo.ravel(order='F')
# print(restartDataTo[4*A[2]+10000])
# print(restartDataTo[10000*A[1] + 5])
# print("AAAAAAAAAAAAAAAAAAAAAAAAAAAA")
restartDataTo.tofile(filename)


WDIR=os.getcwd()
fext=open("restart_Interp.dat","wb")
for f in filenames:
    fo=open(os.path.join(WDIR,f),"rb")
    shutil.copyfileobj(fo, fext)
    fo.close()
    os.remove(os.path.join(WDIR,f))
fext.close()



# # check file
# filename="restart_Interp.dat"
# A = np.fromfile(filename, dtype=np.int32, count=5)
# # print(A)
# for i in range(A[1]):
#     strbits = np.fromfile(filename, dtype=np.int8, count=33, offset=20+33*i)
#     a_string = ''.join([chr(item) for item in strbits])
#     restartFields.append(a_string)

# restartData = np.fromfile(filename, dtype=np.float64, count=A[1]*A[2], offset=20+33*A[1]).reshape((A[1], A[2]), order='F')
# print(restartData[5][10000])
