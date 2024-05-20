import numpy as np
import math
import time
import h5py

def getNeighborsPoints(PointIDsTot, CellTypes, nPoints):

    # Function to construct the neighbors of points structure (aka PointsOfPoint)
    # It will result in an NPoints x nNeighs array.
    # It will be a uniform array such that it can be easily casted
    # into a np.array. Thus, the row associated to point i will not have just nNeigh_i elements
    # but will extend up untile the max(nNeigh) among all of the points
    # All of the elements exceeding the i-th nNeigh will be set to the i-th point.

    # This function will also construct the structure CellsOfPoint, which is similar
    # in shape to the PointsOfPoint. It will also be a uniform array, and the remaining
    # elements of each row will be set to -1

    # Input:
    # - PointIDsTot = array of size nCells x nVertices. For each row i associated to cell i
    #                 it contains the point IDs of the vertices of the cell
    # - CellTypes   = array of size nCells x 1. Each row i associated to cell i
    #                 contains the cell type variable from VTK standards.
    # - NPoints     = integer. Number of points of the mesh

    # Output:
    # - Neighbors   = array of size nPoints x max(nNeighOfPoints). Structure explained above
    # - nNeighbors  = array of size nPoints x 1. For each row i associated to point i
    #                 it contains the real number of neighbors, aka points connected to point i
    #                 via an edge in the mesh
    # - whereNeighbors = boolean array of size nPoints x max(nNeighOfPoints).
    #                    For each row i associated to point i it has True for the first
    #                    nNeigh_i elements, and then all False.
    # - CellsOfPoint = array of size nPoints x max(nCellsOfPoint). Explained above
    # - nCells       = array of size nPoints x 1. For each row i associated to point i
    #                 it contains the real number of cells to which the point belongs

    print("Analyzing ", nPoints, " points...")
    # Start with a list of lists of nPoints x [].
    Neighbors = [[]] * nPoints
    Cells = [[]] * nPoints

    start=time.time()

    print("Arrays created.. Now we store neighbors and cells CellsOfPoint structure")


    # Cycle on each cell
    for i in range(len(CellTypes)):

        if i%10000 == 0:
            print("Cells analyzed: ", i)

        # Extract variable from array
        PointIDs = PointIDsTot[i]
        CellType = CellTypes[i]

        if CellType == 3:
            # Cell is a Line

            Point0 = PointIDs[0]
            Point1 = PointIDs[1]

            PointID = Point0
            Neighbors[PointID] = Neighbors[PointID] + [Point1]
            Cells[PointID] = Cells[PointID] + [i]

            PointID = Point1
            Neighbors[PointID] = Neighbors[PointID] + [Point0]
            Cells[PointID] = Cells[PointID] + [i]

        elif CellType == 5:
            # Cell is a triangle

            Point0 = PointIDs[0]
            Point1 = PointIDs[1]
            Point2 = PointIDs[2]

            PointID = Point0
            Neighbors[PointID] = Neighbors[PointID] + [Point1, Point2]
            Cells[PointID] = Cells[PointID] + [i]

            PointID = Point1
            Neighbors[PointID] = Neighbors[PointID] + [Point0, Point2]
            Cells[PointID] = Cells[PointID] + [i]

            PointID = Point2
            Neighbors[PointID] = Neighbors[PointID] + [Point0, Point1]
            Cells[PointID] = Cells[PointID] + [i]

        elif CellType == 9:
            # Cell is a Quad

            Point0 = PointIDs[0]
            Point1 = PointIDs[1]
            Point2 = PointIDs[2]
            Point3 = PointIDs[3]

            PointID = Point0
            Neighbors[PointID] = Neighbors[PointID] + [Point1, Point3]
            Cells[PointID] = Cells[PointID] + [i]

            PointID = Point1
            Neighbors[PointID] = Neighbors[PointID] + [Point0, Point2]
            Cells[PointID] = Cells[PointID] + [i]

            PointID = Point2
            Neighbors[PointID] = Neighbors[PointID] + [Point1, Point3]
            Cells[PointID] = Cells[PointID] + [i]

            PointID = Point3
            Neighbors[PointID] = Neighbors[PointID] + [Point0, Point2]
            Cells[PointID] = Cells[PointID] + [i]

        elif CellType == 10:
            # Cell is a Tetra

            Point0 = PointIDs[0]
            Point1 = PointIDs[1]
            Point2 = PointIDs[2]
            Point3 = PointIDs[3]

            PointID = Point0
            Neighbors[PointID] = Neighbors[PointID] + [Point1, Point2, Point3]
            Cells[PointID] = Cells[PointID] + [i]

            PointID = Point1
            Neighbors[PointID] = Neighbors[PointID] + [Point0, Point2, Point3]
            Cells[PointID] = Cells[PointID] + [i]

            PointID = Point2
            Neighbors[PointID] = Neighbors[PointID] + [Point0, Point1, Point3]
            Cells[PointID] = Cells[PointID] + [i]

            PointID = Point3
            Neighbors[PointID] = Neighbors[PointID] + [Point0, Point1, Point2]
            Cells[PointID] = Cells[PointID] + [i]

        elif CellType == 12:
            # Cell is a Hexaedron

            Point0 = PointIDs[0]
            Point1 = PointIDs[1]
            Point2 = PointIDs[2]
            Point3 = PointIDs[3]
            Point4 = PointIDs[4]
            Point5 = PointIDs[5]
            Point6 = PointIDs[6]
            Point7 = PointIDs[7]

            PointID = Point0
            Neighbors[PointID] = Neighbors[PointID] + [Point1, Point3, Point4]
            Cells[PointID] = Cells[PointID] + [i]

            PointID = Point1
            Neighbors[PointID] = Neighbors[PointID] + [Point0, Point2, Point3]
            Cells[PointID] = Cells[PointID] + [i]

            PointID = Point2
            Neighbors[PointID] = Neighbors[PointID] + [Point1, Point3, Point6]
            Cells[PointID] = Cells[PointID] + [i]

            PointID = Point3
            Neighbors[PointID] = Neighbors[PointID] + [Point0, Point2, Point7]
            Cells[PointID] = Cells[PointID] + [i]

            PointID = Point4
            Neighbors[PointID] = Neighbors[PointID] + [Point0, Point5, Point7]
            Cells[PointID] = Cells[PointID] + [i]

            PointID = Point5
            Neighbors[PointID] = Neighbors[PointID] + [Point1, Point4, Point6]
            Cells[PointID] = Cells[PointID] + [i]

            PointID = Point6
            Neighbors[PointID] = Neighbors[PointID] + [Point1, Point2, Point5]
            Cells[PointID] = Cells[PointID] + [i]

            PointID = Point7
            Neighbors[PointID] = Neighbors[PointID] + [Point3, Point4, Point6]
            Cells[PointID] = Cells[PointID] + [i]

        elif CellType == 13:
            # Cell is a Wedge (Prism with triangular base)

            Point0 = PointIDs[0]
            Point1 = PointIDs[1]
            Point2 = PointIDs[2]
            Point3 = PointIDs[3]
            Point4 = PointIDs[4]
            Point5 = PointIDs[5]

            PointID = Point0
            Neighbors[PointID] = Neighbors[PointID] + [Point1, Point2, Point3]
            Cells[PointID] = Cells[PointID] + [i]

            PointID = Point1
            Neighbors[PointID] = Neighbors[PointID] + [Point0, Point2, Point4]
            Cells[PointID] = Cells[PointID] + [i]

            PointID = Point2
            Neighbors[PointID] = Neighbors[PointID] + [Point0, Point1, Point5]
            Cells[PointID] = Cells[PointID] + [i]

            PointID = Point3
            Neighbors[PointID] = Neighbors[PointID] + [Point0, Point4, Point5]
            Cells[PointID] = Cells[PointID] + [i]

            PointID = Point4
            Neighbors[PointID] = Neighbors[PointID] + [Point1, Point3, Point5]
            Cells[PointID] = Cells[PointID] + [i]

            PointID = Point5
            Neighbors[PointID] = Neighbors[PointID] + [Point2, Point3, Point4]
            Cells[PointID] = Cells[PointID] + [i]

        elif CellType == 14:
            # Cell is a Pyramid

            Point0 = PointIDs[0]
            Point1 = PointIDs[1]
            Point2 = PointIDs[2]
            Point3 = PointIDs[3]
            Point4 = PointIDs[4]

            PointID = Point0
            Neighbors[PointID] = Neighbors[PointID] + [Point1, Point3, Point4]
            Cells[PointID] = Cells[PointID] + [i]

            PointID = Point1
            Neighbors[PointID] = Neighbors[PointID] + [Point0, Point2, Point4]
            Cells[PointID] = Cells[PointID] + [i]

            PointID = Point2
            Neighbors[PointID] = Neighbors[PointID] + [Point1, Point3, Point4]
            Cells[PointID] = Cells[PointID] + [i]

            PointID = Point3
            Neighbors[PointID] = Neighbors[PointID] + [Point0, Point2, Point4]
            Cells[PointID] = Cells[PointID] + [i]

            PointID = Point4
            Neighbors[PointID] = Neighbors[PointID] + [Point0, Point1, Point2, Point3]

            Cells[PointID] = Cells[PointID] + [i]


        else:
            print(CellType, " is not a known type")

    end = time.time()
    print("Elapsed time ", end-start, "s")

    # Now I need to remove all of the duplicates
    NNeighbors = np.zeros((nPoints,),dtype=int)
    NCells = np.zeros((nPoints,),dtype=int)
    for nPoint in range(nPoints):
        Neighbors[nPoint] = list(dict.fromkeys(Neighbors[nPoint]))
        NNeighbors[nPoint] = len(Neighbors[nPoint])
        Cells[nPoint] = list(dict.fromkeys(Cells[nPoint]))
        NCells[nPoint] = len(Cells[nPoint])

    # And now I construct the uniform list of lists to be casted then into an np.array
    maxNNeigh = np.max(NNeighbors)
    maxNCell = np.max(NCells)
    whereNeighbors = np.full((nPoints,maxNNeigh+1), True)
    for nPoint in range(nPoints):
        Neighbors[nPoint] = Neighbors[nPoint]+ [nPoint]*(maxNNeigh-NNeighbors[nPoint]+1)
        whereNeighbors[nPoint,(NNeighbors[nPoint]+1):] = False
        Cells[nPoint] = Cells[nPoint]+ [-1]*(maxNCell-NCells[nPoint])

    return np.array(Neighbors), NNeighbors, whereNeighbors, np.array(Cells), NCells



def computeCellVolumesOrAreas(PointIDs, CellTypes, coords, options):

    # Function that computes the area (or volume) of the cell.
    # Inputs:
    # - PointIDs = array of size nCells x nVertices. For each row i associated to cell i
    #              it contains the point IDs of the vertices of the cell
    # - CellTypes= array of size nCells x 1. Each row i associated to cell i
    #                 contains the cell type variable from VTK standards.
    # - coords   = array of size nPoints x 3. It contains the coordinates of the
    #              i-th point.
    # - options  = dictionary of user-input options

    # Output:
    # - CellAreasOrVolumes = array of size nCells+1. It contains the area (or volume)
    #                        of each cell. The last element of the array has value 0
    #                        and is used due to the uniformity of the CellsOfPoint array

    # Initialize the vector
    CellAreasOrVolumes = np.zeros((len(CellTypes)+1,))

    print("Computing Cells volumes...")

    start = time.time()
    if options['nDim'] == 1:
        lineElems = np.where(CellTypes==3)[0]
        print("N of lines: ", len(lineElems))
        IDsLines = PointIDs[lineElems, :]

        CellAreasOrVolumes[lineElems] = np.linalg.norm(coords[IDsLines[:, 0],:]-coords[IDsLines[:, 1],:])

    elif options['nDim'] == 2:
        triaElems = np.where(CellTypes==5)[0]
        quadElems = np.where(CellTypes==9)[0]

        print("N of trias: ", len(triaElems))
        print("N of quads: ", len(quadElems))

        # There should always be trias
        IDsTria = PointIDs[triaElems, :]
        CellAreasOrVolumes[triaElems] = computeTriangleArea_Vec(IDsTria[:, 0], IDsTria[:, 1], IDsTria[:, 2], coords)
        if len(quadElems)>0:
            IDsQuads = PointIDs[quadElems, :]
            CellAreasOrVolumes[quadElems] = computeQuadArea_Vec(IDsQuads, coords)

    else:
        tetraElems = np.where(CellTypes==10)[0]
        hexaElems = np.where(CellTypes==12)[0]
        wedgeElems = np.where(CellTypes==13)[0]
        pyraElems = np.where(CellTypes==14)[0]

        print("N of tetras: ", len(tetraElems))
        print("N of hexas: ", len(hexaElems))
        print("N of wedges: ", len(wedgeElems))
        print("N of pyras: ", len(pyraElems))

        # I need this 'cause the formula for the volume of a tetrahedron requires it
        # Now the i-th row of the coords vector is x_i,y_i,z_i,1
        coordsHere = np.append(coords, np.ones((len(coords),1)), axis=1)

        # There should always be tetras
        IDsTetra = PointIDs[tetraElems, :]
        CellAreasOrVolumes[tetraElems] = computeTetraVolume_Vec(IDsTetra[:, 0], IDsTetra[:, 1], IDsTetra[:, 2], IDsTetra[:, 3], coordsHere)

        if len(hexaElems)>0:
            IDsHexa = PointIDs[hexaElems]
            CellAreasOrVolumes[hexaElems] = computeHexaVolume_Vec(IDsHexa, coords)

        if len(wedgeElems)>0:
            IDsWedge = PointIDs[wedgeElems]
            CellAreasOrVolumes[wedgeElems] = computeWedgeVolume_Vec(IDsWedge, coordsHere)

        if len(pyraElems)>0:
            IDsPyra = PointIDs[pyraElems]
            CellAreasOrVolumes[pyraElems] = computePyramidVolume_Vec(IDsPyra, coordsHere)

    end = time.time()
    print("Done! Elapsed time ", end-start, "s")

    # Is it correct?
    CellAreasOrVolumes[len(CellTypes)] = 0.0

    return CellAreasOrVolumes

def computeTriangleArea(p1, p2, p3):

    return 0.5*np.linalg.norm(np.cross(p2-p1, p3-p1))

def computeTriangleArea_Vec(p1IDs, p2IDs, p3IDs, coords):

    p1 = coords[p1IDs]
    p2 = coords[p2IDs]
    p3 = coords[p3IDs]

    return 0.5*np.linalg.norm(np.cross(p2-p1, p3-p1, axis=1), axis=1)


def computeQuadArea(quad, coords):


    Area = computeTriangleArea(coords[quad[0]], coords[quad[1]], coords[quad[2]])
    Area += computeTriangleArea(coords[quad[0]], coords[quad[2]], coords[quad[3]])

    return Area

def computeQuadArea_Vec(quad, coords):


    Area = computeTriangleArea_Vec(quad[:,0], quad[:,1], quad[:,2], coords)
    Area += computeTriangleArea_Vec(quad[:,0], quad[:,2], quad[:,3], coords)

    return Area

def computeTetraVolume(p1, p2, p3, p4):

    p1 = np.append(p1, 1)
    p2 = np.append(p2, 1)
    p3 = np.append(p3, 1)
    p4 = np.append(p4, 1)

    A = np.stack((p1, p2, p3, p4))

    return (1/6.0)*np.abs(np.linalg.det(A))

def computeTetraVolume_Vec(p1IDs, p2IDs, p3IDs, p4IDs, coords):

    p1 = coords[p1IDs]
    p2 = coords[p2IDs]
    p3 = coords[p3IDs]
    p4 = coords[p4IDs]

    A = np.transpose(np.stack((p1, p2, p3, p4)), (1, 0, 2))

    return (1/6.0)*np.abs(np.linalg.det(A))

def computeHexaVolume(hexa, coords):

    Volume = 0.0


    # First approximation: It is just an extrusion at constant base area

    # Create the first wedge

    p1 = coords[hexa[0]]
    p2 = coords[hexa[1]]
    p3 = coords[hexa[4]]
    p4 = coords[hexa[5]]
    p6 = coords[hexa[6]]

    a = np.linalg.norm(p2 - p1)
    b = np.linalg.norm(p3 - p2)
    c = np.linalg.norm(p3 - p1)

    s = (a+b+c)/2

    BaseArea = math.sqrt(s * (s-a) * (s-b) * (s-c))

    a = np.linalg.norm(p4 - p3)
    b = np.linalg.norm(p1 - p4)
    c = np.linalg.norm(p1 - p3)

    s = (a+b+c)/2

    BaseArea += math.sqrt(s * (s-a) * (s-b) * (s-c))

    Height = np.linalg.norm(p6 - p4)

    return BaseArea*Height


def computeHexaVolume_Vec(hexa, coords):

    Volume = 0.0


    # First approximation: It is just an extrusion at constant base area

    # Create the first wedge

    p1 = coords[hexa[:,0]]
    p2 = coords[hexa[:,1]]
    p3 = coords[hexa[:,4]]
    p4 = coords[hexa[:,5]]
    p6 = coords[hexa[:,6]]

    a = np.linalg.norm(p2 - p1, axis=1)
    b = np.linalg.norm(p3 - p2, axis=1)
    c = np.linalg.norm(p3 - p1, axis=1)

    s = (a+b+c)/2

    BaseArea = np.sqrt(s * (s-a) * (s-b) * (s-c))

    a = np.linalg.norm(p4 - p3, axis=1)
    b = np.linalg.norm(p1 - p4, axis=1)
    c = np.linalg.norm(p1 - p3, axis=1)

    s = (a+b+c)/2

    BaseArea += np.sqrt(s * (s-a) * (s-b) * (s-c))

    Height = np.linalg.norm(p6 - p4, axis=1)

    return BaseArea*Height

def computeWedgeVolume(wedge, coords):

    Volume = 0.0


    p1 = coords[wedge[0]]
    p2 = coords[wedge[1]]
    p3 = coords[wedge[2]]
    p4 = coords[wedge[3]]

    Volume += computeTetraVolume(p1, p2, p3, p4)

    p1 = coords[wedge[2]]
    p2 = coords[wedge[3]]
    p3 = coords[wedge[4]]
    p4 = coords[wedge[5]]

    Volume += computeTetraVolume(p1, p2, p3, p4)

    p1 = coords[wedge[1]]
    p2 = coords[wedge[2]]
    p3 = coords[wedge[3]]
    p4 = coords[wedge[4]]

    Volume += computeTetraVolume(p1, p2, p3, p4)

    return Volume

def computeWedgeVolume_Vec(wedge, coords):

    Volume = 0.0


    p1 = wedge[:,0]
    p2 = wedge[:,1]
    p3 = wedge[:,2]
    p4 = wedge[:,3]

    Volume += computeTetraVolume_Vec(p1, p2, p3, p4, coords)

    p1 = wedge[:,2]
    p2 = wedge[:,3]
    p3 = wedge[:,4]
    p4 = wedge[:,5]

    Volume += computeTetraVolume_Vec(p1, p2, p3, p4, coords)

    p1 = wedge[:,1]
    p2 = wedge[:,2]
    p3 = wedge[:,3]
    p4 = wedge[:,4]

    Volume += computeTetraVolume_Vec(p1, p2, p3, p4, coords)

    return Volume

def computePyramidVolume(pyra, coords):

    Volume = 0.0

    p1 = coords[pyra[0]]
    p2 = coords[pyra[1]]
    p3 = coords[pyra[3]]
    p4 = coords[pyra[4]]

    Volume += computeTetraVolume(p1, p2, p3, p4)

    p1 = coords[pyra[1]]
    p2 = coords[pyra[2]]
    p3 = coords[pyra[3]]
    p4 = coords[pyra[4]]

    Volume += computeTetraVolume(p1, p2, p3, p4)

    return Volume

def computePyramidVolume_Vec(pyra, coords):

    Volume = 0.0

    p1 = pyra[:,0]
    p2 = pyra[:,1]
    p3 = pyra[:,3]
    p4 = pyra[:,4]

    Volume += computeTetraVolume_Vec(p1, p2, p3, p4, coords)

    p1 = pyra[:,1]
    p2 = pyra[:,2]
    p3 = pyra[:,3]
    p4 = pyra[:,4]

    Volume += computeTetraVolume_Vec(p1, p2, p3, p4, coords)

    return Volume

def computeSpacing(AdaptationMethod, newSpacings, coords, ActualIndices2Analyze, avg_Length, avg_S, S_thresh, S_Std, toWrite, options, AdaptMethodOptions, iStep=0):

    ## All of the external lists are callable here
    nPoints2Coarse = 0
    nPoints2Refine = 0
    nPointsHittingMax = 0
    nPointsHittingMin = 0

    notRefinedIndices = []

    if AdaptationMethod == "Pointwise_OnlyRef":

        coordsHere = coords[ActualIndices2Analyze]
        spacing = avg_Length[ActualIndices2Analyze]
        avg_S_Here = avg_S[ActualIndices2Analyze]
        indices = avg_S_Here > S_thresh
        indices2Mod = ActualIndices2Analyze[indices]
        toWrite[indices2Mod] = True
        spacing[indices] = avg_Length[indices2Mod]*((S_thresh*np.ones((len(indices2Mod),))/avg_S[indices2Mod])**(1.0/AdaptMethodOptions["p"]))
        # indices = spacing < options["minSpacing"]
        # spacing[indices] = options["minSpacing"]
        # nPointsHittingMin = max(len(np.where(indices)[0])-1, 0)
        nPoints2Refine = len(indices2Mod)
        newSpacings[ActualIndices2Analyze] = spacing[:]

        indices = avg_S_Here <= S_thresh
        notRefinedIndices = ActualIndices2Analyze[indices]

    elif AdaptationMethod == "Pointwise_WithCoarsening":

        coordsHere = coords[ActualIndices2Analyze]
        spacing = avg_Length[ActualIndices2Analyze]
        avg_S_Here = avg_S[ActualIndices2Analyze]
        indices = avg_S_Here > S_thresh
        indices2Mod = ActualIndices2Analyze[indices]
        toWrite[indices2Mod] = True
        spacing[indices] = avg_Length[indices2Mod]*((S_thresh*np.ones((len(indices2Mod),))/avg_S[indices2Mod])**(1.0/AdaptMethodOptions["p"]))
        # indices = spacing < options["minSpacing"]
        # spacing[indices] = options["minSpacing"]
        # nPointsHittingMin = max(len(np.where(indices)[0])-1, 0)
        nPoints2Refine = len(indices2Mod)
        newSpacings[ActualIndices2Analyze] = spacing[:]

        indices = avg_S_Here <= S_thresh
        notRefinedIndices = ActualIndices2Analyze[indices]

        doIt = True
        # I have to coarsen only at the last step of the refinement
        # if refinement by steps is on
        if options['refineBySteps'] and iStep < (options['nStepRefinement']-1):
            doIt = False
            print("Not performing coarsening now")

        if doIt:
            indices = avg_S_Here < (S_thresh/AdaptMethodOptions["CoarseningFactor"])
            indices2Mod = ActualIndices2Analyze[indices]
            toWrite[indices2Mod] = True
            coarsening = (((S_thresh/AdaptMethodOptions["CoarseningFactor"])*np.ones((len(indices2Mod),))/avg_S[indices2Mod])**(1.0/AdaptMethodOptions["p"]))
            indicesOver = coarsening > AdaptMethodOptions["Eps_Coarse"]
            print("N of points requiring more than",AdaptMethodOptions["Eps_Coarse"], "=",len(np.where(indicesOver)[0])-1)
            coarsening[indicesOver] = AdaptMethodOptions["Eps_Coarse"]
            spacing[indices] = avg_Length[indices2Mod]*coarsening
            # indices = spacing > options["maxSpacing"]
            # spacing[indices] = options["maxSpacing"]
            # nPointsHittingMax = max(len(np.where(indices)[0])-1, 0)
            nPoints2Coarse = len(indices2Mod)
            newSpacings[ActualIndices2Analyze] = spacing[:]

    elif AdaptationMethod == "Pointwise_Clusters":

        # THIS ONE HAS TO BE FIXED WITH THE NEW VECTORS, SEE Pointwise_OnlyRef and Pointwise_WithCoarsening

        # Sort the clusters by value of the threshold
        sort_index = np.argsort(centroid)
        # print(sort_index)
        for i in range(len(avg_S)):
            coordsHere = coords[i]
            spacing = avg_Length[i]

            if isInsideTheBox[i]:
                toWrite[i] = True
                currentThreshold = centroid[label[i]]

                # Check if it is lower than the threshold set by the centroid
                if avg_S[i] < currentThreshold:
                    # If it's the one with lowest threshold, then I should just coarsen it
                    if currentThreshold == min(centroid):
                        spacing = avg_Length[i]*Eps_Coarse
                        if spacing > options["maxSpacing"]:
                            spacing = options["maxSpacing"]
                            nPointsHittingMax += 1
                        nPoints2Coarse+=1

                        # XOfCoarsenedNodes.append(coords[i][0])
                        # YOfCoarsenedNodes.append(coords[i][1])

                    # If it's whatever else threshold, than I need to refine it according to the previous one in order of values
                    else:
                        # print("This point has label ", label[i])
                        # print("Its avg_S is ", avg_S[i], " and the corresponding centroid is ", currentThreshold)
                        currentThresholdIndex = np.where(sort_index==label[i])[0][0]
                        previousThresholdIndex = sort_index[currentThresholdIndex-1]
                        # print("Its position in order is n ", currentThresholdIndex, " and the previous centroid is the number ", previousThresholdIndex)
                        spacing = EdgeLengthMaxS[i]*(centroid[previousThresholdIndex]/avg_S[i])**(1.0/AdaptMethodOptions["p"])
                        if spacing < options["minSpacing"]:
                            spacing = options["minSpacing"]
                            nPointsHittingMin += 1
                        nPoints2Refine+=1

                        # XOfRefinedNodes.append(coords[i][0])
                        # YOfRefinedNodes.append(coords[i][1])

                else:
                    # If it's greater than current threshold, then I will refine with this threshold in mind
                    spacing = EdgeLengthMaxS[i]*(currentThreshold/avg_S[i])**(1.0/AdaptMethodOptions["p"])
                    if spacing < options["minSpacing"]:
                        spacing = options["minSpacing"]
                        nPointsHittingMin += 1
                    nPoints2Refine+=1

                    # XOfRefinedNodes.append(coords[i][0])
                    # YOfRefinedNodes.append(coords[i][1])

            else:
                toWrite[i] = False

            newSpacings[i] = spacing

    elif AdaptationMethod == "Re":

        spacing = avg_Length[ActualIndices2Analyze]
        avg_S_Here = avg_S[ActualIndices2Analyze]

        # Refine the points

        S_thresh_1 = S_thresh
        nPoints2Refine = 0
        for i in range(AdaptMethodOptions['NRef']):
            S_thresh_2 = S_thresh + AdaptMethodOptions['RefCoeffs'][i] * S_Std
            print("S_thresh_1 =", S_thresh_1, "S_thresh_2 =", S_thresh_2)
            indices = np.logical_and(avg_S_Here >= S_thresh_1, avg_S_Here < S_thresh_2)
            indices2Mod = ActualIndices2Analyze[indices]
            nPoints2Refine_Here = len(indices2Mod)
            if nPoints2Refine_Here > 0:
                nPoints2Refine += nPoints2Refine_Here
                toWrite[indices2Mod] = True
                spacing[indices] = avg_Length[indices2Mod]/AdaptMethodOptions['RefScaling'][i]

            S_thresh_1 = S_thresh_2

        indices = avg_S_Here >= S_thresh_1
        indices2Mod = ActualIndices2Analyze[indices]
        nPoints2Refine_Here = len(indices2Mod)
        if nPoints2Refine_Here > 0:
            nPoints2Refine += nPoints2Refine_Here
            toWrite[indices2Mod] = True
            spacing[indices] = avg_Length[indices2Mod]/AdaptMethodOptions['RefScaling'][-1]


        # Coarsen the points
        S_thresh_1 = S_thresh
        nPoints2Coarse = 0
        for i in range(AdaptMethodOptions['NCoarse']):
            S_thresh_2 = S_thresh - AdaptMethodOptions['CoarseCoeffs'][i] * S_Std
            print("S_thresh_1 =", S_thresh_1, "S_thresh_2 =", S_thresh_2)
            indices = np.logical_and(avg_S_Here <= S_thresh_1, avg_S_Here > S_thresh_2)
            indices2Mod = ActualIndices2Analyze[indices]
            nPoints2Coarse_Here = len(indices2Mod)
            if nPoints2Coarse_Here > 0:
                nPoints2Coarse += nPoints2Coarse_Here
                toWrite[indices2Mod] = True
                spacing[indices] = avg_Length[indices2Mod]*AdaptMethodOptions['CoarseScaling'][i]

            S_thresh_1 = S_thresh_2

        indices = avg_S_Here <= S_thresh_1
        indices2Mod = ActualIndices2Analyze[indices]
        nPoints2Coarse_Here = len(indices2Mod)
        if nPoints2Coarse_Here > 0:
            nPoints2Coarse += nPoints2Coarse_Here
            toWrite[indices2Mod] = True
            spacing[indices] = avg_Length[indices2Mod]*AdaptMethodOptions['CoarseScaling'][-1]


        newSpacings[ActualIndices2Analyze] = spacing[:]

    elif AdaptationMethod == "Re_OnlyRef":

        spacing = avg_Length[ActualIndices2Analyze]
        avg_S_Here = avg_S[ActualIndices2Analyze]

        indices = avg_S_Here <= S_thresh
        notRefinedIndices = ActualIndices2Analyze[indices]

        # Refine the points

        S_thresh_1 = S_thresh
        nPoints2Refine = 0
        for i in range(AdaptMethodOptions['NRef']):
            S_thresh_2 = S_thresh + AdaptMethodOptions['RefCoeffs'][i] * S_Std
            print("S_thresh_1 =", S_thresh_1, "S_thresh_2 =", S_thresh_2)
            indices = np.logical_and(avg_S_Here >= S_thresh_1, avg_S_Here < S_thresh_2)
            indices2Mod = ActualIndices2Analyze[indices]
            nPoints2Refine_Here = len(indices2Mod)
            if nPoints2Refine_Here > 0:
                nPoints2Refine += nPoints2Refine_Here
                toWrite[indices2Mod] = True
                spacing[indices] = avg_Length[indices2Mod]/AdaptMethodOptions['RefScaling'][i]

            S_thresh_1 = S_thresh_2

        indices = avg_S_Here >= S_thresh_1
        indices2Mod = ActualIndices2Analyze[indices]
        nPoints2Refine_Here = len(indices2Mod)
        if nPoints2Refine_Here > 0:
            nPoints2Refine += nPoints2Refine_Here
            toWrite[indices2Mod] = True
            spacing[indices] = avg_Length[indices2Mod]/AdaptMethodOptions['RefScaling'][-1]


        newSpacings[ActualIndices2Analyze] = spacing[:]




    pointsStats = [nPoints2Coarse, nPoints2Refine]
    return pointsStats, notRefinedIndices, newSpacings
