import numpy as np
import vtk
from matplotlib import pyplot as plt
from scipy.cluster.vq import kmeans2
import math
from vtkmodules.vtkCommonDataModel import vtkTetra
from vtkmodules.vtkCommonDataModel import vtkWedge
import time

def getNeighborsPoints_New(vtkOut):

    print("Analyzing ", str(vtkOut.GetNumberOfPoints()), " points...")
    rows, cols = (vtkOut.GetNumberOfPoints(), 200)
    Neighbors = [[]] * vtkOut.GetNumberOfPoints()
    Cells = [[]] * vtkOut.GetNumberOfPoints()

    # rows= vtkOut.GetNumberOfPoints()
    # skipIndex = np.zeros((rows,), dtype=int)
    # skipIndexCells = np.copy(skipIndex)

    print("Arrays created.. Now we store neighbors and cells CellsOfPoint structure")

    for i in range(vtkOut.GetNumberOfCells()):

        if i%10000 == 0:
            print("Cells analyzed: ", i)

        if vtkOut.GetCellType(i) == 5:
            # Cell is a triangle

            # print()
            # print(vtkOut.GetCell(i).GetPointIds().GetId(0), " ", vtkOut.GetCell(i).GetPointIds().GetId(1), " ", vtkOut.GetCell(i).GetPointIds().GetId(2))
            # print(vtkOut.GetCell(i).GetPointIds().GetId(0), " ", Neighbors[vtkOut.GetCell(i).GetPointIds().GetId(0)])
            PointID = vtkOut.GetCell(i).GetPointIds().GetId(0)
            Neighbors[PointID] = Neighbors[PointID] + [vtkOut.GetCell(i).GetPointIds().GetId(1),
                                       vtkOut.GetCell(i).GetPointIds().GetId(2)]
            Cells[PointID] = Cells[PointID] + [i]

            # print(PointID, " ", " ", skipIndex[PointID])
            PointID = vtkOut.GetCell(i).GetPointIds().GetId(1)
            Neighbors[PointID] = Neighbors[PointID] + [vtkOut.GetCell(i).GetPointIds().GetId(0),
                                                                vtkOut.GetCell(i).GetPointIds().GetId(2)]
            Cells[PointID] = Cells[PointID] + [i]

            PointID = vtkOut.GetCell(i).GetPointIds().GetId(2)
            Neighbors[PointID] = Neighbors[PointID] + [vtkOut.GetCell(i).GetPointIds().GetId(0),
                                                                vtkOut.GetCell(i).GetPointIds().GetId(1)]
            Cells[PointID] = Cells[PointID] + [i]
            # print()

        elif vtkOut.GetCellType(i) == 9:
            # Cell is a Quad

            # print("It is a quadrilater")
            PointID = vtkOut.GetCell(i).GetPointIds().GetId(0)
            Neighbors[PointID] = Neighbors[PointID] + [vtkOut.GetCell(i).GetPointIds().GetId(1),
                                                                vtkOut.GetCell(i).GetPointIds().GetId(3)]
            Cells[PointID] = Cells[PointID] + [i]

            PointID = vtkOut.GetCell(i).GetPointIds().GetId(1)
            Neighbors[PointID] = Neighbors[PointID] + [vtkOut.GetCell(i).GetPointIds().GetId(0),
                                                                vtkOut.GetCell(i).GetPointIds().GetId(2)]
            Cells[PointID] = Cells[PointID] + [i]

            PointID = vtkOut.GetCell(i).GetPointIds().GetId(2)
            Neighbors[PointID] = Neighbors[PointID] + [vtkOut.GetCell(i).GetPointIds().GetId(1),
                                                                vtkOut.GetCell(i).GetPointIds().GetId(3)]
            Cells[PointID] = Cells[PointID] + [i]

            PointID = vtkOut.GetCell(i).GetPointIds().GetId(3)
            Neighbors[PointID] = Neighbors[PointID] + [vtkOut.GetCell(i).GetPointIds().GetId(0),
                                                                vtkOut.GetCell(i).GetPointIds().GetId(2)]
            Cells[PointID] = Cells[PointID] + [i]

        elif vtkOut.GetCellType(i) == 10:
            # Cell is a Tetra

            PointID = vtkOut.GetCell(i).GetPointIds().GetId(0)
            Neighbors[PointID] = Neighbors[PointID] + [vtkOut.GetCell(i).GetPointIds().GetId(1),
                                                       vtkOut.GetCell(i).GetPointIds().GetId(2),
                                                       vtkOut.GetCell(i).GetPointIds().GetId(3)]
            Cells[PointID] = Cells[PointID] + [i]

            PointID = vtkOut.GetCell(i).GetPointIds().GetId(1)
            Neighbors[PointID] = Neighbors[PointID] + [vtkOut.GetCell(i).GetPointIds().GetId(0),
                                       vtkOut.GetCell(i).GetPointIds().GetId(2),
                                       vtkOut.GetCell(i).GetPointIds().GetId(3)]
            Cells[PointID] = Cells[PointID] + [i]

            PointID = vtkOut.GetCell(i).GetPointIds().GetId(2)
            Neighbors[PointID] = Neighbors[PointID] + [vtkOut.GetCell(i).GetPointIds().GetId(0),
                                       vtkOut.GetCell(i).GetPointIds().GetId(1),
                                       vtkOut.GetCell(i).GetPointIds().GetId(3)]
            Cells[PointID] = Cells[PointID] + [i]

            PointID = vtkOut.GetCell(i).GetPointIds().GetId(3)
            Neighbors[PointID] = Neighbors[PointID] + [vtkOut.GetCell(i).GetPointIds().GetId(0),
                                                                vtkOut.GetCell(i).GetPointIds().GetId(1),
                                                                vtkOut.GetCell(i).GetPointIds().GetId(2)]
            Cells[PointID] = Cells[PointID] + [i]

        elif vtkOut.GetCellType(i) == 12:
            # Cell is a Hexaedron

            PointID = vtkOut.GetCell(i).GetPointIds().GetId(0)
            Neighbors[PointID] = Neighbors[PointID] + [vtkOut.GetCell(i).GetPointIds().GetId(1),
                                                                vtkOut.GetCell(i).GetPointIds().GetId(3),
                                                                vtkOut.GetCell(i).GetPointIds().GetId(4)]
            Cells[PointID] = Cells[PointID] + [i]

            PointID = vtkOut.GetCell(i).GetPointIds().GetId(1)
            Neighbors[PointID] = Neighbors[PointID] + [vtkOut.GetCell(i).GetPointIds().GetId(0),
                                                                vtkOut.GetCell(i).GetPointIds().GetId(2),
                                                                vtkOut.GetCell(i).GetPointIds().GetId(3)]
            Cells[PointID] = Cells[PointID] + [i]

            PointID = vtkOut.GetCell(i).GetPointIds().GetId(2)
            Neighbors[PointID] = Neighbors[PointID] + [vtkOut.GetCell(i).GetPointIds().GetId(1),
                                                                vtkOut.GetCell(i).GetPointIds().GetId(3),
                                                                vtkOut.GetCell(i).GetPointIds().GetId(6)]
            Cells[PointID] = Cells[PointID] + [i]

            PointID = vtkOut.GetCell(i).GetPointIds().GetId(3)
            Neighbors[PointID] = Neighbors[PointID] + [vtkOut.GetCell(i).GetPointIds().GetId(0),
                                                                vtkOut.GetCell(i).GetPointIds().GetId(2),
                                                                vtkOut.GetCell(i).GetPointIds().GetId(7)]
            Cells[PointID] = Cells[PointID] + [i]

            PointID = vtkOut.GetCell(i).GetPointIds().GetId(4)
            Neighbors[PointID] = Neighbors[PointID] + [vtkOut.GetCell(i).GetPointIds().GetId(0),
                                                                vtkOut.GetCell(i).GetPointIds().GetId(5),
                                                                vtkOut.GetCell(i).GetPointIds().GetId(7)]
            Cells[PointID] = Cells[PointID] + [i]

            PointID = vtkOut.GetCell(i).GetPointIds().GetId(5)
            Neighbors[PointID] = Neighbors[PointID] + [vtkOut.GetCell(i).GetPointIds().GetId(1),
                                                                vtkOut.GetCell(i).GetPointIds().GetId(4),
                                                                vtkOut.GetCell(i).GetPointIds().GetId(6)]
            Cells[PointID] = Cells[PointID] + [i]

            PointID = vtkOut.GetCell(i).GetPointIds().GetId(6)
            Neighbors[PointID] = Neighbors[PointID] + [vtkOut.GetCell(i).GetPointIds().GetId(1),
                                                                vtkOut.GetCell(i).GetPointIds().GetId(2),
                                                                vtkOut.GetCell(i).GetPointIds().GetId(5)]
            Cells[PointID] = Cells[PointID] + [i]

            PointID = vtkOut.GetCell(i).GetPointIds().GetId(7)
            Neighbors[PointID] = Neighbors[PointID] + [vtkOut.GetCell(i).GetPointIds().GetId(3),
                                                                vtkOut.GetCell(i).GetPointIds().GetId(4),
                                                                vtkOut.GetCell(i).GetPointIds().GetId(6)]
            Cells[PointID] = Cells[PointID] + [i]

        elif vtkOut.GetCellType(i) == 13:
            # Cell is a Wedge (Prism with triangular base)

            PointID = vtkOut.GetCell(i).GetPointIds().GetId(0)
            Neighbors[PointID] = Neighbors[PointID] + [vtkOut.GetCell(i).GetPointIds().GetId(1),
                                                                vtkOut.GetCell(i).GetPointIds().GetId(2),
                                                                vtkOut.GetCell(i).GetPointIds().GetId(3)]
            Cells[PointID] = Cells[PointID] + [i]

            PointID = vtkOut.GetCell(i).GetPointIds().GetId(1)
            Neighbors[PointID] = Neighbors[PointID] + [vtkOut.GetCell(i).GetPointIds().GetId(0),
                                                                vtkOut.GetCell(i).GetPointIds().GetId(2),
                                                                vtkOut.GetCell(i).GetPointIds().GetId(4)]
            Cells[PointID] = Cells[PointID] + [i]

            PointID = vtkOut.GetCell(i).GetPointIds().GetId(2)
            Neighbors[PointID] = Neighbors[PointID] + [vtkOut.GetCell(i).GetPointIds().GetId(0),
                                                                vtkOut.GetCell(i).GetPointIds().GetId(1),
                                                                vtkOut.GetCell(i).GetPointIds().GetId(5)]
            Cells[PointID] = Cells[PointID] + [i]

            PointID = vtkOut.GetCell(i).GetPointIds().GetId(3)
            Neighbors[PointID] = Neighbors[PointID] + [vtkOut.GetCell(i).GetPointIds().GetId(0),
                                                                vtkOut.GetCell(i).GetPointIds().GetId(4),
                                                                vtkOut.GetCell(i).GetPointIds().GetId(5)]
            Cells[PointID] = Cells[PointID] + [i]

            PointID = vtkOut.GetCell(i).GetPointIds().GetId(4)
            Neighbors[PointID] = Neighbors[PointID] + [vtkOut.GetCell(i).GetPointIds().GetId(1),
                                                                vtkOut.GetCell(i).GetPointIds().GetId(3),
                                                                vtkOut.GetCell(i).GetPointIds().GetId(5)]
            Cells[PointID] = Cells[PointID] + [i]

            PointID = vtkOut.GetCell(i).GetPointIds().GetId(5)
            Neighbors[PointID] = Neighbors[PointID] + [vtkOut.GetCell(i).GetPointIds().GetId(2),
                                                                vtkOut.GetCell(i).GetPointIds().GetId(3),
                                                                vtkOut.GetCell(i).GetPointIds().GetId(4)]
            Cells[PointID] = Cells[PointID] + [i]

        elif vtkOut.GetCellType(i) == 14:
            # Cell is a Pyramid

            PointID = vtkOut.GetCell(i).GetPointIds().GetId(0)
            Neighbors[PointID] = Neighbors[PointID] + [vtkOut.GetCell(i).GetPointIds().GetId(1),
                                                                vtkOut.GetCell(i).GetPointIds().GetId(3),
                                                                vtkOut.GetCell(i).GetPointIds().GetId(4)]
            Cells[PointID] = Cells[PointID] + [i]

            PointID = vtkOut.GetCell(i).GetPointIds().GetId(1)
            Neighbors[PointID] = Neighbors[PointID] + [vtkOut.GetCell(i).GetPointIds().GetId(0),
                                                                vtkOut.GetCell(i).GetPointIds().GetId(2),
                                                                vtkOut.GetCell(i).GetPointIds().GetId(4)]
            Cells[PointID] = Cells[PointID] + [i]

            PointID = vtkOut.GetCell(i).GetPointIds().GetId(2)
            Neighbors[PointID] = Neighbors[PointID] + [vtkOut.GetCell(i).GetPointIds().GetId(1),
                                                                vtkOut.GetCell(i).GetPointIds().GetId(3),
                                                                vtkOut.GetCell(i).GetPointIds().GetId(4)]
            Cells[PointID] = Cells[PointID] + [i]

            PointID = vtkOut.GetCell(i).GetPointIds().GetId(3)
            Neighbors[PointID] = Neighbors[PointID] + [vtkOut.GetCell(i).GetPointIds().GetId(0),
                                       vtkOut.GetCell(i).GetPointIds().GetId(2),
                                       vtkOut.GetCell(i).GetPointIds().GetId(4)]
            Cells[PointID] = Cells[PointID] + [i]

            PointID = vtkOut.GetCell(i).GetPointIds().GetId(4)
            Neighbors[PointID] = Neighbors[PointID] + [vtkOut.GetCell(i).GetPointIds().GetId(0),
                                       vtkOut.GetCell(i).GetPointIds().GetId(1),
                                       vtkOut.GetCell(i).GetPointIds().GetId(2),
                                       vtkOut.GetCell(i).GetPointIds().GetId(3)]

            Cells[PointID] = Cells[PointID] + [i]


        else:
            print(vtkOut.GetCellType(i), " is not a known type")

    # Remove dulicates
    # Sbagliato, devo farlo per forza per ogni punto dell'elemento
    # Aggiungere un ciclo for che cicla sui punti dell'elemento
    for nPoint in range(vtkOut.GetNumberOfPoints()):
        Neighbors[nPoint] = list(dict.fromkeys(Neighbors[nPoint]))
        Cells[nPoint] = list(dict.fromkeys(Cells[nPoint]))


    return Neighbors, Cells



def computeAreaOrVolume(vtkOut, AllCellVolumesOrAreas):

    CellAreasOrVolumes = np.zeros((vtkOut.GetNumberOfPoints(),))

    print("Computing Cells volumes...")

    dummyTetra = vtkTetra()

    for i in range(vtkOut.GetNumberOfPoints()):

        if i%10000 == 0:
            print("Analyzing cells connected to point ", i)

        CellAreasOrVolumes[i] = np.sum(AllCellVolumesOrAreas[Cells[i]])/len(Cells[i])

    return CellAreasOrVolumes

def computeCellVolumesOrAreas(vtkOut, coords):

    CellAreasOrVolumes = [0] * vtkOut.GetNumberOfCells()

    print("Computing Cells volumes...")

    dummyTetra = vtkTetra()

    for i in range(vtkOut.GetNumberOfCells()):

        if i%10000 == 0:
            print("Analyzing cells ", i)
        VolumeOrArea = 0

        if vtkOut.GetCellType(i) == 5:
            # Cell is a triangle
            VolumeOrArea = vtkOut.GetCell(i).ComputeArea()
        elif vtkOut.GetCellType(i) == 9:
            # Cell is a Quad
            VolumeOrArea = computeQuadArea(vtkOut.GetCell(i), coords)
        elif vtkOut.GetCellType(i) == 10:
            # Cell is a Tetra
            p1 = coords[vtkOut.GetCell(i).GetPointIds().GetId(0)]
            p2 = coords[vtkOut.GetCell(i).GetPointIds().GetId(1)]
            p3 = coords[vtkOut.GetCell(i).GetPointIds().GetId(2)]
            p4 = coords[vtkOut.GetCell(i).GetPointIds().GetId(3)]
            VolumeOrArea = dummyTetra.ComputeVolume(p1, p2, p3, p4)
        elif vtkOut.GetCellType(i) == 12:
            # Cell is a Hexaedron
            VolumeOrArea = computeHexaVolume(vtkOut.GetCell(i), coords, dummyTetra)
        elif vtkOut.GetCellType(i) == 13:
            # Cell is a Wedge (Prism with triangular base)
            VolumeOrArea = computeWedgeVolume(vtkOut.GetCell(i), coords, dummyTetra)
        elif vtkOut.GetCellType(i) == 14:
            # Cell is a Pyramid
            VolumeOrArea = computePyramidVolume(vtkOut.GetCell(i), coords, dummyTetra)

        CellAreasOrVolumes[i] = VolumeOrArea

    return np.array(CellAreasOrVolumes)

def computeQuadArea(quad, coords):

    CoordsPointID1 = coords[quad.GetPointIds().GetId(0)]
    CoordsPointID2 = coords[quad.GetPointIds().GetId(1)]
    CoordsPointID3 = coords[quad.GetPointIds().GetId(2)]

    a = np.linalg.norm(CoordsPointID2 - CoordsPointID1)
    b = np.linalg.norm(CoordsPointID3 - CoordsPointID2)
    c = np.linalg.norm(CoordsPointID3 - CoordsPointID1)

    s = (a+b+c)/2

    Area = math.sqrt(s * (s-a) * (s-b) * (s-c))

    CoordsPointID1 = coords[quad.GetPointIds().GetId(2)]
    CoordsPointID2 = coords[quad.GetPointIds().GetId(3)]
    CoordsPointID3 = coords[quad.GetPointIds().GetId(0)]

    a = np.linalg.norm(CoordsPointID2 - CoordsPointID1)
    b = np.linalg.norm(CoordsPointID3 - CoordsPointID2)
    c = np.linalg.norm(CoordsPointID3 - CoordsPointID1)

    s = (a+b+c)/2

    Area += math.sqrt(s * (s-a) * (s-b) * (s-c))

    return Area

def computeHexaVolume(hexa, coords, dummyTetra):

    Volume = 0.0

    # First approximation: It is just an extrusion at constant base area

    # Create the first wedge

    p1 = coords[hexa.GetPointIds().GetId(0)]
    p2 = coords[hexa.GetPointIds().GetId(1)]
    p3 = coords[hexa.GetPointIds().GetId(4)]
    p4 = coords[hexa.GetPointIds().GetId(5)]
    p6 = coords[hexa.GetPointIds().GetId(6)]

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

def computeWedgeVolume(wedge, coords, dummyTetra):

    Volume = 0.0

    p1 = coords[wedge.GetPointIds().GetId(0)]
    p2 = coords[wedge.GetPointIds().GetId(1)]
    p3 = coords[wedge.GetPointIds().GetId(2)]
    p4 = coords[wedge.GetPointIds().GetId(3)]

    Volume += dummyTetra.ComputeVolume(p1, p2, p3, p4)

    p1 = coords[wedge.GetPointIds().GetId(2)]
    p2 = coords[wedge.GetPointIds().GetId(3)]
    p3 = coords[wedge.GetPointIds().GetId(4)]
    p4 = coords[wedge.GetPointIds().GetId(5)]

    Volume += dummyTetra.ComputeVolume(p1, p2, p3, p4)

    p1 = coords[wedge.GetPointIds().GetId(1)]
    p2 = coords[wedge.GetPointIds().GetId(2)]
    p3 = coords[wedge.GetPointIds().GetId(3)]
    p4 = coords[wedge.GetPointIds().GetId(4)]

    Volume += dummyTetra.ComputeVolume(p1, p2, p3, p4)

    return Volume

def computePyramidVolume(pyra, coords, dummyTetra):

    Volume = 0.0

    p1 = coords[pyra.GetPointIds().GetId(0)]
    p2 = coords[pyra.GetPointIds().GetId(1)]
    p3 = coords[pyra.GetPointIds().GetId(3)]
    p4 = coords[pyra.GetPointIds().GetId(4)]

    Volume += dummyTetra.ComputeVolume(p1, p2, p3, p4)

    p1 = coords[pyra.GetPointIds().GetId(1)]
    p2 = coords[pyra.GetPointIds().GetId(2)]
    p3 = coords[pyra.GetPointIds().GetId(3)]
    p4 = coords[pyra.GetPointIds().GetId(4)]

    Volume += dummyTetra.ComputeVolume(p1, p2, p3, p4)

    return Volume

def computeSpacing(AdaptationMethod, NotAdaptedComplexity, newSpacings):

    ## All of the external lists are callable here
    nPoints2Coarse = 0
    nPoints2Refine = 0
    nPointsHittingMax = 0
    nPointsHittingMin = 0
    AdaptedComplexity = NotAdaptedComplexity

    if AdaptationMethod == "Pointwise":
        for i in range(len(avg_S)):
            coordsHere = coords[i]
            spacing = avg_length[i]
            if isInsideTheBox[i]:
                toWrite[i] = True
                # print("S_thresh1 = ", S_thresh1)
                if avg_S[i] < S_thresh:
                    spacing = avg_length[i]*Eps_Coarse
                    if spacing > maxSpacing:
                        spacing = maxSpacing
                        nPointsHittingMax += 1
                    nPoints2Coarse+=1

                    # XOfCoarsenedNodes.append(coords[i][0])
                    # YOfCoarsenedNodes.append(coords[i][1])

                elif avg_S[i] > S_thresh:
                    spacing = EdgeLengthMaxS[i]*(S_thresh/avg_S[i])**(1.0/p)
                    if spacing < minSpacing:
                        spacing = minSpacing
                        nPointsHittingMin += 1
                    nPoints2Refine+=1

                    # XOfRefinedNodes.append(coords[i][0])
                    # YOfRefinedNodes.append(coords[i][1])

            else:
                toWrite[i] = False

            newSpacings[i] = spacing
            AdaptedComplexity += CellAreasOrVolumes[i]/ (spacing**2)

    elif AdaptationMethod == "Pointwise_OnlyRef":
        coordsHere = coords[ActualIndices2Analyze]
        spacing = avg_length[ActualIndices2Analyze]
        avg_S_Here = avg_S[ActualIndices2Analyze]
        indices = avg_S_Here > S_thresh
        indices2Mod = ActualIndices2Analyze[indices]
        toWrite[indices2Mod] = True
        spacing[indices] = EdgeLengthMaxS[indices2Mod]*((S_thresh*np.ones((len(indices2Mod),))/avg_S[indices2Mod])**(1.0/p))
        indices = spacing < minSpacing
        spacing[indices] = minSpacing
        nPointsHittingMin = len(np.where(indices))-1
        nPoints2Refine = len(indices2Mod)
        newSpacings[ActualIndices2Analyze] = spacing[:]

        AdaptedComplexity += np.sum(CellAreasOrVolumes[ActualIndices2Analyze] / (newSpacings[ActualIndices2Analyze]**nDim), axis=0)

    elif AdaptationMethod == "Pointwise_Clusters":

        # Sort the clusters by value of the threshold
        sort_index = np.argsort(centroid)
        # print(sort_index)
        for i in range(len(avg_S)):
            coordsHere = coords[i]
            spacing = avg_length[i]

            if isInsideTheBox[i]:
                toWrite[i] = True
                currentThreshold = centroid[label[i]]

                # Check if it is lower than the threshold set by the centroid
                if avg_S[i] < currentThreshold:
                    # If it's the one with lowest threshold, then I should just coarsen it
                    if currentThreshold == min(centroid):
                        spacing = avg_length[i]*Eps_Coarse
                        if spacing > maxSpacing:
                            spacing = maxSpacing
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
                        spacing = EdgeLengthMaxS[i]*(centroid[previousThresholdIndex]/avg_S[i])**(1.0/p)
                        if spacing < minSpacing:
                            spacing = minSpacing
                            nPointsHittingMin += 1
                        nPoints2Refine+=1

                        # XOfRefinedNodes.append(coords[i][0])
                        # YOfRefinedNodes.append(coords[i][1])

                else:
                    # If it's greater than current threshold, then I will refine with this threshold in mind
                    spacing = EdgeLengthMaxS[i]*(currentThreshold/avg_S[i])**(1.0/p)
                    if spacing < minSpacing:
                        spacing = minSpacing
                        nPointsHittingMin += 1
                    nPoints2Refine+=1

                    # XOfRefinedNodes.append(coords[i][0])
                    # YOfRefinedNodes.append(coords[i][1])

            else:
                toWrite[i] = False

            newSpacings[i] = spacing
            AdaptedComplexity += CellAreasOrVolumes[i]/ (spacing**2)







    pointsStats = [nPoints2Coarse, nPoints2Refine, nPointsHittingMax, nPointsHittingMin]
    return AdaptedComplexity, pointsStats, newSpacings

filenames = ["flow2Refine.vtu"]
nDim = 2
# AdaptationMethod = "Pointwise"
AdaptationMethod = "Pointwise_OnlyRef"
# AdaptationMethod = "Pointwise_Clusters"
AdaptationField = "Mach"
# AdaptationMethod = "Barbara"
p = 2
k_R = 1.0
Eps_Ref = 1.2
Eps_Ref_2 = 1.4
Eps_Coarse = 1.2
decay = 0.7
maxSpacing = 0.5
minSpacing = 1e-5
ComplexityRatio = 1.1
ScalingCoeffForThreshold = 3
Smoothing = True

Limits = "Box"
xLimMin = -1.0
xLimMax = 5.0
yLimMin = -1
yLimMax = 1
zLimMin = -4
zLimMax = 4

# Limits = "Sphere"
# Center_X = 0.5
# Center_Y = 0.0
# Center_Z = 0.0
# Radius = 0.6

#xLimMin = -20
#xLimMax = 50
# yLimMin = -20
# yLimMax = 10

if len(filenames) == 1:

    gridreader = vtk.vtkXMLUnstructuredGridReader()
    gridreader.SetFileName(filenames[0])
    gridreader.Update()
    vtkOut = gridreader.GetOutput()
    vtkData = vtkOut.GetPoints().GetData()

    print("Extracting Point coordinates...")
    coords = np.array([vtkData.GetTuple3(x)
                          for x in range(vtkData.GetNumberOfTuples())])
    # print(cell_connectivity_matrix)

    print("Compute cell volumes/areas of all of the cells...")
    AllCellVolumesOrAreas = computeCellVolumesOrAreas(vtkOut, coords)

    # Compute neighbors for each point
    print("Extracting Neighbors points...")
    Neighbors, Cells = getNeighborsPoints_New(vtkOut)

    # npNeighbors = np.array(Neighbors, dtype="object")

    print("Extracting Gradient of ", AdaptationField, " in the x-direction...")
    vtkData = vtkOut.GetPointData().GetArray("Gradient_"+AdaptationField+"_X")
    grad_x = np.array([vtkData.GetTuple1(x)
                          for x in range(vtkData.GetNumberOfTuples())])

    print("Extracting Gradient of ", AdaptationField, " in the y-direction...")
    vtkData = vtkOut.GetPointData().GetArray("Gradient_"+AdaptationField+"_Y")
    grad_y = np.array([vtkData.GetTuple1(x)
                          for x in range(vtkData.GetNumberOfTuples())])

    grad_z = [0 for j in range(vtkOut.GetNumberOfPoints())]
    if nDim == 3:
        print("Extracting Gradient of ", AdaptationField, " in the z-direction...")
        vtkData = vtkOut.GetPointData().GetArray("Gradient_"+AdaptationField+"_Z")
        grad_z = np.array([vtkData.GetTuple1(x)
            for x in range(vtkData.GetNumberOfTuples())])

    grad = np.transpose(np.stack((grad_x, grad_y, grad_z)))

    avg_S = np.zeros((vtkOut.GetNumberOfPoints(),))
    avg_length = np.zeros((vtkOut.GetNumberOfPoints(),))
    # xVec = [0 for j in range(vtkOut.GetNumberOfPoints())]
    # yVec = [0 for j in range(vtkOut.GetNumberOfPoints())]
    EdgeLengthMaxS = np.zeros((vtkOut.GetNumberOfPoints(),))
    isInsideTheBox = False*np.ones((vtkOut.GetNumberOfPoints(),))

    # Now compute the average edge length for each point
    print("Computing the error for each point...")
    if AdaptationMethod == "Pointwise" or AdaptationMethod == "Pointwise_OnlyRef" or AdaptationMethod == "Pointwise_Clusters":

        timeCheckInside = 0.0
        timeComputeInitialMetric = 0.0

        start = time.time()

        isInsideTheBox = []

        print("Computing if a point is inside the prescribed bounding box...")

        if Limits == "Box":
            limitsMin = np.array([xLimMin, yLimMin, zLimMin])
            limitsMax = np.array([xLimMax, yLimMax, zLimMax])
            isInsideTheBox = np.all(np.logical_and(coords <= np.tile(np.transpose(limitsMax), (vtkOut.GetNumberOfPoints(), 1)), coords >= np.tile(np.transpose(limitsMin), (vtkOut.GetNumberOfPoints(), 1))), axis=1)

        elif Limits == "Sphere":
            CenterCoords = np.tile(np.transpose(np.array([Center_X, Center_Y, Center_Z])), (vtkOut.GetNumberOfPoints(), 1))
            isInsideTheBox = np.sum( (coords-CenterCoords)*(coords-CenterCoords), axis = 1) <= Radius*Radius*np.ones((vtkOut.GetNumberOfPoints(),))

        indices = np.array(range(vtkOut.GetNumberOfPoints()))
        ActualIndices2Analyze = indices[isInsideTheBox]
        ActualIndices2NotAnalyze = indices[np.logical_not(isInsideTheBox)]

        end = time.time()

        timeCheckInside+= end-start

        print("Done! ", len(ActualIndices2Analyze), " found inside the box, ", len(ActualIndices2NotAnalyze), " found outside. Elapsed time = ", timeCheckInside, " s")


        for i in range(vtkOut.GetNumberOfPoints()):
            coordsHere = coords[i]

            tot_S = -1e6
            tot_dist = 0

            if i%10000 == 0:
                print("Analyzing point number ", i)


            start = time.time()

            coordsNeigh = coords[Neighbors[i], :]
            diffCoords = coordsNeigh - np.tile(coordsHere, (len(Neighbors[i]), 1))
            dist = np.resize(np.linalg.norm(diffCoords, axis=1), (len(Neighbors[i]), 1))
            h = diffCoords/np.tile(dist, (1,3))
            delta_grad = np.tile(grad[i, :], (len(Neighbors[i]), 1)) - grad[Neighbors[i],:]

            Sij = (dist**p) * np.absolute(np.resize(np.sum(h*delta_grad, axis=1), (len(Neighbors[i]), 1)))
            indexMaximum = np.argmax(Sij)
            tot_S = Sij[indexMaximum][0]
            EdgeLengthMaxS[i] = dist[indexMaximum][0]

            tot_dist = np.sum(dist)

            avg_length[i] = tot_dist/len(Neighbors[i])
            avg_S[i] = tot_S

            # print(tot_dist)

            end = time.time()
            timeComputeInitialMetric+= end-start


        # print(EdgeLengthMaxS)
        centroid, label = kmeans2(avg_S, 4, minit='points')
        print("clustered thresholds = ", centroid)
        print(np.bincount(label))

        # xVecnp = np.array(xVec)
        # yVecnp = np.array(yVec)
        # plt.scatter(xVecnp[label == 0],yVecnp[label == 0], color='red')
        # plt.scatter(xVecnp[label == 1],yVecnp[label == 1], color='blue')
        # plt.scatter(xVecnp[label == 2],yVecnp[label == 2], color='black')
        # plt.scatter(xVecnp[label == 3],yVecnp[label == 3], color='yellow')
        # plt.show()

        S_thresh = np.mean(avg_S)
        print("np.mean(avg_S) = ", np.mean(avg_S))
        print("np.std(avg_S) = ", np.std(avg_S))


        start= time.time()
        CellAreasOrVolumes = computeAreaOrVolume(vtkOut, AllCellVolumesOrAreas)
        end = time.time()
        print("Done! Elapsed time = ", end-start, " s")

        # Compute current complexity
        CurrentComplexity = np.sum(CellAreasOrVolumes / (avg_length**nDim), axis=0)


        print("CurrentComplexity = ", CurrentComplexity )

        # for i in range(vtkOut.GetNumberOfCells()):

        newSpacings = np.zeros((vtkOut.GetNumberOfPoints(),))
        iterMax = 100
        toll = 1e-3
        error = 1
        iter = 0
        modifyCoeff = False

        NotAdaptedComplexity = 0.0
        if AdaptationMethod == "Pointwise_OnlyRef":

            newSpacings = avg_length.copy()
            NotAdaptedComplexity = np.sum(CellAreasOrVolumes[ActualIndices2NotAnalyze]/(avg_length[ActualIndices2NotAnalyze]**nDim), axis=0)

            print("Residual complexity of points not refined = ", NotAdaptedComplexity)

        XOfCoarsenedNodes = []
        YOfCoarsenedNodes = []
        XOfRefinedNodes = []
        YOfRefinedNodes = []

        while error > toll and iter < iterMax:


            XOfCoarsenedNodes = []
            YOfCoarsenedNodes = []
            XOfRefinedNodes = []
            YOfRefinedNodes = []

            toWrite = False * np.ones((vtkOut.GetNumberOfPoints(),))

            AdaptedComplexity, pointsStats, newSpacings = computeSpacing(AdaptationMethod, NotAdaptedComplexity, newSpacings)
            nPoints2Coarse = pointsStats[0]
            nPoints2Refine = pointsStats[1]
            nPointsHittingMax = pointsStats[2]
            nPointsHittingMin = pointsStats[3]


            if Smoothing:
                AdaptedComplexity = 0
                SmoothedSpacings = newSpacings.copy()
                for i in range(vtkOut.GetNumberOfPoints()):
                    # avgSpacing = newSpacings[i] + np.sum(newSpacings[Neighbors[i]], axis=0)
                    # avgSpacing /= len(Neighbors[i])+1

                    avgSpacing = newSpacings[i]
                    for j in Neighbors[i]:
                        avgSpacing += newSpacings[j]
                    avgSpacing /= len(Neighbors[i])+1
                    SmoothedSpacings[i] = avgSpacing

                AdaptedComplexity += np.sum(CellAreasOrVolumes/ (SmoothedSpacings**nDim), axis=0)
                newSpacings = SmoothedSpacings.copy()


            # print("AdaptedComplexity = ", AdaptedComplexity )

            previousThreshold = S_thresh

            error = abs(AdaptedComplexity - ComplexityRatio*CurrentComplexity) / (ComplexityRatio*CurrentComplexity)
            if AdaptedComplexity > ComplexityRatio*CurrentComplexity:
                if modifyCoeff:
                    ScalingCoeffForThreshold = 1+(ScalingCoeffForThreshold-1)/1.05
                if AdaptationMethod == "Pointwise" or AdaptationMethod == "Pointwise_OnlyRef":
                    S_thresh *= ScalingCoeffForThreshold
                elif AdaptationMethod == "Pointwise_Clusters":
                    centroid *= ScalingCoeffForThreshold

            else :
                modifyCoeff = True
                ScalingCoeffForThreshold = 1+(ScalingCoeffForThreshold-1)/1.05
                if AdaptationMethod == "Pointwise" or AdaptationMethod == "Pointwise_OnlyRef":
                    S_thresh /= ScalingCoeffForThreshold
                elif AdaptationMethod == "Pointwise_Clusters":
                    centroid /= ScalingCoeffForThreshold

            print("Iter ", iter, " previous S_threshold = ", previousThreshold, " current S_threshold = ", S_thresh, " ScalingCoeffForThreshold = ", ScalingCoeffForThreshold)
            print("nPoints2Coarse = ", nPoints2Coarse, " nPoints2Refine = ", nPoints2Refine)
            print("CurrentComplexity = ", CurrentComplexity, " AdaptedComplexity = ", AdaptedComplexity, " error = ", error)
            print("Desired ComplexityRatio = ", ComplexityRatio, " actual ComplexityRatio = ", AdaptedComplexity/CurrentComplexity)
            print()

            iter += 1

        print("Final AdaptedComplexity = ", AdaptedComplexity )
        print("N of iterations = ", iter )
        print("Final threshold", S_thresh )
        print("N of points hitting max spacing = ", nPointsHittingMax, "  N of points hitting min spacing = ", nPointsHittingMin )

            # print(spacing)

        # plt.scatter(XOfCoarsenedNodes,YOfCoarsenedNodes,color='red')
        # plt.scatter(XOfRefinedNodes,YOfRefinedNodes,color='blue')
        # plt.show()

        # cm = plt.cm.get_cmap('RdYlBu')
        # sc = plt.scatter(xVec,yVec, c=avg_S, vmin=4e-6, vmax=4e-5, s=35, cmap=cm)
        # plt.colorbar(sc)
        # plt.show()

        fid = open("PointCloud.dat", 'w')
        for i in range(0,  len(avg_S)):
            if toWrite[i]:
                coordsHere = coords[i]
                fid.write(str(coordsHere[0]) + " " + str(coordsHere[1]) + " " + str(coordsHere[2]) + " " + str(newSpacings[i]) + " " + str(decay) + "\n")

        fid.close()

        print("time for computing initial metric = ", timeComputeInitialMetric, " s")



else :

    print("Unsteady adaptation")
