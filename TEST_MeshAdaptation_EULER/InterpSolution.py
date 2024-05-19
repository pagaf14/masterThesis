from InterpSolution_Fun import *
from InterpSolution_Params import options


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
    restartFields.append(a_string.rstrip('\x00'))


restartData = np.fromfile(filenameStart, dtype=np.float64, count=A[1]*A[2], offset=20+33*A[1])
restartData = restartData.reshape((A[1], A[2]), order='F')
# print(restartData)
end = time.time()
print("Done! Elapsed time ", end-start, "s")

if options['interpOnlyNecessaryFields']:

    print("Reading Standard restart...")
    start = time.time()
    # First read the restart file to interp from
    filenameStd="restart_flow_Std.dat"
    A_Std = np.fromfile(filenameStd, dtype=np.int32, count=5)
    restartFields_Std = []
    for i in range(A_Std[1]):
        strbits = np.fromfile(filenameStd, dtype=np.int8, count=33, offset=20+33*i)
        a_string = ''.join([chr(item) for item in strbits])
        # print(i, " ", a_string)
        restartFields_Std.append(a_string.rstrip('\x00'))

    print("Done! Elapsed time ", time.time()-start, "s")

    print("Now extract the field that I am really interested into...")

    interestingFields = []
    for i in range(len(restartFields_Std)):

        interestingFields.append(np.where(np.array(restartFields) == restartFields_Std[i])[0][0])

    restartFields = list(np.array(restartFields)[np.array(interestingFields, dtype=int)])
    restartData = restartData[np.array(interestingFields, dtype=int), :]

    print("Done! Elapsed time ", time.time()-start, "s")


print("Reading final mesh... ")
start = time.time()

if options['MeshType'] == "SU2":
    Points2InterpTo = readMesh_SU2(options['MeshName'], options)
elif options['MeshType'] == "CGNS":
    Points2InterpTo = readMesh_CGNS(options['MeshName'], options)
else:
    print("ERROR! Mesh Format ", options['MeshType'], "is unknown!")
    exit()

end = time.time()
print("Done! Elapsed time ", end-start, "s")



print("Interpolating from mesh with ", len(restartData[0]), " points into the new mesh with ", len(Points2InterpTo), " points...")

# Interpolate the restart
restartDataTo = []
if options['InterpMethod'] == "Barycentric":
    restartDataTo = LinInterp_Mine(restartData, Points2InterpTo, options)
elif options['InterpMethod'] == "NN":
    restartDataTo = NNInterp_Mine(restartData, Points2InterpTo, options)
else:
    print("ERROR! Interpolation method", options['InterpMethod'], "is unknown!")
    exit(2)

if options['enforceFieldsGreaterThanZero']:
    restartDataTo = adjustFields(restartDataTo, np.array(restartFields), options)

# Finally store the new restart.
# First read the restart file to interp from

print("Write interpolated restart...")
start = time.time()
with open("restart_Interp.dat", "wb") as f:

    A[1] = len(restartDataTo)
    A[2] = len(Points2InterpTo)
    A.tofile(f)
    for i in range(A[1]):
        strbits = np.fromfile(filenameStart, dtype=np.int8, count=33, offset=20+33*i)
        strbits.tofile(f)

    restartDataTo = restartDataTo.ravel(order='F')
    restartDataTo.tofile(f)


    end = time.time()
    print("Done! Elapsed time", end-start,"s")


# # check file
# filename="restart_Interp.dat"
# A = np.fromfile(filename, dtype=np.int32, count=5)
# print(A)
# for i in range(A[1]):
#     strbits = np.fromfile(filename, dtype=np.int8, count=33, offset=20+33*i)
#     a_string = ''.join([chr(item) for item in strbits])
#     # print(a_string)
#     restartFields.append(a_string)
#
# print(restartDataTo.reshape((A[1], A[2]), order='F')[5][10000])
# restartData = np.fromfile(filename, dtype=np.float64, count=A[1]*A[2], offset=20+33*A[1]).reshape((A[1], A[2]), order='F')
# print(restartData[5][10000])
