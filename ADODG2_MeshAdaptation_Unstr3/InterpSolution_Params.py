nDim= 2
nNeigh= 1000

MeshName= "Mesh_2InterpTo.su2"
MeshType= "SU2"

# MeshName= "Mesh_2InterpTo.cgns"
# MeshType= "CGNS"

# InterpMethod= "NN"

InterpMethod= "Barycentric"
triaToll= 1e-12 # Tollerance for the check of a point inside a triangle. Keep it really low
nOfCores= 6  # Used only for performance improvements. Number of cores used for the search in the KDTree
             # Change accordingly to how many cores you can use
skipPoints= int(50e3) # Used only for performance improvements. Batch of points to analyze.
                      # Do not change or change if you have enough or low RAM
skipElements= 1 # Used only for performance improvements. Bath of elements to analyze.
                # Do not change or change if you have enough or low RAM
nMinComputations= 10e3 # Used only for performance improvements. Bath of elements to analyze.
                    # Do not change or change if you have enough or low RAM

enforceFieldsGreaterThanZero= True
Fields2Adjust= ['Density', 'Energy', 'Turb_Kin_Energy', 'Eddy_Viscosity']

interpOnlyNecessaryFields= True

options = {'nDim': nDim,
           'nNeigh': nNeigh,
           'InterpMethod': InterpMethod,
           'nOfCores': nOfCores,
           'triaToll': triaToll,
           'skipPoints': skipPoints,
           'nMinComputations': nMinComputations,
           'skipElements': skipElements,
           'MeshName': MeshName,
           'MeshType': MeshType,
           'enforceFieldsGreaterThanZero': enforceFieldsGreaterThanZero,
           'Fields2Adjust': Fields2Adjust,
           'interpOnlyNecessaryFields': interpOnlyNecessaryFields}
