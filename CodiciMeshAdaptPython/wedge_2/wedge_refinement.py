# OBLIQUE SHOCK test case
import os
import pandas as pd
import numpy as np
import gmsh
from scipy.interpolate import LinearNDInterpolator

# Find neighbours of a node of gmsh-index id
def neighboring_nodes(id):
    # We look for id in triNodeTags
    neighbours_rows = np.where(triNodeTags==id)
    neighbours_rows = neighbours_rows[0]
    neighbours_tags = triNodeTags[neighbours_rows,:]
    # neighbours_tags is a nX3 matrix which contains the gmsh-indices of the vertices of the traingles surrounding point id.
    # We can observe that neighbours_tags contains in each row id and repeated indices. We extract the non-repeated indices: these are the neighbours of id.
    neighbours_ids = neighbours_tags.reshape((-1,1))[:,0]
    neighbours_ids = [int(x) for x in neighbours_ids]
    neighbours_ids = list(set(neighbours_ids))
    n_list = []
    for nn in neighbours_ids:
        if nn != id:
            n_list.append(nn)
    neighbours_ids = n_list
    # neighbours_list contains the gmsh-indices of the neighbours of id.
    neighbours_ids_py = [int(x-1) for x in neighbours_ids] # py-indices
    # Neighbours coorinates
    neighboursCoord = nodeCoords[neighbours_ids_py,:]
    return neighbours_ids, neighboursCoord

def callback(dim, tag, x, y, z):

    len = len_interp(x,y)
    return len

# Number of adaptation iterations
itmax = 8
# Maximum mesh size
l_max = 0.04
# Minimum mesh size
l_min = 0.0025
# Refinement iteration recap
recap = np.zeros((itmax, 3))
# Compute the initial mesh
os.system("gmsh -2 wedge.geo")
os.system("gmsh -2 -format su2 wedge.geo")
os.system("mv wedge.su2 mesh.su2")

for iter in range(0,itmax):
    if iter==0:
        # Run the simulation without restart
        os.system("mpirun -n 6 SU2_CFD inv_wedge_HLLC.cfg")
    else:
        # Run the simulation with restart
        os.system("mpirun -n 6 SU2_CFD inv_wedge_HLLC_restart.cfg")

    print("-------------------------------------------------------------------------")
    print("-------------------------------------------------------------------------")
    print("--------------------------AUTOMATIC ADAPTATION---------------------------")
    print("-------------------------------------------------------------------------")
    print("-------------------------------------------------------------------------")
    print("Refinement iteration:", iter+1)
    # Run the Paraview Python trace-script on pvpython (terminal)
    print("Importing SU2 solution...")
    os.system("pvpython paraview_trace.py")
    print("Done importing SU2 solution")

    # Read the .csv file and extract density gradient columns
    data = pd.read_csv(r"prova.csv")
    df_X = pd.DataFrame(data, columns=["Gradient_Density_X"])
    grad_rho_x = df_X.to_numpy()
    grad_rho_x = grad_rho_x[:,0]

    df_Y = pd.DataFrame(data, columns=["Gradient_Density_Y"])
    grad_rho_y = df_Y.to_numpy()
    grad_rho_y = grad_rho_y[:,0]

    print("Reading current mesh")
    gmsh.initialize()
    gmsh.open("wedge.msh")

    # Import mesh nodes
    dim = -1
    tag = -1
    nodeTags, nodeCoords, nodeParamCoord = gmsh.model.mesh.getNodes(dim, tag)
    nodeCoords = nodeCoords.reshape((-1,3))
    # The nodes' tags start from 1. Python uses 0-based indices (firt row's index is 0)
    nodeTags_py = [int(x-1) for x in nodeTags]
    # We'll consider two types of indices:
    #   - gmsh-index : starting from 1 onwards (1-based index)
    #   - py-index   : starting from 0 onwards (0-based index)

    # Import mesh elements (triangles)
    eletype = 2 # triangles
    tag = -1 # get them all
    triTags, triNodeTags = gmsh.model.mesh.getElementsByType(eletype, tag)
    # triNodeTags is a gmsh-index
    triNodeTags = triNodeTags.reshape((-1,3)) # reshape into NX3 matrix
    # Each row of triNodeTags contains the 3 gmsh-indices of the vertices of the triangle
    triNodeTags_mod = np.zeros( (len(triNodeTags[:,0]),4) )
    triNodeTags_mod = triNodeTags_mod.astype(int)
    for jj in range(0,len(triNodeTags_mod[:,0])):
        triNodeTags_mod[jj,0] = int(triNodeTags[jj,0])
        triNodeTags_mod[jj,1] = int(triNodeTags[jj,1])
        triNodeTags_mod[jj,2] = int(triNodeTags[jj,2])
        triNodeTags_mod[jj,3] = int(-1)

    # Import mesh elements (quadrilaterals)
    eletype = 3 # quadrilaterals
    tag = -1 # get them all
    quadTags, quadNodeTags = gmsh.model.mesh.getElementsByType(eletype, tag)
    quadNodeTags = quadNodeTags.reshape((-1,4))

    gmsh.finalize()

    # Compute the averaged edge length on each node
    print("Computing new mesh sizes...")
    nodeAvgDist = []
    # Compute the error on each node
    err_nodes = np.zeros((len(nodeTags_py),1))
    for ii in nodeTags_py:
        neighbours_ids, neighboursCoord = neighboring_nodes(ii+1)
        nodeDist = []
        for jj in range(0, len(neighboursCoord[:,0])):
            nodeDist.append( ( np.sqrt( (neighboursCoord[jj,0] - nodeCoords[ii,0])**2 + (neighboursCoord[jj,1] - nodeCoords[ii,1])**2 + (neighboursCoord[jj,2] - nodeCoords[ii,2])**2 ) ) )

        nodeAvgDist.append( np.mean(nodeDist) )
        err_nodes[ii] = np.sqrt(grad_rho_x[ii]**2 + grad_rho_y[ii]**2)*np.mean(nodeDist)

    nodeDist = np.array(nodeDist)
    nodeAvgDist = np.array(nodeAvgDist)

    # Compute the error on each node
    #err_nodes = np.zeros((len(nodeTags_py),1))

    #for ii in nodeTags_py:
    #    err_nodes[ii] = np.sqrt(grad_rho_x[ii]**2 + grad_rho_y[ii]**2)

    # Compute the mean and the standard deviation of the error on each node
    mean_nodes = np.mean(err_nodes)
    std_dev_nodes = np.std(err_nodes)

    # Compute threshold on nodes
    tau_r2_nodes = mean_nodes + std_dev_nodes
    tau_r1_nodes = mean_nodes + std_dev_nodes/2
    tau_c_nodes = mean_nodes

    # Mesh refinement
    newDist = np.zeros((len(nodeTags_py), 1))
    for ii in nodeTags_py:
        if err_nodes[ii] < tau_c_nodes:
            newDist[ii] = min( [1.5*nodeAvgDist[ii], l_max] )
        if err_nodes[ii] < tau_r2_nodes and err_nodes[ii] > tau_r1_nodes:
            newDist[ii] = max( [nodeAvgDist[ii]/1.2, l_min] )
        elif err_nodes[ii] > tau_r2_nodes:
            newDist[ii] = max( [nodeAvgDist[ii]/1.5, l_min] )
        else:
            newDist[ii] = nodeAvgDist[ii]

    len_interp = LinearNDInterpolator(nodeCoords[:,(0,1)], newDist)

    print("Done computing new mesh sizes")
    print("Generating new mesh")
    gmsh.initialize()
    gmsh.merge('wedge.geo')
    gmsh.model.occ.synchronize()

    gmsh.model.mesh.setSizeCallback(callback)

    gmsh.model.mesh.MeshSizeFromPoints = 0
    gmsh.model.mesh.MeshSizeFromCurvature = 0
    gmsh.model.mesh.MeshSizeExtendFromBoundary = 0
    gmsh.model.mesh.generate(2)
    gmsh.write("adapt.msh")
    gmsh.write("adapt.su2")
    gmsh.finalize()

    # Iteration recap
    recap[iter,0] = len(newDist)    # number of nodes
    recap[iter,1] = max(newDist)[0] # max mesh size
    recap[iter,2] = min(newDist)[0] # min mesh size

    # Run the simulation
    os.system("mv adapt.msh wedge.msh")
    os.system("mv adapt.su2 mesh.su2")

# Run the simulation with restart
os.system("mpirun -n 6 SU2_CFD inv_wedge_HLLC_restart.cfg")
#os.system("paraview flow.vtu")
# Export refinement iterations recap
recap_export = pd.DataFrame(recap)
recap_export.to_csv("RefinementRecap.csv")

print("Coarsening threshold:", tau_c_nodes)
print("Intermediate refinement threshold:", tau_r1_nodes)
print("Full refinement threshold:", tau_r2_nodes)
print("Max mesh size:", max(newDist)[0])
print("Min mesh size:", min(newDist)[0])
print("Number of nodes:", len(newDist))