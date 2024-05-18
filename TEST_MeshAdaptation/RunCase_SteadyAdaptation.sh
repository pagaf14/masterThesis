#/bin/bash

set -e

function createLink() {
  if [[ -L $2 ]]; then
    unlink $2
  fi
  if [[ -f $2 ]]; then
    rm $2
  fi
  ln -s $1 $2
}

source RunCase_SteadyAdaptation.cfg

sed -i 's/^\s*AdaptField=.*$/AdaptField= "'${AdaptField}'"/' exportFieldsAndIndices_Params.py
sed -i 's/^\s*PressureAdimFactor=.*$/PressureAdimFactor= '${PressureAdimFactor}'/' exportFieldsAndIndices_Params.py
sed -i 's/^\s*AdaptField=.*$/AdaptField= "'${AdaptField}'"/' CreatePointCloud_Params.py
sed -i 's/^\s*nDim=.*$/nDim= '${nDim}'/' exportFieldsAndIndices_Params.py
sed -i 's/^\s*nDim=.*$/nDim= '${nDim}'/' CreatePointCloud_Params.py
sed -i 's/^\s*nDim=.*$/nDim= '${nDim}'/' InterpSolution_Params.py
sed -i 's/^\s*set nDim.*$/set nDim '${nDim}'/' CreateMesh_Params.glf
sed -i 's/^\s*actualNDim=.*$/actualNDim= '${nDim}'/' CreatePointCloud_Params.py
sed -i 's/^\s*ErrorComputationMethod=.*$/ErrorComputationMethod= "'${ErrorComputationMethod}'"/' CreatePointCloud_Params.py

sed -i 's/^\s*set OutputFormat .*$/set OutputFormat "'${MeshOutput}'"/' CreateMesh_Params.glf

mkdir -p Flows

if [[ ${BLTreatment} == 1 ]]
then

  sed -i 's/^\s*wallDistanceMin=.*$/wallDistanceMin= '${wallDistanceMin}'/' CreatePointCloud_Params.py
  sed -i 's/^\s*set BLTreatment .*$/set BLTreatment 1/' CreateMesh_Params.glf

else

  sed -i 's/^\s*set BLTreatment .*$/set BLTreatment 0/' CreateMesh_Params.glf
  sed -i 's/^\s*BLTreatment=.*$/BLTreatment= False/' exportFieldsAndIndices_Params.py

  sed -i 's/^\s*fileNameIn=.*$/fileNameIn= "flow2Refine.vtu"/' exportFieldsAndIndices_Params.py
  sed -i 's/^\s*fileNameOut=.*$/fileNameOut= "flow2Refine.cgns"/' exportFieldsAndIndices_Params.py
  sed -i 's/^\s*AdaptField=.*$/AdaptField= "'${AdaptField}'"/' exportFieldsAndIndices_Params.py
  sed -i 's/^\s*BLTreatment=.*$/BLTreatment= True/' exportFieldsAndIndices_Params.py
  sed -i 's/^\s*fileFormatOut=.*$/fileFormatOut= "CGNS"/' exportFieldsAndIndices_Params.py

  sed -i 's/^\s*fileFormatIn=.*$/fileFormatIn= "CGNS"/' CreatePointCloud_Params.py
  sed -i 's/^\s*fileNameIn=.*$/fileNameIn= "flow2Refine.cgns"/' CreatePointCloud_Params.py
  sed -i 's/^\s*fileNameOut=.*$/fileNameOut= "PointCloud.dat"/' CreatePointCloud_Params.py
  sed -i 's/^\s*maxSpacing=.*$/maxSpacing= '${volumeMaxSpacing}'/' CreatePointCloud_Params.py
  sed -i 's/^\s*minSpacing=.*$/minSpacing= '${volumeMinSpacing}'/' CreatePointCloud_Params.py
  sed -i 's/^\s*BLTreatment=.*$/BLTreatment= False/' CreatePointCloud_Params.py

fi

if [[ ${InterpolateOnlyRequiredFields} -eq 1 ]]
then
  sed -i 's/^\s*interpOnlyNecessaryFields=.*$/interpOnlyNecessaryFields= True/' InterpSolution_Params.py
  # Now create the restart with just the necessary fields

  mkdir -p 0_Std
  cd 0_Std

  cp ../FirstStep.cfg .

  if [[ ! -f restart_flow.dat ]]
  then
    createLink ../Meshes/Mesh_0.su2 Mesh.su2
    sed -i 's/^\s*ITER=.*$/ITER= 0/' FirstStep.cfg
    sed -i 's/^\s*OUTPUT_FILES=.*$/OUTPUT_FILES= (RESTART)/' FirstStep.cfg
    sed -i 's/^\s*VOLUME_OUTPUT=.*$/VOLUME_OUTPUT= (SOLUTION)/' FirstStep.cfg
    mpirun -n 6 $SU2_FOLDER/SU2_CFD -t 1 FirstStep.cfg
  fi

  cd ..
  createLink 0_Std/restart_flow.dat restart_flow_Std.dat

else
  sed -i 's/^\s*interpOnlyNecessaryFields=.*$/interpOnlyNecessaryFields= False/' InterpSolution_Params.py
fi




# Start iterative refinement...
for (( i = startRefIter; i < NRefinements; i++ )); do

  if [[ $i -gt $skipCFDIter ]]
  then

    mkdir -p $i
    cd $i
    createLink ../Meshes/Mesh_${i}.su2 Mesh.su2
    cp ../*.cfg .

    echo $i

    if [[ $i -gt 0 ]]
    then
      sed -i 's/^\s*RESTART_SOL=.*$/RESTART_SOL= YES/' FirstStep.cfg
      sed -i 's/^\s*ITER=.*$/ITER= 5000/' SecondStep.cfg
      sed -i 's/^\s*ITER=.*$/ITER= 5000/' ThirdStep.cfg
      createLink ../restart_Interp.dat restart_flow.dat
    fi

    mpirun -n 6 $SU2_FOLDER/SU2_CFD -t 1 FirstStep.cfg  | tee -a logFirstStep.log
    mpirun -n 6 $SU2_FOLDER/SU2_CFD -t 1 SecondStep.cfg  | tee -a logSecondStep.log
    mpirun -n 6 $SU2_FOLDER/SU2_CFD -t 1 ThirdStep.cfg  | tee -a logThirdStep.log

    cd ..

  fi

  createLink $i/flow.vtu flow2Refine.vtu
  cd Flows
  createLink ../$i/flow.vtu flow_$(printf '%05d' $i).vtu
  createLink ../$i/surface_flow.vtu surface_flow_$(printf '%05d' $i).vtu
  cd ..

  if [[ ${BLTreatment} -eq 1 ]]
  then
    # If NS Simulation

    createLink $i/surface_flow.vtu surface_flow2Refine.vtu
    sed -i 's/^\s*fileNameIn=.*$/fileNameIn= "surface_flow2Refine.vtu"/' exportFieldsAndIndices_Params.py
    if [[ $nDim -eq 2 ]]
    then
    sed -i 's/^\s*fileNameOut=.*$/fileNameOut= "PointFields.csv"/' exportFieldsAndIndices_Params.py
    sed -i 's/^\s*fileFormatOut=.*$/fileFormatOut= "CSV"/' exportFieldsAndIndices_Params.py
    else
      sed -i 's/^\s*fileNameOut=.*$/fileNameOut= "surface_flow2Refine.cgns"/' exportFieldsAndIndices_Params.py
      sed -i 's/^\s*fileFormatOut=.*$/fileFormatOut= "CGNS"/' exportFieldsAndIndices_Params.py
    fi
    sed -i 's/^\s*AdaptField=.*$/AdaptField= "'${SurfAdaptField}'"/' exportFieldsAndIndices_Params.py

    # Set the BL treatment variables
    sed -i 's/^\s*BLTreatment=.*$/BLTreatment= False/' exportFieldsAndIndices_Params.py
    ${PARAVIEW_FOLDER}/pvpython exportFieldsAndIndices.py

    if [[ $nDim -eq 2 ]]
    then
      sed -i 's/^\s*fileNameIn=.*$/fileNameIn= "PointFields.csv"/' CreatePointCloud_Params.py
      sed -i 's/^\s*fileFormatIn=.*$/fileFormatIn= "CSV"/' CreatePointCloud_Params.py
    else
      sed -i 's/^\s*fileNameIn=.*$/fileNameIn= "surface_flow2Refine.cgns"/' CreatePointCloud_Params.py
      sed -i 's/^\s*fileFormatIn=.*$/fileFormatIn= "CGNS"/' CreatePointCloud_Params.py
    fi
    sed -i 's/^\s*fileNameOut=.*$/fileNameOut= "surface_PointCloud.dat"/' CreatePointCloud_Params.py
    sed -i 's/^\s*AdaptField=.*$/AdaptField= "'${SurfAdaptField}'"/' CreatePointCloud_Params.py
    sed -i 's/^\s*nDim=.*$/nDim= '$(($nDim-1))'/' CreatePointCloud_Params.py
    # Remove BL treatment when dealing with the surface
    sed -i 's/^\s*AdaptationMethod=.*$/AdaptationMethod= "'${SurfAdaptationMethod}'"/' CreatePointCloud_Params.py
    sed -i 's/^\s*maxSpacing=.*$/maxSpacing= '${surfMaxSpacing}'/' CreatePointCloud_Params.py
    if [[ $variableMinSurfSpacing -eq 1 ]]
    then
      sed -i 's/^\s*minSpacing=.*$/minSpacing= '${MinSurfSpacings[$i]}'/' CreatePointCloud_Params.py
    else
      sed -i 's/^\s*minSpacing=.*$/minSpacing= '${surfMinSpacing}'/' CreatePointCloud_Params.py
    fi
    sed -i 's/^\s*BLTreatment=.*$/BLTreatment= False/' CreatePointCloud_Params.py
    sed -i 's/^\s*WriteAllPoints=.*$/WriteAllPoints= True/' CreatePointCloud_Params.py
    python3 -u CreatePointCloud.py > Surf_PointCloud.out 2> Surf_PointCloud.err

    # Re-assign standard names and values
    sed -i 's/^\s*fileNameIn=.*$/fileNameIn= "flow2Refine.vtu"/' exportFieldsAndIndices_Params.py
    sed -i 's/^\s*fileNameOut=.*$/fileNameOut= "flow2Refine.cgns"/' exportFieldsAndIndices_Params.py
    sed -i 's/^\s*fileFormatOut=.*$/fileFormatOut= "CGNS"/' exportFieldsAndIndices_Params.py
    sed -i 's/^\s*AdaptField=.*$/AdaptField= "'${AdaptField}'"/' exportFieldsAndIndices_Params.py
    sed -i 's/^\s*BLTreatment=.*$/BLTreatment= True/' exportFieldsAndIndices_Params.py

    sed -i 's/^\s*fileNameIn=.*$/fileNameIn= "flow2Refine.cgns"/' CreatePointCloud_Params.py
    sed -i 's/^\s*fileFormatIn=.*$/fileFormatIn= "CGNS"/' CreatePointCloud_Params.py
    sed -i 's/^\s*fileNameOut=.*$/fileNameOut= "PointCloud.dat"/' CreatePointCloud_Params.py
    sed -i 's/^\s*AdaptField=.*$/AdaptField= "'${AdaptField}'"/' CreatePointCloud_Params.py
    sed -i 's/^\s*maxSpacing=.*$/maxSpacing= '${volumeMaxSpacing}'/' CreatePointCloud_Params.py
    if [[ $variableMinVolSpacing -eq 1 ]]
    then
      sed -i 's/^\s*minSpacing=.*$/minSpacing= '${MinVolSpacings[$i]}'/' CreatePointCloud_Params.py
    else
      sed -i 's/^\s*minSpacing=.*$/minSpacing= '${volumeMinSpacing}'/' CreatePointCloud_Params.py
    fi
    sed -i 's/^\s*WriteAllPoints=.*$/WriteAllPoints= False/' CreatePointCloud_Params.py
    sed -i 's/^\s*nDim=.*$/nDim= '${nDim}'/' CreatePointCloud_Params.py

    # Set the BL treatment variables
    sed -i 's/^\s*BLTreatment=.*$/BLTreatment= True/' CreatePointCloud_Params.py
    sed -i 's/^\s*AdaptationMethod=.*$/AdaptationMethod= "'${AdaptationMethod}'"/' CreatePointCloud_Params.py

  fi

  ${PARAVIEW_FOLDER}/pvpython exportFieldsAndIndices.py

  python3 -u CreatePointCloud.py > Vol_PointCloud.out 2> Vol_PointCloud.err


  #fi


  sed -i 's/^\s*set iter .*$/set iter '$(($i+1))'/' CreateMesh_Params.glf
  if [[ $variableBLGR -eq 1 ]]
  then
    sed -i 's/^\s*set GRBL .*$/set GRBL '${BLGRs[$i]}'/' CreateMesh_Params.glf
  fi
  if [[ $variableBLSpacing -eq 1 ]]
  then
    sed -i 's/^\s*set BLSpacing .*$/set BLSpacing '${BLSpacings[$i]}'/' CreateMesh_Params.glf
  fi

  ${POINTWISE_FOLDER}/pointwise -b CreateMesh.glf

  createLink flow2Refine.vtu flow2InterpolateFrom.vtu
  ${PARAVIEW_FOLDER}/pvpython exportTetrahedralized.py

  createLink $i/restart_flow.dat restart_2InterpFrom.dat
  sed -i 's/^\s*MeshType=.*$/MeshType= "'${MeshOutput}'"/' InterpSolution_Params.py

  if [[ $MeshOutput -eq "SU2" ]]
  then
    sed -i 's/^\s*MeshName=.*$/MeshName= "Mesh_2InterpTo.su2"/' InterpSolution_Params.py
    createLink Meshes/Mesh_$(($i+1)).su2 Mesh_2InterpTo.su2
  else
    sed -i 's/^\s*MeshName=.*$/MeshName= "Mesh_2InterpTo.cgns"/' InterpSolution_Params.py
    createLink Meshes/Mesh_$(($i+1)).cgns Mesh_2InterpTo.cgns
  fi

  python3 -u InterpSolution.py  > InterpSolution.out 2> InterpSolution.err

done
