
set -e

function createLink() {
  if [[ -L $2 ]]; then
    unlink $2
  fi
  ln -s $1 $2
}

NRefinements=6
declare -a ComplexityRations=("1.0" "1.1" "1.1" "1.2" "1.3" "1.5")
SU2_FOLDER=/media/rausa/4TB/SU2_Versions/SU2_withGradients/bin
PARAVIEW_FOLDER=/home/rausa/Software/ParaView-5.11.0-RC1-MPI-Linux-Python3.9-x86_64/bin


for (( i = 0; i < NRefinements; i++ )); do

  #if [[ $i > 0 ]]
  #then
  mkdir -p $i
  cd $i
  createLink ../Meshes/Mesh_${i}.su2 Mesh.su2
  cp ../*.cfg .

  echo $i

  if [[ $i > 0 ]]
  then
    sed -i 's/^\s*RESTART_SOL=.*$/RESTART_SOL= YES/' FirstStep.cfg
    sed -i 's/^\s*ITER=.*$/ITER= 500/' SecondStep.cfg
    sed -i 's/^\s*ITER=.*$/ITER= 500/' ThirdStep.cfg
    createLink ../restart_Interp.dat restart_flow.dat
  fi

  mpirun -n 6 $SU2_FOLDER/SU2_CFD -t 1 FirstStep.cfg  | tee -a logFirstStep.log
  mpirun -n 6 $SU2_FOLDER/SU2_CFD -t 1 SecondStep.cfg  | tee -a logSecondStep.log
  mpirun -n 6 $SU2_FOLDER/SU2_CFD -t 1 ThirdStep.cfg  | tee -a logThirdStep.log

  cd ..



  createLink $i/flow.vtu flow2Refine.vtu

  sed -i 's/^\s*ComplexityRatio =.*$/ComplexityRatio = '${ComplexityRations[$i]}'/' CreatePointCloud.py
  python3 CreatePointCloud.py

  #fi


  sed -i 's/^\s*set iter .*$/set iter '$(($i+1))'/' prova.glf

  pointwise -b prova.glf

  createLink flow2Refine.vtu flow2InterpolateFrom.vtu
  ${PARAVIEW_FOLDER}/pvpython tetrahedralize.py

  createLink $i/restart_flow.dat restart_2InterpFrom.dat
  createLink Meshes/Mesh_$(($i+1)).su2 Mesh_2InterpTo.su2

  python3 readBinary.py

done
