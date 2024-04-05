clear
#set -e # interrompe se codice di uscita !=0

function createLink() {
  if [[ -L $2 ]]; then
    unlink $2
  fi # chiude il blocco di istruzioni aperto con "if"
  ln -s $1 $2
}


NRefinements=1

declare -a ComplexityRations=("1.0" "1.1" "1.1" "1.2" "1.3" "1.5") #dichiarazione array in bash

SU2_FOLDER=/home/pagaf/Desktop/Codes/SU2_bin/bin
PARAVIEW_FOLDER=/home/pagaf/ParaView/bin


for (( i = 0; i < NRefinements; i++ )); do

  mkdir -p $i
  cd $i
  createLink ../Meshes/Mesh_${i}.su2 Mesh.su2
  cp ../*.cfg .

  echo '-------------------------------------------'
  echo '╦╔╦╗╔═╗╦═╗╔═╗╔╦╗╦╔═╗╔╗╔  ╔╗╔╦ ╦╔╦╗╔╗ ╔═╗╦═╗
║ ║ ║╣ ╠╦╝╠═╣ ║ ║║ ║║║║  ║║║║ ║║║║╠╩╗║╣ ╠╦╝
╩ ╩ ╚═╝╩╚═╩ ╩ ╩ ╩╚═╝╝╚╝  ╝╚╝╚═╝╩ ╩╚═╝╚═╝╩╚═'
  echo '-------------------------------------------'
  echo '------------------>' $i '<--------------------'
  echo '-------------------------------------------'


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

  echo '-----------------------------------------------------'
  echo '╔═╗╦═╗╔═╗╔═╗╔╦╗╔═╗╔═╗╔═╗╦╔╗╔╔╦╗╔═╗╦  ╔═╗╦ ╦╔╦╗ ╔═╗╦ ╦
║  ╠╦╝║╣ ╠═╣ ║ ║╣ ╠═╝║ ║║║║║ ║ ║  ║  ║ ║║ ║ ║║ ╠═╝╚╦╝
╚═╝╩╚═╚═╝╩ ╩ ╩ ╚═╝╩  ╚═╝╩╝╚╝ ╩ ╚═╝╩═╝╚═╝╚═╝═╩╝o╩   ╩ '
  echo '-----------------------------------------------------'

  sed -i 's/^\s*ComplexityRatio =.*$/ComplexityRatio = '${ComplexityRations[$i]}'/' CreatePointCloud.py
  python3 CreatePointCloud.py


  sed -i 's/^\s*set iter .*$/set iter '$(($i+1))'/' prova.glf
  
  echo '--------------------------'
  echo '╔═╗╦═╗╔═╗╦  ╦╔═╗ ╔═╗╦  ╔═╗
╠═╝╠╦╝║ ║╚╗╔╝╠═╣ ║ ╦║  ╠╣ 
╩  ╩╚═╚═╝ ╚╝ ╩ ╩o╚═╝╩═╝╚  '
  echo '--------------------------'

  pointwise -b prova.glf

  createLink flow2Refine.vtu flow2InterpolateFrom.vtu

  echo '-----------------------------------------------'
  echo '╔╦╗╔═╗╔╦╗╦═╗╔═╗╦ ╦╔═╗╔╦╗╦═╗╔═╗╦  ╦╔═╗╔═╗ ╔═╗╦ ╦
 ║ ║╣  ║ ╠╦╝╠═╣╠═╣║╣  ║║╠╦╝╠═╣║  ║╔═╝║╣  ╠═╝╚╦╝
 ╩ ╚═╝ ╩ ╩╚═╩ ╩╩ ╩╚═╝═╩╝╩╚═╩ ╩╩═╝╩╚═╝╚═╝o╩   ╩ '
  echo '-----------------------------------------------'

  ${PARAVIEW_FOLDER}/pvpython tetrahedralize.py

  createLink $i/restart_flow.dat restart_2InterpFrom.dat
  createLink Meshes/Mesh_$(($i+1)).su2 Mesh_2InterpTo.su2

  echo '----------------------------------'
  echo '╦═╗╔═╗╔═╗╔╦╗╔╗ ╦╔╗╔╔═╗╦═╗╦ ╦╔═╗╦ ╦
╠╦╝║╣ ╠═╣ ║║╠╩╗║║║║╠═╣╠╦╝╚╦╝╠═╝╚╦╝
╩╚═╚═╝╩ ╩═╩╝╚═╝╩╝╚╝╩ ╩╩╚═ ╩o╩   ╩ '
  echo '----------------------------------'

  python3 readBinary.py

done
