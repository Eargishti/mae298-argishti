#!/usr/bin/bash
#set -x
homeFolder=/home/edgar
parentFolder=$homeFolder/mae298
target1=$parentFolder/gear
target4=$parentFolder/liftdragcells
target2=$homeFolder/OpenFOAM-dev/tutorials/incompressibleFluid/EdgarFoil/constant/polyMesh
#target3=$homeFolder/myFoamUtils/printCelldata




mkdir -p $target1
mkdir -p $target4

cp $parentFolder/faces ./gear/
cp $parentFolder/points ./gear/
cp $parentFolder/cells ./gear/
cp $parentFolder/owner ./gear/
cp $parentFolder/neighbour ./gear/
cp $parentFolder/boundary ./gear/
cp $parentFolder/aerocells ./liftdragcells
cp $parentFolder/cellID ./liftdragcells

cp -r $parentFolder/gear/* $target2/

#cp -r $parentFolder/liftdragcells/* $target3/
