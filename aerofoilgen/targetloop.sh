#!/bin/bash
calc_sin() {
	echo "s($1)" | bc -l
}
calc_cos() {
	echo "c($1)" | bc -l
}

PI=3.1415926535
LOWER_RANGE=0
UPPER_RANGE=120
homeFolder=/home/edgar
projectFolder=OpenFOAM-dev/tutorials/incompressibleFluid/EdgarFoil
fileName=U
UFile=$homeFolder/$projectFolder/0/$fileName
mkdir -p graphs

currentDateTime=`date +%Y_%m_%d_%H_%M_%S`
fileCd="./graphs/aoaCd_$currentDateTime"
fileCl="./graphs/aoaCl_$currentDateTime"
fileCm="./graphs/aoaCm_$currentDateTime"
fileL="./graphs/aoaL_$currentDateTime"
fileD="./graphs/aoaD_$currentDateTime"


touch $fileCd
touch $fileCl
touch $fileCm
touch $fileL
touch $fileD



num='[+-]?\d+\.?\d*(?:[eE][+-]?\d+)?'

Umag=26

for ((i=LOWER_RANGE; i<=UPPER_RANGE; i++)); do
	angle_radians=$(echo "scale=10; $i * $PI / (180*4)" | bc -l)

	sine_val=$(calc_sin "$angle_radians")
	cos_val=$(calc_cos "$angle_radians")
	Ux=$(echo "scale=10; $Umag * $cos_val" | bc -l)
	Uy=$(echo "scale=10; $Umag * $sine_val" | bc -l)
    echo "i=$i  angle(rad)=$angle_radians  Ux=$Ux  Uy=$Uy"


	sed -i "s/uniform[^;]*;/uniform ($Ux $Uy 0);/" $UFile
	"$homeFolder/$projectFolder/Allrun"
	DataString=$(foamPostProcess -solver incompressibleFluid -latestTime)

	PressureString=$(echo "$DataString" | grep -A2 "sum of forces" | grep "pressure")
	ViscousString=$(echo "$DataString" | grep -A2 "sum of forces" | grep " viscous")
	read Px Py Pz <<< $(echo "$PressureString" | perl -ne 'print "$+{Px} $+{Py} $+{Pz}" if /\((?<Px>'"$num"')\s+(?<Py>'"$num"')\s+(?<Pz>'"$num"')/')

	read Vcx Vcy Vcz <<< $(echo "$ViscousString" | perl -ne 'print "$+{Vcx} $+{Vcy} $+{Vcz}" if /\((?<Vcx>'"$num"')\s+(?<Vcy>'"$num"')\s+(?<Vcz>'"$num"')/')
    L=$(echo "scale=10; $Py + $Vcy" | bc -l)
	D=$(echo "scale=10; $Px + $Vcx" | bc -l)


    Cmoment=$(echo "$DataString" | grep -A3 "forceCoeffs write"  | grep "Cm" | cut -d'=' -f2 | cut -d' ' -f2)
    Cdrag=$(echo "$DataString" | grep -A3 "forceCoeffs write"  | grep "Cd" | cut -d'=' -f2 | cut -d' ' -f2)
    Clift=$(echo "$DataString" | grep -A3 "forceCoeffs write" | grep "Cl" | cut -d'=' -f2 | cut -d' ' -f2)
	echo -e "$angle_radians\t$L" >> $fileL
	echo -e "$angle_radians\t$D" >> $fileD
	echo -e "$angle_radians\t$Clift" >> $fileCl
	echo -e "$angle_radians\t$Cdrag" >> $fileCd
	echo -e "$angle_radians\t$Cmoment" >> $fileCm
    


    echo "Iteration number $i"
	"$homeFolder/$projectFolder/Allclean"
done


#forces forces write:
#    sum of forces:
#        pressure : (-1.321 18.641 -1.8273e-18)
#        viscous  : (0.1467 0.012719 -2.2262e-21)
#        porous   : (0 0 0)
#    sum of moments:
#        pressure : (-0.46602 -0.033026 5.91)
#        viscous  : (-0.00031797 0.0036674 -0.00098073)
#        porous   : (0 0 0)
#
#forceCoeffs forceCoeffs write:
#    Cm    = 0.12885
#    Cd    = -0.025606
#    Cl    = 0.40675
#    Cl(f) = 0.33222
#    Cl(r) = 0.074524





#sed -i 's/uniform[^;]*;/uniform (25.75, 3.66, 0);/' U
#sed MagUinf
#sed Aref

#./b 0 SimpleMeshInput
#foamPostProcess -solver incompressibleFluid -latestTime

#str="pressure : (-1.321 10.34 -1.222e-3)"
#read first second third <<< $(echo "$str" | perl -ne 'print "$+{first} $+{second} $+{third}" if /\((?<first>'"$num"')\s+(?<second>'"$num"')\s+(?<third>'"$num"')/')
#echo "1st=$first; 2nd=$second; 3rd=$third" 1st=-1.321; 2nd=10.34; 3rd=-1.222e-3



#grep -A2 "sum of forces" shfshf.txt | grep "pressure"
