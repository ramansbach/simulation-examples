#!/bin/bash
S=$(($(date +%s%N)/1000000)) #setting a random seed
WDIR=`pwd` #current working directory
CPFOLD=`dirname "$WDIR"` #base directory with various input files and topologies
GFOLD=/home/rachael/Analysis_and_run_code/simulations/gromacs_input_and_automation_scripts/DFAG
INDINPUT=NDIinput.txt
FRAMES=500 #has to do with # of steps, should be equivalent to #ps
DXXX=$7 #should be the chemistry from the 1-letter code, eg dfmi
GROF=NDI_stripped.gro
TOP=NDI_stripped.top
RUNFOLDER=/home/rachael/cluster/scratch/ramansbach/aaPMFS/$DXXX/
SUBFOLDER=/home/rachael/cluster/home/mansbac2/coarsegraining/AApmfs/$DXXX
NMOLS=184
GROUPI=15
########
GROUP1=15
GROUP2=16

if [ $1 == 1 ]
then 

mkdir 1_solvate
cd 1_solvate

cp $CPFOLD/opv3Files/charmm27.ff -r .
cp $CPFOLD/opv3Files/residuetypes.dat .
cp $CPFOLD/opv3Files/$DXXX\-opv3.gro ./$DXXX\.gro
cp $CPFOLD/opv3Files/$DXXX\-opv3.itp ./$DXXX\.itp
cp $CPFOLD/opv3Files/$DXXX\-opv3.top ./$DXXX\.top
	editconf -f $DXXX\.gro -o $DXXX\1.gro -align 0 0 1 -center 6 6 5.85 -box 12 12 12 -bt triclinic
	editconf -f $DXXX\.gro -o $DXXX\2.gro -align 0 0 1 -center 6 6 6.15 -box 12 12 12 -bt triclinic
	cat $DXXX\1.gro $DXXX\2.gro > $DXXX\_dimer.gro
	echo "fix $DXXX\_dimer.gro by removing extra middle lines"
	echo "add second protein to .top file"

cd ..
fi
if [ $2 == 1 ]
then 

cd 1_solvate
genbox -cp $DXXX\_dimer.gro -o solvated.gro -p $DXXX\.top -cs spc216.gro
cd ..

mkdir 2_em
cd 2_em

cp ../1_solvate/solvated.gro .
cp ../1_solvate/$DXXX\.top .
cp ../1_solvate/$DXXX\.itp .
cp ../1_solvate/charmm27.ff -r .
cp ../1_solvate/residuetypes.dat .
cp $CPFOLD/em.mdp .

grompp -f em.mdp -c solvated.gro -p $DXXX\.top -o em.tpr 
mdrun -v -s em.tpr -o em.trr -c after_em.gro -g em.log -e em.edr -x em.xtc

cd ..
fi

if [ $3 == 1 ]
then

mkdir 3_eq
cd 3_eq
cp ../2_em/after_em.gro .
cp ../2_em/$DXXX\.top .
cp ../2_em/$DXXX\.itp .
cp $CPFOLD/make_input.py .
cp $CPFOLD/eq.mdp .
cp $CPFOLD/opv3Files/charmm27.ff -r .
cp $CPFOLD/opv3Files/residuetypes.dat .

python make_input.py $NMOLS $GROUPI
(cat indinput.txt) | make_ndx -f after_em.gro -o mono.ndx

sed -i "s/SS/$S/g" eq.mdp

grompp -f eq.mdp -c after_em.gro -p $DXXX\.top -o eq.tpr -n mono.ndx
mdrun -v -s eq.tpr -o eq.trr -c after_eq.gro -g eq.log -e eq.edr -x eq.xtc -cpo eq.cpt

rm \#*

cd ..
fi

if [ $4 == 1 ]
then

mkdir 4_pull
cd 4_pull
cp ../3_eq/$DXXX\.top .
cp ../3_eq/$DXXX\.itp .
cp ../3_eq/mono.ndx .
cp ../3_eq/after_eq.gro .
cp ../3_eq/eq.cpt .
cp $CPFOLD/pull.mdp ./pull.mdp
cp $CPFOLD/umbrella.sge .
cp ../3_eq/charmm27.ff -r .
cp ../3_eq/residuetypes.dat .

sed -i "s/CHEM/$DXXX/g" umbrella.sge

#grompp -f pull.mdp -c after_eq.gro -p $DXXX\.top -o pull.tpr -t eq.cpt -n mono.ndx
#mdrun -v -s pull.tpr -o pull.trr -c after_pull.gro -g pull.log -e pull.edr -x pull.xtc 
cd ..
fi

if [ $5 == 1 ]
then

cd 4_pull

#echo 1 | trjconv -f after_eq.gro -o after_eq_whole.gro -s pull.tpr -n mono.ndx -pbc whole
#echo 1 | trjconv -f pull.xtc -o pull_whole.xtc -s pull.tpr -n mono.ndx -pbc whole

rm \#*

cp /home/rachael/coarsegraining/orient_restrain_mon_PMF/setupUmbrella.py .
python setupUmbrella.py summary_distances.txt 0.1 $GROUP1 $GROUP2 $FRAMES 1 mono.ndx umbrella.sge
cp frame-*.sge $SUBFOLDER
cd ..

fi

if [ $6 == 1 ]
then
mkdir 5_umbrella
cd 5_umbrella

mv ../4_pull/conf* .
cp ../4_pull/$DXXX\.top .
cp ../4_pull/$DXXX\.itp .
cp ../3_eq/mono.ndx .
cp ../4_pull/charmm27.ff -r .
cp ../4_pull/residuetypes.dat .
cp $CPFOLD/umbrella.mdp .
cp $CPFOLD/npt_umbrella.mdp .

cd ..
cp 5_umbrella -r $RUNFOLDER
fi
