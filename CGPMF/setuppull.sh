#!/bin/bash
S=$(($(date +%s%N)/1000000))
WDIR=`pwd`
CPFOLD=`dirname "$WDIR"`/python_scripts
GFOLD=`dirname "$WDIR"`/simulation_input_files
INDINPUT=indinput.txt
NQ=2
Q=0
NION=0
NSTEPS=100000
FRAMES=500 #has to do with # of steps, should be equivalent to the # of ps the pull simulation is run for
QS=N
STRUCT=coildfmi
SCODE=CDFMI
GROF=dimer.gro
RUNFOLDER=/home/rachael/cluster/scratch/ramansbach/coarsegraining/martini-assembly/pmf2ndStructCheck/$STRUCT
SUBFOLDER=/home/rachael/cluster/home/mansbac2/coarsegraining/active_learning/martini-assembly/pmf2ndStructCheck/$STRUCT
LFOLD=/home/rachael/coarsegraining/CG/active_learning/mappingMartini/patchyMap
RFOLD=/home/rachael/coarsegraining/CG/active_learning/martini-assembly
########
GROUP1=14
GROUP2=15

if [ $1 == 1 ]
then
	cp $CPFOLD/fixdimer.py .
	editconf -f dfmi.gro -o dfmi1.gro -align 0 0 1 -center 6 6 5.85 -box 12 12 12 -bt triclinic
	editconf -f dfmi.gro -o dfmi2.gro -align 0 0 1 -center 6 6 6.15 -box 12 12 12 -bt triclinic
	cat dfmi1.gro dfmi2.gro > dfmi_dimer.gro
	python fixdimer.py dfmi_dimer.gro

fi

if [ $2 == 1  ]
then

mkdir 1_solvate
cd 1_solvate

cp $GFOLD/martini.itp .
cp $GFOLD/water.gro .
cp ../dimer.gro .
cp $GFOLD/table* .
cp $GFOLD/residuetypes.dat .
cp $CPFOLD/fixtop.py ./fixtop.py
cp ../dfmi.itp ./DFAG$SCODE\.itp 

genbox -cp dimer.gro -cs water.gro -o DFAG_solv.gro &> genbox.out
python fixtop.py 2 DFAG$SCODE

cd ..
fi

if [ $3 == 1 ]
then 

mkdir 2_em
cd 2_em

cp ../1_solvate/*.itp .
cp ../1_solvate/table* .
cp ../1_solvate/CG_DFAG$SCODE\.top .
cp ../1_solvate/DFAG_solv.gro . 
cp ../1_solvate/residuetypes.dat .
cp $GFOLD/em.mdp .

grompp -f em.mdp -c DFAG_solv.gro -p CG_DFAG$SCODE\.top -o em.tpr 
mdrun -v -s em.tpr -o em.trr -c after_em.gro -g em.log -e em.edr -x em.xtc

rm \#*


cd ..
fi

if [ $4 == 1 ]
then

mkdir 3_eq
cd 3_eq
cp ../1_solvate/martini.itp .
cp ../1_solvate/*.itp .
cp ../1_solvate/table* .
cp ../1_solvate/CG_DFAG$SCODE\.top .
cp ../1_solvate/residuetypes.dat .
cp ../2_em/after_em.gro .
cp $GFOLD/eq.mdp ./eq.mdp
cp  ../indinputdfmi.txt ./indinput.txt

(cat indinput.txt) | make_ndx -f after_em.gro -o index.ndx

sed -i "s/SS/$S/g" eq.mdp

grompp -f eq.mdp -c after_em.gro -p CG_DFAG$SCODE\.top -o eq.tpr -n index.ndx
mdrun -v -s eq.tpr -o eq.trr -c after_eq.gro -g eq.log -e eq.edr -x eq.xtc -cpo eq.cpt


cd ..
fi

if [ $5 == 1 ]
then

mkdir 4_pull
cd 4_pull
cp ../1_solvate/*.itp .
cp ../1_solvate/table* .
cp ../1_solvate/CG_DFAG$SCODE\.top .
cp ../1_solvate/residuetypes.dat .
cp ../3_eq/index.ndx .
cp ../3_eq/after_eq.gro .
cp ../3_eq/eq.cpt .
cp $GFOLD/pull.mdp ./pull.mdp
cp $GFOLD/umbrella.sge .

sed -i "s/CCC/$STRUCT/g" umbrella.sge

sed -i "s/NSTEPS/$NSTEPS/g" pull.mdp
grompp -f pull.mdp -c after_eq.gro -p CG_DFAG$SCODE\.top -o pull.tpr -t eq.cpt -n index.ndx
mdrun -v -s pull.tpr -o pull.trr -c after_pull.gro -g pull.log -e pull.edr -x pull.xtc 
cd ..
fi

if [ $6 == 1 ]
then

cd 4_pull

echo 1 | trjconv -f after_eq.gro -o after_eq_whole.gro -s pull.tpr -n index.ndx -pbc whole
echo 1 | trjconv -f pull.xtc -o pull_whole.xtc -s pull.tpr -n index.ndx -pbc whole

rm \#*

cp $CPFOLD/setupUmbrella.py .
python setupUmbrella.py summary_distances.txt 0.1 $GROUP1 $GROUP2 $FRAMES 1 index.ndx umbrella.sge
cp frame-*.sge $SUBFOLDER
cd ..

fi

if [ $7 == 1 ]
then
mkdir 5_umbrella
cd 5_umbrella

mv ../4_pull/conf* .
cp ../1_solvate/*.itp .
cp ../1_solvate/table* .
cp ../1_solvate/CG_DFAG$SCODE\.top ./CG_DFAG.top
cp ../1_solvate/residuetypes.dat .
cp ../3_eq/index.ndx .
cp $CPFOLD/umbrella.mdp ./umbrella.mdp
cp $CPFOLD/npt_umbrella.mdp ./npt_umbrella.mdp


cd ..
cp 5_umbrella -r $RUNFOLDER
fi
