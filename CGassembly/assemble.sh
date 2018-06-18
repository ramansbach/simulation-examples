#!/bin/bash
S=$(($(date +%s%N)/1000000)) #setting a random seed
DXXX=$7 #chemistry
WDIR=`pwd` #current working directory
GFOLD=`dirname "$WDIR"`/simulation_input_files #just a fancy way of finding this folder because we don't know the absolute file path
LFOLD=`dirname "$WDIR"`/python_scripts
STEPS=400000 #number of steps to simulate
########

if [ $1 == 1 ]
then

cp $LFOLD/createMartiniModel.py .
cp $LFOLD/martini22_ff.py .
cp $LFOLD/IO.py .
cp $LFOLD/MAP.py .
cp $LFOLD/SS.py .
cp $LFOLD/FUNC.py .
cp $LFOLD/DFAG.* .

cp $LFOLD/createPeptideTop.py .
python createPeptideTop.py $DXXX

rm createMartiniModel.py
rm martini22_ff.py 
rm IO.py
rm MAP.py
rm SS.py
rm FUNC.py
fi

if [ $2 == 1 ]
then
mkdir 1_solvate
cd 1_solvate

cp $GFOLD/martini.itp .
cp $GFOLD/water.gro .

cp $GFOLD/table* .
cp $GFOLD/residuetypes.dat .
cp $LFOLD/fixtop.py ./fixtop.py
cp ../$DXXX\.top .
cp ../$DXXX\.itp .
cp ../$DXXX\.gro .

editconf -f $DXXX\.gro -o DXXX_box.gro -box 18.652 18.652 18.652 -bt triclinic
genbox -cp DXXX_box.gro -o DXXX_box100.gro -ci $DXXX\.gro -nmol 99 
genbox -cp DXXX_box100.gro -cs water.gro -o DXXX_solv.gro &> genbox.out
python fixtop.py 100 $DXXX

cd ..
fi

if [ $3 == 1 ]
then

mkdir 2_em
cd 2_em

cp ../1_solvate/martini.itp .
cp ../1_solvate/$DXXX\.itp .
cp ../1_solvate/table* .
cp ../1_solvate/CG_$DXXX\.top .
cp ../1_solvate/residuetypes.dat .
cp ../1_solvate/DXXX_solv.gro .
cp $GFOLD/em5.mdp ./em.mdp

gmx grompp -f em.mdp -c DXXX_solv.gro -p CG_$DXXX\.top -o em.tpr 
gmx mdrun -v -s em.tpr -o em.trr -c after_em.gro -g em.log -e em.edr -x em.xtc -tableb table*
rm \#*

cd ..

fi

if [ $4 == 1 ]
then

mkdir 3_eq
cd 3_eq
cp ../1_solvate/martini.itp .
cp ../1_solvate/$DXXX\.itp .
cp ../1_solvate/table* .
cp ../1_solvate/CG_$DXXX\.top .
cp ../1_solvate/residuetypes.dat .
cp ../2_em/after_em.gro .
cp $GFOLD/eq5.mdp ./eq.mdp
cp $WDIR/indinput.txt .

(cat indinput.txt) | gmx make_ndx -f after_em.gro -o index.ndx

sed -i "s/SS/$S/g" eq.mdp

gmx grompp -f eq.mdp -c after_em.gro -p CG_$DXXX\.top -o eq.tpr -n index.ndx
gmx mdrun -v -s eq.tpr -o eq.trr -c after_eq.gro -g eq.log -e eq.edr -x eq.xtc -cpo eq.cpt -tableb table*

rm \#*
cd ..
fi

if [ $5	== 1 ]
then

mkdir 4_production
cd 4_production

cp ../1_solvate/martini.itp .
cp ../1_solvate/$DXXX\.itp .
cp ../1_solvate/table* .
cp ../1_solvate/CG_$DXXX\.top .
cp ../1_solvate/residuetypes.dat .
cp ../3_eq/after_eq.gro .
cp ../3_eq/eq.cpt .
cp ../3_eq/index.ndx .
cp $GFOLD/md5.mdp ./md.mdp

sed -i "s/STEPS/$STEPS/g" md.mdp
fi

if [ $6 == 1 ]
then

gmx grompp -f md.mdp -c after_eq.gro -p CG_$DXXX.\top -o md.tpr -n index.ndx -t eq.cpt
gmx mdrun -v -s md.tpr -o md.trr -c after_md.gro -g md.log -e md.edr -x md.xtc -cpo md.cpt -tableb table*

fi
