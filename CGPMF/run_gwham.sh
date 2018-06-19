#!/bin/bash
#bash script that calls g_wham repeatedly after slicing up a trajectory into windows of the given size
echo "Usage: ./run_gwham.sh [interval] [start time] [end time] [tpr names] [pullf names] [pullx names]"
INT=$1 #length of each g_wham call
INIT=$2 #where to begin, may not be 0
TOT=$3 #total length (time) of the profile
TPR=$4 #name of tpr-files.dat
PULLF=$5 #name of pullf-files.dat
PULLX=$6 #name of pullx-files.dat

#run full profile
g_wham -it $TPR -if $PULLF -o prof.xvg -hist hist.xvg -unit kT -b $INIT -e $TOT
for ((int=$INIT;int<$TOT;int+=$INT))
do
	START=$int
	let END=$START+$INT
	#echo "start is $START"
	#echo "end is $END"
	g_wham -it $TPR -if $PULLF -o prof$START-$END\.xvg -hist hist$START-$END\.xvg -unit kT -b $START -e $END
done
