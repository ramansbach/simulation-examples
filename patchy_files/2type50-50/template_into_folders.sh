#!/bin/bash
#simple bash script that puts the correct script in each folder
#for runs 1-5 with changing AA/SC well depth and ERAD (SC LJ sigma)

for A in 250 500 750
do
	for SC in 09
	do
		for SIG in 100 125 150 175
		do
			if [ ! -d "$A-$SC-$SIG" ]; then
				mkdir $A-$SC-$SIG
			fi
			for RUN in 1 2 3 4 5
			do
				python template_patterning.py $A $SC $SIG $RUN
				mv production_template_temp.py $A-$SC-$SIG/production_run$RUN.py
				sed "s/RUN/$RUN/g" production_template.pbs > $A-$SC-$SIG/production_run$RUN.pbs
				sed "s/AA/$A/g" -i $A-$SC-$SIG/production_run$RUN.pbs
				sed "s/SC/$SC/g" -i $A-$SC-$SIG/production_run$RUN.pbs
				sed "s/SIG/$SIG/g" -i $A-$SC-$SIG/production_run$RUN.pbs
				cp rundebughelpers.py $A-$SC-$SIG
			done		
		done

	done

done
