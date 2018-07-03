#!/bin/bash
#simple bash script that puts the correct script in each folder
#for runs 1-5 with changing AA/SC well depth and ERAD (SC LJ sigma)
BFRAC=0.1
for A in 250
do
	for SCA in 02
	do
		for SCB in 0 001 02 09
		do
		for SIGA in 150
		do
		for SIGB in 100 125 150
			do
			if [ ! -d "$A-$SCA-$SCB-$SIGA-$SIGB" ]; then
				mkdir $A-$SCA-$SCB-$SIGA-$SIGB
			fi
			for RUN in 1 2 3 4 5
			do
				python template_patterning.py $A $SCA $SIGA $SCB $SIGB $RUN $BFRAC
				mv production_template_temp.py $A-$SCA-$SCB-$SIGA-$SIGB/production_run$RUN.py
				sed "s/RUN/$RUN/g" production_template.pbs > $A-$SCA-$SCB-$SIGA-$SIGB/production_run$RUN.pbs
				sed "s/AA/$A/g" -i $A-$SCA-$SCB-$SIGA-$SIGB/production_run$RUN.pbs
				sed "s/SCA/$SCA/g" -i $A-$SCA-$SCB-$SIGA-$SIGB/production_run$RUN.pbs
				sed "s/SCB/$SCB/g" -i $A-$SCA-$SCB-$SIGA-$SIGB/production_run$RUN.pbs
				sed "s/SIGA/$SIGA/g" -i $A-$SCA-$SCB-$SIGA-$SIGB/production_run$RUN.pbs
				sed "s/SIGB/$SIGB/g" -i $A-$SCA-$SCB-$SIGA-$SIGB/production_run$RUN.pbs
				cp rundebughelpers.py $A-$SCA-$SCB-$SIGA-$SIGB
			done	
			done	
		done
		done

	done

done
