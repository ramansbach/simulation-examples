source active py35
for SCB in 100 125 150 175; do
	cd 250-09-09-100-$SCB
		for RUN in 1 2 3 4 5; do
			python production_run$RUN.py
		done
	cd ..
done
