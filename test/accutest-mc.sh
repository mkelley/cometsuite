#!/bin/bash
for n in 1 2 3; do
    ../src/rundynamics accutest${n}.par --program="make comet" \
	--pfunc="simple_coma 0.0 7305.0 -999 -999" --nparticles=200 \
	--planetlookup=no
    ../src/xyzinfo accutest${n}.xyz
done
../src/accutest

