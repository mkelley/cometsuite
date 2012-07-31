#!/bin/bash
for n in 1 2 3; do
    ../src/rundynamics accutest${n}.par
    ../src/xyzinfo accutest${n}.xyz
done
../src/accutest

