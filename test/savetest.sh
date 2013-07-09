#!/bin/bash
../src/rundynamics savetest.par --nparticles=10 -o savetestall.xyz
sed -i 's/^SAVE.*/SAVE: radius graindensity beta age origin v_ej r_i v_i t_i r_f v_f t_f/' savetest.par
../src/rundynamics savetest.par --nparticles=10 -o savetestmost.xyz
sed -i 's/^SAVE.*/SAVE: radius graindensity beta age origin r_i v_ej r_f/' savetest.par
../src/rundynamics savetest.par --nparticles=10 -o savetestdefault.xyz
sed -i 's/^SAVE.*/SAVE: beta age origin r_i v_ej r_f/' savetest.par
../src/rundynamics savetest.par --nparticles=10 -o savetestdefault050.xyz
sed -i 's/^SAVE.*/SAVE: graindensity age origin radius label r_i v_i r_f v_f beta v_ej t_i t_f/' savetest.par
../src/rundynamics savetest.par --nparticles=10 -o savetestmixedup.xyz
sed -i 's/^SAVE.*/SAVE: radius graindensity beta age origin v_ej r_i v_i t_i r_f v_f t_f label/' savetest.par
../src/xyzinfo savetestall.xyz | grep -E '(SAVE|Average|Success|^$$)'
../src/xyzinfo savetestmost.xyz | grep -E '(SAVE|Average|Success|^$$)'
../src/xyzinfo savetestdefault.xyz | grep -E '(SAVE|Average|Success|^$$)'
../src/xyzinfo savetestdefault050.xyz | grep -E '(SAVE|Average|Success|^$$)'
../src/xyzinfo savetestmixedup.xyz | grep -E '(SAVE|Average|Success|^$$)'
