#!/bin/bash
../src/rundynamics mctest.par
../src/xyzinfo mctest.xyz
../src/xyz2psd mctest.xyz -o mctest.psd --aper=-1 -b10 -a0
