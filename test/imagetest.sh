#!/bin/bash
../src/xyz2fits -p60 syntest.xyz -o syntest.fits
../src/syn2ds9 syntest.xyz -o syntest.reg
