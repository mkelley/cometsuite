#!/bin/bash
gnuplot scripts/plot_all.gpl -
ds9 syntest.fits -regions syntest.reg
