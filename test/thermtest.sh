#!/bin/bash
echo > therm.dat
for w in 7 8 9 10 11 12; do \
    ../src/xyz2psd --aper=-1 -a0 -t$w comp-g.xyz -o therm-g-$w.dat;
    ../src/xyz2psd --aper=-1 -a0 -t$w comp-ac.xyz -o therm-ac-$w.dat;
    ../src/xyz2psd --aper=-1 -a0 -t$w comp-ol50.xyz -o therm-ol50-$w.dat;
    echo -n "$w " >> therm.dat;
    awk '/^[^#]/{sum += $2;} END {printf "%e ", sum;}' \
	therm-g-$w.dat >> therm.dat;
    awk '/^[^#]/{sum += $2;} END {printf "%e ", sum;}' \
	therm-ac-$w.dat >> therm.dat;
    awk '/^[^#]/{sum += $2;} END {printf "%e\n", sum;}' \
	therm-ol50-$w.dat >> therm.dat;
done
