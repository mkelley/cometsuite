
mc-extra = --nparticles=200 --planetlookup=no

test: syntest mctest accutest jettest imagetest comptest thermtest plots
test-noplot: syntest mctest accutest jettest imagetest comptest thermtest
test-mc: accutest-mc plots

syntest:
	../src/rundynamics syntest.par
	../src/xyzinfo syntest.xyz

mctest:
	../src/rundynamics mctest.par
	../src/xyzinfo mctest.xyz
	../src/xyz2psd mctest.xyz -o mctest.psd --aper=-1 -b 10 -a0

accutest:
	../src/rundynamics accutest1.par
	../src/xyzinfo accutest1.xyz
	../src/rundynamics accutest2.par
	../src/xyzinfo accutest2.xyz
	../src/rundynamics accutest3.par
	../src/xyzinfo accutest3.xyz
	../src/accutest

jettest:
	../src/jettest

rbvtest:
	../src/rbvtest

comptest:
	../src/rundynamics comptest.par -o comp-g.xyz \
	  --pfunc="simple_coma 0.3 3 -1 4; composition g"

	../src/rundynamics comptest.par -o comp-ac.xyz \
	  --pfunc="simple_coma 0.3 3 -1 4; composition ac"
	../src/rundynamics comptest.par -o comp-ac-d2857.xyz \
	  --pfunc="simple_coma 0.3 3 -1 4; composition ac; fractaldim 2.857"
	../src/rundynamics comptest.par -o comp-ac-d2727.xyz \
	  --pfunc="simple_coma 0.3 3 -1 4; composition ac; fractaldim 2.727"
	../src/rundynamics comptest.par -o comp-ac-d2609.xyz \
	  --pfunc="simple_coma 0.3 3 -1 4; composition ac; fractaldim 2.609"
	../src/rundynamics comptest.par -o comp-ac-d2500.xyz \
	  --pfunc="simple_coma 0.3 3 -1 4; composition ac; fractaldim 2.5"

	../src/rundynamics comptest.par -o comp-ol50.xyz \
	  --pfunc="simple_coma 0.3 3 -1 4; composition ol50"
	../src/rundynamics comptest.par -o comp-ol50-d2857.xyz \
	  --pfunc="simple_coma 0.3 3 -1 4; composition ol50; fractaldim 2.857"
	../src/rundynamics comptest.par -o comp-ol50-d2727.xyz \
	  --pfunc="simple_coma 0.3 3 -1 4; composition ol50; fractaldim 2.727"
	../src/rundynamics comptest.par -o comp-ol50-d2609.xyz \
	  --pfunc="simple_coma 0.3 3 -1 4; composition ol50; fractaldim 2.609"
	../src/rundynamics comptest.par -o comp-ol50-d2500.xyz \
	  --pfunc="simple_coma 0.3 3 -1 4; composition ol50; fractaldim 2.5"

thermtest:
	echo > therm.dat
	for w in 7 8 9 10 11 12; do \
	  ../src/xyz2psd --aper=-1 -a0 -t$$w comp-g.xyz -o therm-g-$$w.dat; \
	  ../src/xyz2psd --aper=-1 -a0 -t$$w comp-ac.xyz -o therm-ac-$$w.dat; \
	  ../src/xyz2psd --aper=-1 -a0 -t$$w comp-ol50.xyz -o therm-ol50-$$w.dat; \
	  echo -n "$$w " >> therm.dat; \
	  awk '/^[^#]/{sum += $$2;} END {printf "%e ", sum;}' therm-g-$$w.dat >> therm.dat; \
	  awk '/^[^#]/{sum += $$2;} END {printf "%e ", sum;}' therm-ac-$$w.dat >> therm.dat; \
	  awk '/^[^#]/{sum += $$2;} END {printf "%e\n", sum;}' therm-ol50-$$w.dat >> therm.dat; \
	done

imagetest:
	../src/xyz2fits -p60 syntest.xyz -o syntest.fits
	../src/syn2ds9 syntest.xyz -o syntest.reg

accutest-mc:
	../src/rundynamics accutest1.par --program="Make Comet" \
	  --pfunc="simple_coma 0.0 7305.0 -999 -999" $(mc-extra)
	../src/xyzinfo accutest1.xyz
	../src/rundynamics accutest2.par --program="Make Comet" \
	  --pfunc="simple_coma 0.0 7305.0 -999 -999" $(mc-extra)
	../src/xyzinfo accutest2.xyz
	../src/rundynamics accutest3.par --program="Make Comet" \
	  --pfunc="simple_coma 0.0 7305.0 -999 -999" $(mc-extra)
	../src/xyzinfo accutest3.xyz
	./accutest

vtest:
	../src/rundynamics vtest1.par

savetest:
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

plots:
	gnuplot scripts/plot_all.gpl -
	ds9 syntest.fits -regions syntest.reg

clean:
	rm -f *~ core.* mctest.xyz syntest.xyz accutest[123].dat	\
		accutest[123].xyz dynamicslog-*.txt *.xyz.log		\
		jettest*.dat mctest.psd syntest.{fits,reg} comp-*.xyz	\
		rbvtest*.dat vtest-[123].dat therm*.dat savetest*.xyz	\
		randtest*.dat syntest.reg syntest.fits

