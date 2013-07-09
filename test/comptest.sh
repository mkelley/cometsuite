#!/bin/bash
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
