#!/bin/awk -f

NF == 4 {
  rbin = log($3) / log(10.0);
  vbin = log($4) / log(10.0);

  if (rbin < 0) {
    rbin = int(rbin - 1);
  } else {
    rbin = int(rbin);
  }

  if (vbin < 0) {
    vbin = int(vbin - 1);
  } else {
    vbin = int(vbin);
  }

  rcount[rbin]++;
  vcount[vbin]++;
}

END {
  for (bin=-2; bin<8; bin++) {
    print bin, rcount[bin], vcount[bin];
  }
}
