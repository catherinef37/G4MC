#!/bin/csh -fb

set infile = ($1)
set outfile = ($2)
echo "###X(mm)    Y(mm)    Z(mm) Ekin(MeV) Mom(MeV) Theta(deg) Phi(deg)" >! $outfile
grep '   0 ' $infile | awk '{printf("%8.2f %8.2f %8.2f %8.3f %8.3f %8.2f %8.2f\n", $2, $3, $4, $5, $7, $8, $9)}' >> $outfile
