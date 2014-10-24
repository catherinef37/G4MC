#!/bin/csh -fb
#this csh script is trying to red G4 txt output and convert it 
#to G4 input file
#usage: $0:t <infile> <outfile> <pdgcode>
set infile = ($1)
set outfile = ($2)
set pdgcode = ($3)

#echo "###X(mm)    Y(mm)    Z(mm)  Mom(MeV) Theta(deg) Phi(deg)" >! $outfile
#grep '   0 ' $infile | awk '{printf("%8.2f %8.2f %8.2f %8.3f %8.2f %8.2f\n", $2, $3, $4, $7, $8, $9)}' >> $outfile

echo "/mydet/particleNum 1 " > $outfile
echo "/mydet/particle1/particlePDGCode $pdgcode" >> $outfile
grep '   0 ' $infile | awk '{printf("/mydet/position3V %8.2f %8.2f %8.2f mm \n/mydet/particle1/momentum %8.3f MeV \n/mydet/particle1/theta %8.2f deg \n/mydet/particle1/phi %8.2f deg \n/run/beamOn \n\n", $2, $3, $4, $7, $8, $9)}' >> $outfile

