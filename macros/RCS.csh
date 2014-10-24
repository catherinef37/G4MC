#!/bin/csh

#this script is use to run the G4MC with various beam energies and detector angle
#usage:  runG4MC.csh <beam> <PhotonAngle> <ProtonAngle> [option]
set DEBUG = echo
set DEBUG = ""

if ( $#argv < 3 ) then
    echo "usage: $0:t <beam> <PhotonAngle> <ProtonAngle> [option]"
    echo '       will launch this cmd:'
    echo '       G4MC -b $beam -o nt_g${PhotonAngle}_p${ProtonAngle}_E${beam}.root 0  -m 1 this_wacs_TCS.mac $option'
    exit
endif

set beam = ($1)
set PhotonAngle = ($2)
set ProtonAngle = ($3)
cat _Detector.ini.good >! Detector.ini
echo "VDRotYAngle=${PhotonAngle};" >> Detector.ini
cat _Detector_SBS.ini.good >! Detector_SBS.ini
echo "SuperBigBiteAngle=${ProtonAngle};" >> Detector_SBS.ini
cat _Detector_HMS.ini.good >! Detector_HMS.ini
echo "HMSAngle=${ProtonAngle};" >> Detector_HMS.ini

set thisPhotonAngle = ($PhotonAngle)
if ( $PhotonAngle > 180 ) @ thisPhotonAngle = $PhotonAngle - 360
cat _wacs_TCS.mac.good >! this_wacs_TCS.mac
echo "/mydet/particle1/detectorAngle $thisPhotonAngle deg" >> this_wacs_TCS.mac
if ( $PhotonAngle > 180 ) then
  echo "/run/beamOn 100000" >> this_wacs_TCS.mac
else 
  echo "/run/beamOn 100000" >> this_wacs_TCS.mac
endif

$DEBUG G4MC -b $beam -o nt_g${PhotonAngle}_p${ProtonAngle}_E${beam}.root 0 -m 1 this_wacs_TCS.mac $4

#get the last ouput file 
set rootfile = (`ls -1t  nt_g${PhotonAngle}_p${ProtonAngle}_E${beam}_??.root | head -n 1`)
set gpfile =  nt_gp_g${PhotonAngle}_p${ProtonAngle}_E${beam}.root
if !( -f $gpfile ) set gpfile =  nt_gp_g${PhotonAngle}_p${ProtonAngle}_E${beam}.0.root


if (-f "$rootfile" ) then  
  if (-f ../macros/AnaWACS.cc) then
    set scriptpath = ../macros;
  else
    set scriptpath = /media/DISK500G/work/RunMC/macros
  endif
  $DEBUG root -b -q $rootfile $scriptpath/CheckAngleAcc.C
  $DEBUG root -b -q $rootfile $scriptpath/SolidAcc.C+
  $DEBUG root -b -q $rootfile $scriptpath/AnaWACS.cc+
  if !( -f $gpfile ) set gpfile =  nt_gp_g${PhotonAngle}_p${ProtonAngle}_E${beam}.0.root
  if ( -f $gpfile ) $DEBUG root -b -q $gpfile $scriptpath/gpAngleAcc.C
endif
