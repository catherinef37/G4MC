######################################################
#use -m64 sompiler option for all code
export COMPILER_BIT=64
#export g2psoft /work/halla/g2p/disk1/centos62
#export g4soft /site/12gev_phys/production/Linux_CentOS6.2-x86_64-gcc4.4.6
export g4soft=/home/Nickie/JLab/Software/geant4.9.6.p03
######################################################
###!set up cernlib

#export CERN /usr/lib64/cernlib
#export CERN_LEVEL 2006
#export CERN_ROOT $CERN/$CERN_LEVEL
#export PATH ${CERN_ROOT}/bin:${PATH}
#if !($?LD_LIBRARY_PATH) then
    #export LD_LIBRARY_PATH ${CERN_ROOT}/lib:/usr/lib64:/usr/local/lib64
#else
    #export LD_LIBRARY_PATH ${CERN_ROOT}/lib:${LD_LIBRARY_PATH}
#endif
export LD_LIBRARY_PATH=/usr/lib64:/usr/local/lib64:$LD_LIBRARY_PATH
######################################################
###!set up root

#export ROOTSYS ${g2psoft}/root-5.34.09-x86_64
#export ROOTLIB $ROOTSYS/lib
#export PATH ${ROOTSYS}/bin:${PATH}
#export LD_LIBRARY_PATH ${ROOTSYS}/lib:${LD_LIBRARY_PATH}
######################################################
###!set up XERCES
export XERCESCROOT=/home/Nickie/JLab/Software/xerces-c-src_2_8_0/
export PATH=$XERCESCROOT/bin:$PATH
export LD_LIBRARY_PATH=$XERCESCROOT/lib:$LD_LIBRARY_PATH
#export XERCESCROOT=${g4soft}/xercesc/3.1.1
#export XERCESCROOT=${g4soft}/xercesc/2.8.0
#export PATH=${XERCESCROOT}/bin:${PATH}
#export LD_LIBRARY_PATH=${XERCESCROOT}/lib:${LD_LIBRARY_PATH}
######################################################
###!set up CLHEP
#export  CLHEP_BASE_DIR/home/Nickie/Downloads/geant4.9.6.p03/source/externals/clhep/
#export  CLHEP_BASE_DIR=${g4soft}/clhep/2.1.3.1
#export  CLHEP_INCLUDE_DIR=${CLHEP_BASE_DIR}/include
#export  CLHEP_LIB_DIR=${CLHEP_BASE_DIR}/lib
#export  CLHEP_LIB=CLHEP

######################################################
###!set up Geant4

export  G4VERSION=9.6.3
export  G4SYSTEM=Linux-g++
#export  G4INSTALL  ${g4soft}/geant4/4.9.6.p03
#export  G4INCLUDE  ${G4INSTALL}/include/Geant4
#export  G4LIB      ${G4INSTALL}/lib64/Geant4-${G4VERSION}

#export G4DATADIR /u/site/12gev_phys/noarch/Geant4-9.6.2/data
#export G4NEUTRONHPDATA   ${G4DATADIR}/data/G4NDL4.2
#export G4LEDATA          ${G4DATADIR}/data/G4EMLOW6.32    
#export G4LEVELGAMMADATA  ${G4DATADIR}/data/PhotonEvaporation2.3
#export G4RADIOACTIVEDATA ${G4DATADIR}/data/RadioactiveDecay3.6
#export G4NEUTRONXSDATA   ${G4DATADIR}/data/G4NEUTRONXS1.2       
#export G4PIIDATA         ${G4DATADIR}/data/G4PII1.3                   
#export G4REALSURFACEDATA ${G4DATADIR}/data/RealSurface1.0     
#export G4SAIDXSDATA      ${G4DATADIR}/data/G4SAIDDATA1.1   

#export G4UI_USE_TCSH  1

#export G4UI_BUILD_XAW_SESSION  1
#export G4UI_USE_XAW  1
#export XAWFLAGS  ""
#export XAWLIBS " -lXaw "

#export G4UI_BUILD_XM_SESSION  1
#export G4UI_USE_XM  1
#export XMFLAGS  ""
#export XMLIBS  " -lXm -lXpm "

export G4UI_BUILD_QT_SESSION=1
export G4UI_USE_QT=1

#export G4VIS_BUILD_DAWN_DRIVER=1
#export G4VIS_USE_DAWN=1

export G4VIS_BUILD_OPENGLX_DRIVER=1
export G4VIS_USE_OPENGLX=1

#export G4VIS_BUILD_OPENGLXM_DRIVER=1
#export G4VIS_USE_OPENGLXM=1

#export G4VIS_BUILD_RAYTRACERX_DRIVER=1
#export G4VIS_USE_RAYTRACERX=1

#export G4VIS_BUILD_VRML_DRIVER=1
#export G4VIS_USE_VRML=1

export G4VIS_BUILD_OPENGLQT_DRIVER=1
export G4VIS_USE_OPENGLQT=1

 
export G4LIB_BUILD_GDML=1
#export G4LIB_BUILD_ZLIB=1
#export G4LIB_BUILD_STATIC=1

#use the default qt of centos6.2
#export QTDIR /usr
#export QTMOC ${QTDIR}/bin/moc-qt4

#export QTDIR ${g4soft}/qt/4.8.5
export QTDIR=/usr
export QTMOC=${QTDIR}/bin/moc
export QTFLAGS="-I${QTDIR}/include -I${QTDIR}/include/QtCore -I${QTDIR}/include/QtGui  -I${QTDIR}/include/QtOpenGL"
export QTLIBS="-L${QTDIR}/lib64 -lQtCore_debug -lQtGui"
export GLQTLIBS="-L${QTDIR}/lib64 -lQtCore_debug -lQtGui -lQtOpenGL"

#export  G4WORKDIR=/u/scratch/$USER/g4work_${G4VERSION}_${OS_NAME}
#if !(-e  $G4WORKDIR) mkdir -p $G4WORKDIR
#if !(-d  $G4WORKDIR) then 
    #export G4WORKDIR=$HOME/.g4work_${G4VERSION}_${OS_NAME}
#endif
#set G4BIN in home dir such that every pc can access, if set it in scratch
#some batch farm node will not see it since they do not share the same scratch disk 
#export G4BIN ${G4WORKDIR}/bin
export G4BIN=$HOME/.g4work_${G4VERSION}_${OS_NAME}/bin

#SET G4TMP so one can recompile the lib with minimum rebuild
export G4TMP=${G4WORKDIR}/tmp

export PATH=${G4BIN}/${G4SYSTEM}:${PATH}
#export LD_LIBRARY_PATH=${CLHEP_LIB_DIR}:${QTDIR}/lib64:${G4LIB}/${G4SYSTEM}:${LD_LIBRARY_PATH}
export LD_LIBRARY_PATH=${QTDIR}/lib64:${G4LIB}/${G4SYSTEM}:${LD_LIBRARY_PATH}
