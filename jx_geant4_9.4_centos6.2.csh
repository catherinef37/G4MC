######################################################
#use -m64 sompiler option for all code
setenv COMPILER_BIT 64
setenv g2psoft /work/halla/g2p/disk1/centos62
######################################################
###!set up cernlib

if !($?CERN_ROOT) then
    #setenv CERN ${g2psoft}/cernlib
    setenv CERN /usr/lib64/cernlib
    setenv CERN_LEVEL 2006
    setenv CERN_ROOT $CERN/$CERN_LEVEL
    setenv PATH ${CERN_ROOT}/bin:${PATH}
    if !($?LD_LIBRARY_PATH) then
	setenv LD_LIBRARY_PATH ${CERN_ROOT}/lib:/usr/lib64:/usr/local/lib64
    else
	setenv LD_LIBRARY_PATH ${CERN_ROOT}/lib:${LD_LIBRARY_PATH}
    endif
endif
######################################################
###!set up root

if !($?ROOTSYS) then
    setenv ROOTSYS ${g2psoft}/root_v5.34_x86
    setenv ROOTLIB $ROOTSYS/lib
    setenv PATH ${ROOTSYS}/bin:${PATH}
    setenv LD_LIBRARY_PATH ${ROOTSYS}/lib:${LD_LIBRARY_PATH}
endif
######################################################
###!set up XERCES

setenv XERCESCROOT ${g2psoft}/XercesC_v3.1.1_x86
setenv PATH ${XERCESCROOT}/bin:${PATH}
setenv LD_LIBRARY_PATH ${XERCESCROOT}/lib:${LD_LIBRARY_PATH}
######################################################
###!set up CLHEP

setenv  CLHEP_BASE_DIR     ${g2psoft}/CLHEP_v2.1.0.1_x86
setenv  CLHEP_INCLUDE_DIR  ${CLHEP_BASE_DIR}/include
setenv  CLHEP_LIB_DIR      ${CLHEP_BASE_DIR}/lib
setenv  CLHEP_LIB          CLHEP

######################################################
###!set up Geant4

setenv  G4SYSTEM   Linux-g++
setenv  G4INSTALL  ${g2psoft}/geant4.9.4_centos6.2
setenv  G4INCLUDE  ${G4INSTALL}/include
setenv  G4LIB      ${G4INSTALL}/lib

setenv  G4LEVELGAMMADATA  ${G4INSTALL}/data/PhotonEvaporation2.1
setenv  G4RADIOACTIVEDATA  ${G4INSTALL}/data/RadioactiveDecay3.3
setenv  G4LEDATA  ${G4INSTALL}/data/G4EMLOW6.19
setenv  G4NEUTRONHPDATA  ${G4INSTALL}/data/G4NDL3.14
setenv  G4ABLADATA  ${G4INSTALL}/data/G4ABLA3.0
setenv  G4REALSURFACEDATA  ${G4INSTALL}/data/RealSurface1.0
setenv  G4NEUTRONXSDATA  ${G4INSTALL}/data/G4NEUTRONXS1.0
setenv  G4PIIDATA  ${G4INSTALL}/data/G4PII1.2

setenv G4UI_USE_TCSH  1

setenv G4UI_BUILD_XAW_SESSION  1
setenv G4UI_USE_XAW  1
setenv XAWFLAGS  ""
setenv XAWLIBS " -lXaw "

setenv G4UI_BUILD_XM_SESSION  1
setenv G4UI_USE_XM  1
setenv XMFLAGS  ""
setenv XMLIBS  " -lXm -lXpm "

setenv G4UI_BUILD_QT_SESSION  1
setenv G4UI_USE_QT  1

setenv G4VIS_BUILD_DAWN_DRIVER  1
setenv G4VIS_USE_DAWN  1

setenv G4VIS_BUILD_OPENGLX_DRIVER  1
setenv G4VIS_USE_OPENGLX  1

setenv G4VIS_BUILD_OPENGLXM_DRIVER  1
setenv G4VIS_USE_OPENGLXM  1

setenv G4VIS_BUILD_RAYTRACERX_DRIVER  1
setenv G4VIS_USE_RAYTRACERX  1

setenv G4VIS_BUILD_VRML_DRIVER  1
setenv G4VIS_USE_VRML  1

setenv G4VIS_BUILD_OPENGLQT_DRIVER  1
setenv G4VIS_USE_OPENGLQT  1

setenv G4LIB_BUILD_GDML  1
setenv G4LIB_BUILD_ZLIB  1
setenv G4LIB_BUILD_STATIC  1

setenv OGLHOME  /usr
setenv OGLFLAGS " -I${OGLHOME}/include"
setenv OGLLIBS  " -L${OGLHOME}/lib64 -lGL -lGLU"

setenv OIVHOME ${g2psoft}/G4Vis/Coin_v3.1.3_x86
setenv OIVFLAGS " -I${OIVHOME}/include -I${OIVHOME}/include/Inventor/annex -D_REENTRANT" 
setenv OIVLIBS  " -L${OIVHOME}/lib -lSoXt -lCoin -lGL -lXext -lSM -lICE -lX11 -ldl -lpthread -lXt -lXp -lXi -lXmu -lXpm"


setenv QTDIR            ${g2psoft}/qt_v4.7.1_centos6.2
setenv QTHOME           $QTDIR
setenv PKG_CONFIG_PATH  $QTDIR/lib/pkgconfig
setenv QTMOC            $QTDIR/bin/moc
setenv QTFLAGS  " -I${QTDIR}/include -I${QTDIR}/include/QtCore -I${QTDIR}/include/QtGui  -I${QTDIR}/include/QtOpenGL"
setenv QTLIBS   " -L${QTDIR}/lib -lQtCore -lQtGui -lQtOpenGL"
setenv GLQTLIBS " -L${QTDIR}/lib -lQtCore -lQtGui -lQtOpenGL"

if !($?OS_NAME) setenv OS_NAME Linux64RH6
setenv  G4WORKDIR  /u/scratch/$USER/g4work_${OS_NAME}
#/scratch of batch farm is diff from interactive farm, need to change this one 
if !(-e  $G4WORKDIR) mkdir -p $G4WORKDIR
if !(-d  $G4WORKDIR) then 
setenv G4WORKDIR  $HOME/.g4work_${OS_NAME}
endif
setenv  G4TMP   ${G4WORKDIR}/tmp
#set G4BIN in home dir such that every pc can access, if set it in scratch
#some batch farm node will not see it since they do not share the same scratch disk
setenv  G4BIN   $HOME/.g4work_${OS_NAME}/bin

setenv PATH ${G4BIN}/${G4SYSTEM}:${OIVHOME}/bin:${PATH}
setenv LD_LIBRARY_PATH ${OGLHOME}/lib:${OIVHOME}/lib:${QTDIR}/lib:${CLHEP_LIB_DIR}:${G4LIB}/${G4SYSTEM}:${LD_LIBRARY_PATH}

