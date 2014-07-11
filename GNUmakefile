# --------------------------------------------------------------
# $Id: GNUmakefile,v 1.02 2010/12/25  HRSMC Exp $
#
# GNUmakefile for HRSMC 
# --------------------------------------------------------------

name := G4MC
G4TARGET := $(name)
G4EXLIB := true
#show the cmd line
#CPPVERBOSE = 1

ifndef G4INSTALL
  G4INSTALL = ../../../..
endif

.PHONY: all

all: HRSLIB XSLIB lib FORLIB bin 

#include the libHRSTransport.a
HRSLIBDIR := obj.$(shell arch)#put this # here to avoid space
EXTRALIBS += -l$(G4TARGET)_FOR -LHRSTransport/$(HRSLIBDIR) -lHRSTransport -LXSModel/$(HRSLIBDIR) -lXSModel 

#Need to link against libgfortran if fortran code used
#in 64-bit machine it is libgfortran, in 32-bit machine it is libg2c
ifeq ($(G4SYSTEM),Linux-g++) 
EXTRALIBS +=  /usr/lib64/libgfortran.so.3
endif
ifeq ($(G4SYSTEM),Darwin-g++) 
#make sure you have created a soft link of libgfortran.a in /usr/lib
#EXTRALIBS += /usr/local/lib/libgfortran.3.dylib
EXTRALIBS += -lgfortran
#EXTRALIBS += -L/sw/lib -lgfortran
endif

#add root package
ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)
ROOTGLIBS    := $(shell root-config --glibs)

CPPFLAGS  += $(ROOTCFLAGS) -IHRSTransport -IXSModel 

EXTRALIBS += $(ROOTLIBS) 

#using g4 compiling framework

include $(G4INSTALL)/config/binmake.gmk

include forrules.gmk


HRSLIB:
	@echo Trying to compile HRSTransport/$(HRSLIBDIR)/libHRSTransport.a ...
	@make -C HRSTransport lib 
	@if [ -e HRSTransport/$(HRSLIBDIR)/libHRSTransport.a ]; then echo "Finish creating lib" HRSTransport/$(HRSLIBDIR)/libHRSTransport.a ; fi
	@if [ ! -e HRSTransport/${HRSLIBDIR}/libHRSTransport.a ]; then echo "***error***! libHRSTransport.a not found!"; exit 1; fi


XSLIB:
	@echo Trying to compile  XSModel/$(HRSLIBDIR)/libXSModel.a ...
	@make -C XSModel lib 
	@if [ -e XSModel/$(HRSLIBDIR)/libXSModel.a ]; then echo "Finish creating lib" XSModel/$(HRSLIBDIR)/libXSModel.a ; fi
	@if [ ! -e XSModel/${HRSLIBDIR}/libXSModel.a ]; then echo "***error***! libXSModel.a not found!"; exit 1; fi

help:
	@echo G4TARGET=[$(G4TARGET)]
	@echo HRSLIBDIR=[$(HRSLIBDIR)]
	@echo CPPFLAGS=[$(CPPFLAGS)] 
	@echo objects=[$(objects)]
	@echo LDFLAGS=[$(LDFLAGS)]
	@echo LDLIBS=[$(LDLIBS)]
	@echo G4TMPDIR=[$(G4TMPDIR)]
	@echo FSOURCE=[$(FSOURCE)] 
	@echo FOBJS=[$(FOBJS)] 
	@echo FFLAGS=[$(FFLAGS)] 
	@echo EXTRALIBS=[$(EXTRALIBS)] 

distclean:
	@make -C XSModel clean
	@make -C HRSTransport clean
	@make clean
