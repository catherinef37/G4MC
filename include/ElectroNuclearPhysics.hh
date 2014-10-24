//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: ElectroNuclearPhysics.hh,v 1.4 2010/06/03 10:42:44 gunter Exp $
// GEANT4 tag $Name: geant4-09-04-beta-01 $
//
//---------------------------------------------------------------------------
//
// ClassName:   ElectroNuclearPhysics
//
// Author: 2014 Jixie Zhang
//
// Call G4ElectroNuclearBuilder to build electro- or photo- nuclear physics
//----------------------------------------------------------------------------
//
#ifndef ElectroNuclearPhysics_h
#define ElectroNuclearPhysics_h 1

#include "globals.hh"
#include "G4ios.hh"

#include "G4VPhysicsConstructor.hh"
#include "G4ElectroNuclearBuilder.hh"


class ElectroNuclearPhysics : public G4VPhysicsConstructor
{
  public: 
    ElectroNuclearPhysics(G4int verbose = 1);
    ElectroNuclearPhysics(const G4String& name, int verbose = 1);
    virtual ~ElectroNuclearPhysics();

  public: 
    virtual void ConstructParticle();
    virtual void ConstructProcess();


  private:
    void CreateModels();

    G4ElectroNuclearBuilder * theElectroNuclearBuilder;
    
};

#endif

