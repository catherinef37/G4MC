// ********************************************************************
//
// $Id: HRSRunAction.hh,v 1.0, 2010/12/26   HRS Exp $
// GEANT4 tag $Name: geant4-09-04 $
//

#ifndef HRSRunAction_h
#define HRSRunAction_h 1

#include "globals.hh"
#include "G4UserRunAction.hh"

class G4Timer;
class G4Run;

class HRSRunAction : public G4UserRunAction
{
  public:
    HRSRunAction();
   ~HRSRunAction();

  public:
    void BeginOfRunAction(const G4Run* aRun);
    void EndOfRunAction(const G4Run* aRun);
    bool bRunActionVerbose;

  private:
    G4Timer* timer;
};


#endif /*HRSRunAction_h*/
