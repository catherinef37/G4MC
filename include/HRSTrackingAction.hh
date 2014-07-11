// ********************************************************************
// $Id: HRSTrackingAction.hh,v 1.0, 2010/12/26 HRS Exp $
// --------------------------------------------------------------
//

#ifndef HRSTrackingAction_h
#define HRSTrackingAction_h 1

#include "globals.hh"
#include "G4UserTrackingAction.hh"
class HRSTrackingActionMessenger;

class HRSTrackingAction : public G4UserTrackingAction
{
  public:
    HRSTrackingAction();
    virtual ~HRSTrackingAction();

    virtual void PreUserTrackingAction(const G4Track*);
    virtual void PostUserTrackingAction(const G4Track*);

	inline void		SetOutTrackId(G4int val){iOutTrackId=val;strOutParticle="all";}; 
	inline G4int	GetOutTrackId(){return iOutTrackId;}; 
	inline void		SetNoSecondary(G4int val){iNoSecondary=val;}; 
	inline G4int	GetNoSecondary(){return iNoSecondary;};
	inline void		SetOutParticleName(G4String val){strOutParticle=val;iOutTrackId=-1;}; 
	inline G4String	GetOutParticleName(){return strOutParticle;}; 
 
  private:
	  G4int		iOutTrackId;	//default value -1
	  G4int		iNoSecondary;	//default value 1
	  G4String	strOutParticle;	//default value "all"
	  
	  HRSTrackingActionMessenger *messenger;
};

#endif
