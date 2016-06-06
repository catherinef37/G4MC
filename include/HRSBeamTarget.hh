#ifndef __HRSBEAMTARGET_HH
#define __HRSBEAMTARGET_HH

#include "HRStypes.hh"
#include "HRSglobs.hh"
#include "HRSVertex.hh"
#include "HRSDatabase.hh"
#include "G4ThreeVector.hh"
#include <vector>

/*!
     Class that contains information on 
     the beam and target.  It needs to be
     aware of and consistant with what is
     in the geometry.

     This is responsible for:
         Rastering, arbitrary beam angle
	 Sampling along the target
	 Pre-vertex multiple scattering
	 External radiative effects
	 Luminosity calculations
 
     This is implemented in the singleton model

*/

class G4VPhysicalVolume;
class G4Material;
class HRSMultScatt;

class HRSBeamTarget {
    private: 
	static HRSBeamTarget *gSingleton;
	HRSBeamTarget();
	HRSBeamTarget(G4double);

    public:
	static HRSBeamTarget *GetBeamTarget();
	~HRSBeamTarget();

	G4double GetEffLumin();

	void Reset(){ fTargVols.clear(); fMother = NULL; UpdateInfo(); }
	void SetMotherVolume( G4VPhysicalVolume *v ){ fMother = v; UpdateInfo(); }
	void AddVolume( G4VPhysicalVolume *v ){ fTargVols.push_back(v);  UpdateInfo(); }
	void SetTargetPos(G4double z);
	void SetTargetLen(G4double l);
	void SetScatAngle(G4double th);

	void SetBeamCurrent(G4double I){ fBeamCurr = I; }

  HRSDatabase* fDatabase;
  double InterpolateNickie(double);
  double InterpolateNickieCa48(double);

	HRSVertex SampleVertex(G4double);

  G4double fTh;
	G4double fBeamE;
	G4double fBeamCurr;
	G4double fBeamPol;

  G4int mSnakeModel;
  G4double mLHRSMomentum;
  G4double mRHRSMomentum;

  G4double mLSeptumAngle;
  G4double mRSeptumAngle;

	std::vector <G4VPhysicalVolume *> GetTargVols(){ return fTargVols; }

	HRSMultScatt *fMS;

    private:
	std::vector <G4VPhysicalVolume *> fTargVols;
	G4VPhysicalVolume *fMother;

	void UpdateInfo();


	G4double fTotalLength;
	G4double fLH2Length, fZpos, fLH2pos;

	G4Material *fDefaultMat;

	bool fAlreadyWarned;

	
    public:
	// Base position, angle *sampled* info
	G4ThreeVector fVer, fDir;
	G4double fSampE, fRadLen, fSampLen;
	G4double fTravLen;
	G4double fEcut, fEffMatLen;

	// Base position/angle sampling info
	G4double fRasterX, fRasterY;
	G4double fX0, fY0;
	G4double fTh0, fPh0;
	G4double fdTh, fdPh, fCorrTh, fCorrPh;

};


#endif//__HRSBEAMTARGET_HH

