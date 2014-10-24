// ********************************************************************
// $Id: HRSPrimaryGeneratorAction.hh,v 1.0, 2010/12/26 HRS Exp $
// --------------------------------------------------------------
//

#ifndef HRSPrimaryGeneratorAction_h
#define HRSPrimaryGeneratorAction_h 1

#include "globals.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "HRSParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4ThreeVector.hh"
#include "G4VUserPrimaryGeneratorAction.hh"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fstream>
#include <iostream>

#include "BField_Septum.hh"
//#include "Randomize.hh"    //changed by jixie. Use my old Random function
#include "HRSRand.hh"

class HRSPrimaryRootEvent;
class HRSPrimaryRootHisto;
class HRSPrimaryGeneratorMessenger;
using namespace std;

//this number must equal to the one in HRSRootTree
#ifndef MaxPrimaryNum
#define MaxPrimaryNum 8
#endif

class HRSPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
	//the following has been defined for G2P only
	double kTargetXOffset,kTargetYOffset,kTargetZOffset;
	double kLSeptumAngle,kRSeptumAngle;
	double kRasterR;	//=10.0*mm;			//raster size
	double kZLow;		//=-15.0*mm+offset;	//target vertex lower limit
	double kZHigh;		//=15.0*mm+offset;	//target vertex higher limit
	double kZLowTrig;	//=-100.0*cm+offset;//the trigger to use random vertex 
	double kZHighTrig;	//=300.0*cm+offset;	//the trigger to use random vertex 	
	double kHelmCurrentRatio;
	double kHelmXOffset,kHelmYOffset,kHelmZOffset;
	double kEMFieldAtTg[6];

public:
	HRSPrimaryGeneratorAction();
	virtual ~HRSPrimaryGeneratorAction();

	virtual void GeneratePrimaries(G4Event*);

private:
	void   GetMomentum(int index);	//the result will be stored at momentum3V, ekin and momentum_direction;
	void   GetPosition();  //the result will be stored at postion3V
	void   BoNuSProtonEngine(int index);
	void   FastProtonEngine(int index);
	void   HRSElasElectronEngine(int index);
	void   HRSElasNucleusEngine(int index);
	void   HRSQuasiElasElectronEngine(int index);
	void   HRSQuasiElasNucleonEngine(int index);
	void   HRSNtupleEngine(int index);
	void   HRSHistoEngine(int index);
	void   H90UserFitEngine(int index);
	void   ComptonEngine(int index);
	void   BdLEngine(int index);
	void   TwoBodyEngine(int index);
	
	//calculate the BdL (tesla.m) from target to 805mm away along ray with angle of pEndPlaneAngle  
	void  GetHelmBdL(G4ThreeVector &V3BdL, double pEndPlaneAngle);

	void   RandHCSThetaPhi(int i, double &pTheta, double &pPhi);
	void   RandTCSThetaPhi(int i, double &pTheta, double &pPhi);

	int    GenPrimariesFromHisto(G4Event* anEvent);
	int    GenPrimariesFromNtuple(G4Event* anEvent);

	//return the elastic momentum for target X
	double GetElasMomentumX(double pE0, double pTheta, double pMtg);
	double GetElasMomentumX_old(double pE0, double pTheta, double pMtg);
	
	//return the compton momentum for target X
	void GetComptonMomentumX(double pE0, double pTheta_gamma, double pMtg,
		double &pE_gamma, double &pP_x, double &pTheta_x);

	double GetEffectiveTheta(double pZ, double pTheta, double pPhi);
	//return the effective scattering angle if beam is tilted

	//input: y = E_gamma / E_el
	double  GetBremPhotonPol(double y);

private:
	HRSRand mRand;
	HRSPrimaryGeneratorMessenger* gunMessenger;
	HRSParticleGun* particleGun;

	G4ParticleDefinition* electron;
	G4ParticleDefinition* positron;
	G4ParticleDefinition* muonplus;
	G4ParticleDefinition* pion0;
	G4ParticleDefinition* pionplus;
	G4ParticleDefinition* pionminus;
	G4ParticleDefinition* kaon0;
	G4ParticleDefinition* kaonplus;
	G4ParticleDefinition* kaonminus;
	G4ParticleDefinition* proton;

	G4ParticleDefinition* particle[MaxPrimaryNum];

	G4double momentumLow[MaxPrimaryNum],momentumHigh[MaxPrimaryNum],momentum[MaxPrimaryNum];
	G4double thetaLow[MaxPrimaryNum],thetaHigh[MaxPrimaryNum],thetaAngle[MaxPrimaryNum];
	G4double phiLow[MaxPrimaryNum],phiHigh[MaxPrimaryNum],phiAngle[MaxPrimaryNum];
	G4double sigmaMomentum[MaxPrimaryNum];
	G4double sigmaTheta[MaxPrimaryNum];
	G4double sigmaPhi[MaxPrimaryNum];

	//put right arm to be negative or larger than 180 degrees
	G4double detectorAngle[MaxPrimaryNum];  
	//<=0 for HCS, >0  for TCS, 1 is rectangle shape, 2 is eplise
	G4int    randomizeInTCS[MaxPrimaryNum]; 
	G4double outPlaneAngleHigh[MaxPrimaryNum];
	G4double outPlaneAngleLow[MaxPrimaryNum];
	G4double inPlaneAngleHigh[MaxPrimaryNum];
	G4double inPlaneAngleLow[MaxPrimaryNum];

	//available engine are: HRSElasEl,HRSElasNucleus,HRSQuasiElasEl,HRSQusiElasNucleon,
	//BonusProton ,Compton
	G4int    coupleToPrimary[MaxPrimaryNum];  
	G4String primaryEngine[MaxPrimaryNum];  
	G4double beamEnergy;
	G4double targetMass,beamTiltedAngle;  
	G4double targetAtomicNumber; 
	G4double leftHRSMomentum,rightHRSMomentum;
	G4double HMSMomentum;
	G4double incidentEnergy;	 //for photon energy or radiated electron energy
	G4int    helicity;			 //store the helicity status
	G4int	 nThrown;			 //store the number of events that have been thrown 
	G4double polarization[MaxPrimaryNum];  
	
	G4double photonEFracMax,photonEFracMin;  //for compton engine, specify the energy range of the photon

	G4ThreeVector  momentum3V[MaxPrimaryNum];
	G4bool         useMom3V[MaxPrimaryNum];

	G4bool randomizePrimary;
	G4int  particleNum;

	G4int          rasterMode;	
	G4int          vertexFlag;	 //3:fixed point 3V; 2: Random BeamLine 1: Ramdom RZ
	G4double       gunZLow,gunZHigh,gunZ;
	G4double       gunRLow,gunRHigh;
	G4ThreeVector  position3V;
	G4ThreeVector  fixedPointBL3V; //keep the fixed point(x0,y0,z0) in mm
	G4ThreeVector  slopeBL3V;//A,B and Z; if Z out of range then use random Z

	G4bool bUseFastProtonGenerator;

	//Once the user use the following cmds, bUseRootEvent will be set to false;
	//   /mydet/particle#/theta  /mydet/particle#/phi /mydet/particle#/momentum /mydet/particle#/momentum3V  
	G4bool bUseRootEvent;   
	G4bool bUseRootHisto;   
	G4int  iRootEventType;

	HRSPrimaryRootEvent *mRootEvent;
	HRSPrimaryRootHisto *mRootHisto;

	//use to rotate to get the direction ofr elastic electron 
	//need to set the rotation as module member to speed up (avoid reinstance)
	//G4RotationMatrix mRotHRSElasEl;

	//this 2 vector will be used to calculate the effective scattering angle
	//I put it here to avoid contruction and deconstruction, which will improve speed
	G4ThreeVector mBeamV3,mParticleV3;

	//this array will be filled into root tree
	double effectiveTheta[MaxPrimaryNum];

public:

	HRSPrimaryRootHisto* GetRootHisto()
	{
		return mRootHisto;
	};
	HRSPrimaryRootEvent* GetRootEvent()
	{
		return mRootEvent;
	};
	inline void SetMomentumLow(G4int i,G4double val)
	{
		momentumLow[i] = val;
	}
	inline void SetMomentumHigh(G4int i,G4double val)
	{
		momentumHigh[i] = val;
	}
	inline void SetTotalMomentum(G4int i,G4double val)
	{
		momentum[i] = val; useMom3V[i]=false; bUseRootEvent=false;
	}
	inline G4double GetMomentumLow(G4int i) const
	{
		return momentumLow[i]/GeV;
	}
	inline G4double GetMomentumHigh(G4int i) const
	{
		return momentumHigh[i]/GeV;
	}
	inline G4double GetTotalMomentum(G4int i) const
	{
		return momentum[i]/GeV;
	}

	inline void SetSigmaMomentum(G4int i,G4double val)
	{
		sigmaMomentum[i] = val;
	}
	inline G4double GetSigmaMomentum(G4int i) const
	{
		return sigmaMomentum[i]/GeV;
	}

	inline void SetSigmaTheta(G4int i,G4double val)
	{
		sigmaTheta[i] = val;
	}
	inline G4double GetSigmaTheta(G4int i) const
	{
		return sigmaTheta[i]/deg;
	}

	inline void SetSigmaPhi(G4int i,G4double val)
	{
		sigmaPhi[i] = val;
	}
	inline G4double GetSigmaPhi(G4int i) const
	{
		return sigmaPhi[i]/deg;
	}

	inline void SetThetaLow(G4int i,G4double val)
	{
		thetaLow[i] = val;
	}
	inline void SetThetaHigh(G4int i,G4double val)
	{
		thetaHigh[i] = val;
	}
	inline void SetTheta(G4int i,G4double val)
	{
		thetaAngle[i] = val; useMom3V[i]=false; bUseRootEvent=false;
	}
	inline G4double GetThetaLow(G4int i) const
	{
		return thetaLow[i]/deg;
	}
	inline G4double GetThetaHigh(G4int i) const
	{
		return thetaHigh[i]/deg;
	}
	inline G4double GetTheta(G4int i) const
	{
		return thetaAngle[i]/deg;
	}

	inline void SetPhiLow(G4int i,G4double val)
	{
		phiLow[i] = val;
	}
	inline void SetPhiHigh(G4int i,G4double val)
	{
		phiHigh[i] = val;
	}
	inline void SetPhi(G4int i,G4double val)
	{
		phiAngle[i] = val; useMom3V[i]=false; bUseRootEvent=false;
	}
	inline G4double GetPhiLow(G4int i) const
	{
		return phiLow[i]/deg;
	}
	inline G4double GetPhiHigh(G4int i) const
	{
		return phiHigh[i]/deg;
	}
	inline G4double GetPhi(G4int i) const
	{
		return phiAngle[i]/deg;
	}

	//randomize the particle
	inline void SetRandomize(G4bool val)
	{
		randomizePrimary = val;
	}
	inline G4bool GetRandomize() const
	{
		return randomizePrimary;
	}

	//particle num
	inline void SetParticleNum(int val)
	{
		particleNum=val;
	};
	inline G4int GetParticleNum() const
	{
		return particleNum;
	};

	inline void SetGunZLow(G4double val)
	{
		gunZLow = val;
	}
	inline void SetGunZHigh(G4double val)
	{
		gunZHigh = val;
	}
	inline void SetGunZ(G4double val)
	{
		gunZ = val;
	}
	inline G4double GetGunZLow() const
	{
		return gunZLow/mm;
	}
	inline G4double GetGunZHigh() const
	{
		return gunZHigh/mm;
	}
	inline G4double GetGunZ() const
	{
		return gunZ/mm;
	}
	
	inline void SetVertexMode(int val)
	{
		vertexFlag=val;
	}
	inline int GetVertexMode() const
	{
		return vertexFlag;
	}

	inline void SetGunRLow(G4double val)
	{
		gunRLow = val;
	}
	inline G4double GetGunRLow() const
	{
		return gunRLow/mm;
	}
	inline void SetGunRHigh(G4double val)
	{
		gunRHigh = val;
	}
	inline G4double GetGunRHigh() const
	{
		return gunRHigh/mm;
	}

	//engine type
	inline void SetPrimaryEngine(int i,G4String val)
	{
		primaryEngine[i]=val;
	};
	inline G4String GetPrimaryEngine(int i) const
	{
		return primaryEngine[i];
	};
	
	inline void SetRandomizeInTCS(int i,G4int val)
	{
		randomizeInTCS[i]=val;
	};
	inline G4int GetRandomizeInTCS(int i) const
	{
		return randomizeInTCS[i];
	};

	inline void SetDetectorAngle(int i,G4double val)
	{
		detectorAngle[i]=val;
	};
	inline G4double GetDetectorAngle(int i) const
	{
		return detectorAngle[i];
	};

	//theta_tr
	inline void SetOutPlaneAngleHigh(int i,G4double val)
	{
		outPlaneAngleHigh[i]=val;
	};
	inline G4double GetOutPlaneAngleHigh(int i) const
	{
		return outPlaneAngleHigh[i];
	};
	inline void SetOutPlaneAngleLow(int i,G4double val)
	{
		outPlaneAngleLow[i]=val;
	};
	inline G4double GetOutPlaneAngleLow(int i) const
	{
		return outPlaneAngleLow[i];
	};

	//phi_tr
	inline void SetInPlaneAngleHigh(int i,G4double val)
	{
		inPlaneAngleHigh[i]=val;
	};
	inline G4double GetInPlaneAngleHigh(int i) const
	{
		return inPlaneAngleHigh[i];
	};
	inline void SetInPlaneAngleLow(int i,G4double val)
	{
		inPlaneAngleLow[i]=val;
	};
	inline G4double GetInPlaneAngleLow(int i) const
	{
		return inPlaneAngleLow[i];
	};

	//Get particle definition
	inline G4ParticleDefinition* GetParticleDefinition(int i)
	{
		return particle[i];
	};
	//set particle by particle name
	inline void SetParticleDefinition(int i,G4String val)
	{
		particle[i]=G4ParticleTable::GetParticleTable()->FindParticle(val);
	};
	inline G4String GetParticleName(int i) const
	{
		return particle[i]->GetParticleName();
	};
	//set particle by pdg code
	inline void SetParticleDefinition(int i,G4int val)
	{
		particle[i]=G4ParticleTable::GetParticleTable()->FindParticle(val);
	}
	//Get pdg code
	inline G4int GetParticlePDGCode(int i) const
	{
		return particle[i]->GetPDGEncoding();
	}
	//get pdg charge
	inline G4double GetPDGCharge(int i) const
	{
		return particle[i]->GetPDGCharge();
	}
	inline G4double GetPDGMass(int i) const
	{
		return particle[i]->GetPDGMass();
	}

	// momentum 3Vector
	inline void SetMomentum3V(int i,G4ThreeVector val3V)
	{
		momentum3V[i]=val3V;useMom3V[i]=true;bUseRootEvent=false; 
	};
	inline G4ThreeVector GetMomentum3V(int i) const
	{
		return momentum3V[i];
	};

	// Beamline position 3Vector and slope a & b
	inline void SetFixedPointBL3V(G4ThreeVector val3V)
	{
		fixedPointBL3V=val3V;vertexFlag=2;
	};
	inline G4ThreeVector GetFixedPointBL3V() const
	{
		return fixedPointBL3V;
	};
	inline void SetSlopeBL3V(G4ThreeVector val3V)
	{
		slopeBL3V=val3V;vertexFlag=2;
	};
	inline G4ThreeVector GetSlopeBL3V() const
	{
		return slopeBL3V;
	};
	inline void SetPosition3V(G4ThreeVector val3V)
	{
		position3V=val3V;vertexFlag=3;
	};
	inline G4ThreeVector GetPosition3V() const
	{
		return position3V;
	};

	inline void SetRasterMode(int val)
	{
		rasterMode=val;
	};
	inline int GetRasterMode() const
	{
		return rasterMode;
	};

	inline void SetFastPsGenFlag(G4bool val)
	{
		bUseFastProtonGenerator = val;
	}
	inline G4bool GetFastPsGenFlag() const
	{
		return bUseFastProtonGenerator;
	}


	//////////////////////////////////////////////////////
	inline void SetBeamEnergy(G4double val)
	{ 
		beamEnergy=val;
	};
	inline G4double GetBeamEnergy()
	{
		return beamEnergy;
	};
	
	inline void SetTargetMass(G4double val)
	{ 
		targetMass=val;
	};
	inline G4double GetTargetMass()
	{
		return targetMass;
	};
	
	inline void SetTargetAtomicNumber(G4double val)
	{
		targetAtomicNumber=val;
	};
	inline G4double GetTargetAtomicNumber()
	{
		return targetAtomicNumber;
	};

	// momentum 3Vector
	inline void SetBeamNTarget3V(G4ThreeVector val3V)
	{
		beamEnergy=val3V.x()*GeV;targetMass=val3V.y()*GeV;targetAtomicNumber=val3V.z(); 
	};
	inline G4ThreeVector GetBeamNTarget3V()
	{
		return G4ThreeVector(beamEnergy/GeV,targetMass/GeV,targetAtomicNumber); 
	};

	//HRS Momenum
	inline void SetLeftHRSMomentum(G4double val)
	{ 
		leftHRSMomentum=val;
		BField_Septum* pSeptum=BField_Septum::GetInstance();
		pSeptum->SetMomentumL(leftHRSMomentum/GeV);
	};
	inline G4double GetLeftHRSMomentum()
	{
		return leftHRSMomentum;
	};
	inline void SetRightHRSMomentum(G4double val)
	{ 
		rightHRSMomentum=val;
		BField_Septum* pSeptum=BField_Septum::GetInstance();
		pSeptum->SetMomentumR(rightHRSMomentum/GeV);
	};
	inline G4double GetRightHRSMomentum()
	{
		return rightHRSMomentum;
	};

	//read or calculate the effective theta angle for the case that the beam is tilted
	//if it has been calculated, just read the value
	inline double GetThetaEff(int i)
	{
		//in order to speed up this pogram, I put this here
		if(effectiveTheta[i]<-9.99)
		{
			effectiveTheta[i]=GetEffectiveTheta(position3V.z(),momentum3V[i].theta(),momentum3V[i].phi());
		}
		return effectiveTheta[i];
	};
	
	inline void SetCoupleToPrimary(int i, int val)
	{ 
		if(val>i || val<=0) 
		{
			cout<<"*** Cancel coupling for parimary "<<i+1<<" since given parimary ordinal number ="<<val<<" ***"<<endl;
			coupleToPrimary[i]=0;
		}
		else coupleToPrimary[i]=val;
	};
	inline G4double GetCoupleToPrimary(int i)
	{
		return  coupleToPrimary[i];
	};

	inline void SetHMSMomentum(G4double val)
	{ 
		HMSMomentum=val;
		//I am not sure if HMS will use septum field or not
		//should add this feture later
		//BField_Septum* pSeptum=BField_Septum::GetInstance();
		//pSeptum->SetMomentumR(rightHRSMomentum/GeV);
	};
	inline G4double GetHMSMomentum()
	{
		return HMSMomentum;
	};
	//Get polarization 
	inline G4double GetPolarization(int i) 
	{
		return polarization[i];
	}

	inline double GetIncidentEnergy() {return incidentEnergy;};
	inline void SetIncidentEnergy(double v) { incidentEnergy = v;};

	inline double GetHelicity() {return helicity;};
	inline void SetHelicity(double v) { helicity = v;};

	inline double GetPhotonEFracMax(){return photonEFracMax;};
	inline void SetPhotonEFracMax(double v){photonEFracMax=v;};
	
	inline double GetPhotonEFracMin(){return photonEFracMin;};
	inline void SetPhotonEFracMin(double v){photonEFracMin=v;};
};

#endif
