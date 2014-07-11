// ********************************************************************
// $Id: HRSPrimaryGeneratorMessenger.hh,v 1.0, 2010/12/26  HRS Exp $
// --------------------------------------------------------------
//
#ifndef HRSPrimaryGeneratorMessenger_h
#define HRSPrimaryGeneratorMessenger_h 1
//MaxPrimaryNum has been defined by HRSPrimaryGeneratorAction.hh
#include "HRSPrimaryGeneratorAction.hh"

class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithABool;
class G4UIcmdWith3Vector;
class G4UIcmdWith3VectorAndUnit;

class G4UIcmdWithoutParameter;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADouble;
class G4UIdirectory;

#include "G4UImessenger.hh"
#include "globals.hh"

class HRSPrimaryGeneratorMessenger: public G4UImessenger
{
public:
	HRSPrimaryGeneratorMessenger(HRSPrimaryGeneratorAction* mpga);
	~HRSPrimaryGeneratorMessenger();
public:
	void SetNewValue(G4UIcommand * command,G4String newValues);
	G4String GetCurrentValue(G4UIcommand * command);
private:
	HRSPrimaryGeneratorAction * target;
private: //commands
	G4UIdirectory*              mydetDir;
	G4UIdirectory*              gunDir[MaxPrimaryNum];
	//particle number cmds
	G4UIcmdWithAnInteger*       particleNumCmd;

	//////////////////////momentum cmd/////////////////////////////////
	//random
	G4UIcmdWithABool*           fastPsGenCmd;
	G4UIcmdWithADoubleAndUnit*  momentumLowCmd[MaxPrimaryNum];
	G4UIcmdWithADoubleAndUnit*  momentumHighCmd[MaxPrimaryNum];
	G4UIcmdWithADoubleAndUnit*  thetaLowCmd[MaxPrimaryNum];
	G4UIcmdWithADoubleAndUnit*  thetaHighCmd[MaxPrimaryNum];
	G4UIcmdWithADoubleAndUnit*  phiLowCmd[MaxPrimaryNum];
	G4UIcmdWithADoubleAndUnit*  phiHighCmd[MaxPrimaryNum];
	G4UIcmdWithADoubleAndUnit*  sigmaMomCmd[MaxPrimaryNum];
	G4UIcmdWithADoubleAndUnit*  sigmaThetaCmd[MaxPrimaryNum];
	G4UIcmdWithADoubleAndUnit*  sigmaPhiCmd[MaxPrimaryNum];
	//fixed momentum 3v
	G4UIcmdWithADoubleAndUnit*  momentumCmd[MaxPrimaryNum];
	G4UIcmdWithADoubleAndUnit*  thetaCmd[MaxPrimaryNum];
	G4UIcmdWithADoubleAndUnit*  phiCmd[MaxPrimaryNum];
	G4UIcmdWith3VectorAndUnit*  momentum3VCmd[MaxPrimaryNum];

	/////////////////////particle type cmds////////////////////////////
	G4UIcmdWithABool*           randomCmd;
	G4UIcmdWithAString*         particleNameCmd[MaxPrimaryNum];
	G4UIcmdWithAnInteger*       particlePDGCodeCmd[MaxPrimaryNum];

	////////////////position/vertex cmds//////////////////////////////
	//random
	G4UIcmdWithADoubleAndUnit*  gunZLowCmd;
	G4UIcmdWithADoubleAndUnit*  gunZHighCmd;
	G4UIcmdWithADoubleAndUnit*  gunZCmd;
	G4UIcmdWithADoubleAndUnit*  gunRLowCmd;
	G4UIcmdWithADoubleAndUnit*  gunRHighCmd;
	G4UIcmdWith3VectorAndUnit*  fixedPointBL3VCmd;
	G4UIcmdWith3Vector*         slopeBL3VCmd;
	G4UIcmdWithAnInteger*       rasterModeCmd;
	G4UIcmdWithAnInteger*       vertexModeCmd;

	//fixed
	G4UIcmdWith3VectorAndUnit*  position3VCmd;

	///////////////random engine cmds //////////////////////////////////
	//to specify the event generator engine types
	//available engine are: HRSElasEl,HRSElasNucleus,HRSQuasiElasEl,HRSQusiElasNucleon,BonusProton
	G4UIcmdWithAString*         engineTypeCmd[MaxPrimaryNum];
	//to randomize in TCS or HCS
	G4UIcmdWithADoubleAndUnit*	detectorAngleCmd[MaxPrimaryNum];
	G4UIcmdWithAnInteger*		randmizeInTCSCmd[MaxPrimaryNum];
	G4UIcmdWithADoubleAndUnit*  outPlaneAngleHighCmd[MaxPrimaryNum];
	G4UIcmdWithADoubleAndUnit*  outPlaneAngleLowCmd[MaxPrimaryNum];
	G4UIcmdWithADoubleAndUnit*  inPlaneAngleHighCmd[MaxPrimaryNum];
	G4UIcmdWithADoubleAndUnit*  inPlaneAngleLowCmd[MaxPrimaryNum];
	G4UIcmdWithAnInteger*		coupleToPrimaryCmd[MaxPrimaryNum];

	//to specify the beam energy and target for some engines
	G4UIcmdWithADoubleAndUnit*  beamEnergyCmd;
	G4UIcmdWithADoubleAndUnit*  tgMassCmd;
	G4UIcmdWithADouble*         tgAtomicNumberCmd;
	//use a 3 vector to specify the beam energy and target for some engines
	G4UIcmdWith3Vector*         beamNtarget3VCmd;
	G4UIcmdWithADoubleAndUnit*  leftHRSMomentumCmd;
	G4UIcmdWithADoubleAndUnit*  rightHRSMomentumCmd;

	G4UIcmdWithADoubleAndUnit*  HMSMomentumCmd;

	G4UIcmdWithADouble *PhotonEFracMinCmd;
	G4UIcmdWithADouble *PhotonEFracMaxCmd;
};

#endif


