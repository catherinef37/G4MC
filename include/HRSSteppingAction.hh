// ********************************************************************
//
// $Id: HRSSteppingAction.hh,v2.0 2007/01/01 HRS Exp $
//
//..............................................................................

#ifndef HRSSteppingAction_h
#define HRSSteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "globals.hh"
#include <ostream>
#include <fstream>
#include <vector>
using namespace std;

//..............................................................................
class HRSSteppingActionMessenger;
class HRSSteppingAction : public G4UserSteppingAction
{
public:
	HRSSteppingAction();
	~HRSSteppingAction();

	void UserSteppingAction(const G4Step*);

	inline void SetVerbose(G4int val) { verboseLevel = val; }
	inline G4int GetVerbose() const { return verboseLevel; }

	inline void EmptyPrintList() { vPrintList.clear(); }
	inline void Add2PrintList(G4String val) { vPrintList.push_back(val); }

private:
	void InitOutTxt();
	void DoPrint(const G4Step*);
	void PrintHead(const G4Step* theStep, ostream& pOut=G4cout);
	void PrintStep(const G4Step*, ostream& pOut=G4cout);
	void FillRootArray(const G4Step* theStep);

private: 
	HRSSteppingActionMessenger* messenger;

	vector <G4String> vPrintList;
	bool CreateTxt;
	bool CreateRootNt;

	G4int verboseLevel;
	ofstream OutTxt;
	bool PrintHeadFlag;

};

//..............................................................................

#endif
