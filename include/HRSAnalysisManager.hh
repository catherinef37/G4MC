//HRS Exp, v 1.0, 2010/12/26 , By Jixie Zhang
// ********************************************************************
//
// $Id: HRSAnalysisManager.hh,v 1.0, 2010/12/26   HRS Exp $
// GEANT4 tag $Name: geant4-09-04 $
//

#ifndef HRSAnalysisManager_h
#define HRSAnalysisManager_h 1

#include "globals.hh"
#include <vector>
#include "G4ThreeVector.hh"

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"


class HRSAnalysisManager {
public:
	HRSAnalysisManager();
	virtual ~HRSAnalysisManager();
	static HRSAnalysisManager* GetAnalysisManager();

private:
	static HRSAnalysisManager* instance;

};

#endif

