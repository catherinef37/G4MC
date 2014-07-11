//HRS Exp, v 1.0, 2010/12/26 , By Jixie Zhang
// ********************************************************************
//
// $Id: HRSAnalysisManager.cc,v 1.0, 2010/12/26   HRS Exp $
// GEANT4 tag $Name: geant4-09-04 $
//
// **********************************************************************
//By Jixie: This class is currently not used. I plan to use this to replace
//or take care part of what HRSRootTree is doing.

#include "UsageManager.hh"
#include "GlobalDebuger.hh"
#include "HRSAnalysisManager.hh"


extern UsageManager* gConfig;

using namespace std;

HRSAnalysisManager* HRSAnalysisManager::instance = 0;

HRSAnalysisManager::HRSAnalysisManager()
{
	std::cout<<"HRSAnalysisManager() construction done!"<<std::endl;
}

HRSAnalysisManager* HRSAnalysisManager::GetAnalysisManager()
{
	if (!instance) instance = new HRSAnalysisManager();
	return instance;
}


HRSAnalysisManager::~HRSAnalysisManager()
{	
	if (instance)
	{
		delete instance;
		instance = 0;
	}
	std::cout<<"delete HRSAnalysisManager ... done!"<<std::endl;
}

