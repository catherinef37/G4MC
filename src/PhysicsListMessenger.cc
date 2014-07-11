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
//
// $Id: PhysicsListMessenger.cc, v 1.0, 2010/12/26  HRS Exp $
// GEANT4 tag $Name: geant4-09-04-beta-01 $
//
//
/////////////////////////////////////////////////////////////////////////


#include "PhysicsListMessenger.hh"

#include "PhysicsList.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UImanager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsListMessenger::PhysicsListMessenger(PhysicsList* pPhys)
:pPhysicsList(pPhys)
{   
	gammaCutCmd = new G4UIcmdWithADoubleAndUnit("/hadronPhys/gammaCut",this);  
	gammaCutCmd->SetGuidance("Set gamma cut.");
	gammaCutCmd->SetParameterName("GammaCut",false);
	gammaCutCmd->SetUnitCategory("Length");
	gammaCutCmd->SetRange("GammaCut>=0.0");
	gammaCutCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

	electCutCmd = new G4UIcmdWithADoubleAndUnit("/hadronPhys/electronCut",this);  
	electCutCmd->SetGuidance("Set electron cut.");
	electCutCmd->SetParameterName("ElCut",false);
	electCutCmd->SetUnitCategory("Length");
	electCutCmd->SetRange("ElCut>=0.0");
	electCutCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

	posCutCmd = new G4UIcmdWithADoubleAndUnit("/hadronPhys/positronCut",this);
	posCutCmd->SetGuidance("Set positron cut.");
	posCutCmd->SetParameterName("PositronCut",false);
	posCutCmd->SetUnitCategory("Length");
	posCutCmd->SetRange("PositronCut>=0.0");
	posCutCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

	pCutCmd = new G4UIcmdWithADoubleAndUnit("/hadronPhys/protonCut",this);
	pCutCmd->SetGuidance("Set proton cut.");
	pCutCmd->SetParameterName("ProtonCut",false);
	pCutCmd->SetUnitCategory("Length");
	pCutCmd->SetRange("ProtonCut>=0.0");
	pCutCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

	allCutCmd = new G4UIcmdWithADoubleAndUnit("/hadronPhys/allCut",this);
	allCutCmd->SetGuidance("Set cut for all.");
	allCutCmd->SetParameterName("AllCut",false);
	allCutCmd->SetUnitCategory("Length");
	allCutCmd->SetRange("AllCut>=0.0");
	allCutCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

	addPhysCmd = new G4UIcmdWithAString("/hadronPhys/add",this);
	addPhysCmd->SetGuidance("Add one modula physics list. You can find all the available physics model");
	addPhysCmd->SetGuidance(" with this command: /hadronPhys/list. Here are some of them:");
	addPhysCmd->SetGuidance(" FTFP_BERT FTFP_BERT_EMV FTFP_BERT_EMX FTF_BIC ");
	addPhysCmd->SetGuidance(" LHEP LHEP_EMV QBBC QGS_BIC QGSC_BERT ");
	addPhysCmd->SetGuidance(" QGSP QGSP_BERT QGSP_BERT_EMV QGSP_BIC_EMY ");
	addPhysCmd->SetGuidance(" QGSP_BERT_EMX QGSP_BERT_HP QGSP_BIC QGSP_BIC_HP");
	addPhysCmd->SetGuidance(" For details about these models, please refer to ");
	addPhysCmd->SetGuidance(" http://geant4.cern.ch/support/proc_mod_catalog/physics_lists/referencePL.shtml");
	addPhysCmd->SetParameterName("PhysList",false);
	addPhysCmd->AvailableForStates(G4State_PreInit);

	listCmd = new G4UIcmdWithoutParameter("/hadronPhys/list",this);
	listCmd->SetGuidance("List all the available physics models. ");
	listCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsListMessenger::~PhysicsListMessenger()
{
	delete gammaCutCmd;
	delete electCutCmd;
	delete posCutCmd;
	delete pCutCmd;
	delete allCutCmd;
	delete addPhysCmd;
	delete listCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsListMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
	G4UImanager* UI = G4UImanager::GetUIpointer();
	if( command == gammaCutCmd ) 
	{
		if(pPhysicsList) 
		{
			pPhysicsList->SetCutForGamma(gammaCutCmd->GetNewDoubleValue(newValue));
		} 
		else 
		{
			UI->ApplyCommand("/run/setCutForAGivenParticle gamma " + newValue);
		}

	} 
	else if( command == electCutCmd ) 
	{
		if(pPhysicsList) 
		{
			pPhysicsList->SetCutForElectron(electCutCmd->GetNewDoubleValue(newValue));
		} 
		else 
		{
			UI->ApplyCommand("/run/setCutForAGivenParticle e- " + newValue);
		}

	} 
	else if( command == posCutCmd ) 
	{
		if(pPhysicsList) 
		{
			pPhysicsList->SetCutForPositron(posCutCmd->GetNewDoubleValue(newValue));
		} 
		else 
		{
			UI->ApplyCommand("/run/setCutForAGivenParticle e+ " + newValue);
		}

	} 
	else if( command == pCutCmd ) 
	{
		if(pPhysicsList) 
		{
			pPhysicsList->SetCutForProton(pCutCmd->GetNewDoubleValue(newValue));
		} 
		else 
		{
			UI->ApplyCommand("/run/setCutForAGivenParticle proton " + newValue);
		}

	} 
	else if( command == allCutCmd ) 
	{
		if(pPhysicsList) 
		{
			G4double cut = allCutCmd->GetNewDoubleValue(newValue);
			pPhysicsList->SetCutForGamma(cut);
			pPhysicsList->SetCutForElectron(cut);
			pPhysicsList->SetCutForPositron(cut);
			pPhysicsList->SetCutForProton(cut);
		} 
		else 
		{
			UI->ApplyCommand("/run/setCut " + newValue);
		}

	} 
	else if( command == addPhysCmd ) 
	{
		if(pPhysicsList) 
		{
			pPhysicsList->AddPhysicsList(newValue);
		} 
		else 
		{
			G4cout << "### PhysicsListMessenger WARNING: "
				<< " /hadronPhys/Physics UI command is not available "
				<< "for reference Physics List" << G4endl;
		}

	} else if( command == listCmd ) 
	{
		if(pPhysicsList) 
		{
			pPhysicsList->List();
		} 
		else 
		{ 
			G4cout << "### PhysicsListMessenger WARNING: "
				<< " /hadronPhys/ListPhysics UI command is not available "
				<< "for reference Physics List" << G4endl;
		}
	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
