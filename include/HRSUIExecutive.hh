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
// $Id: HRSUIExecutive.hh, v 1.0, 2010/12/26 HRS Exp $
// GEANT4 tag $Name: geant4-09-03 $
//
// ====================================================================
//   HRSUIExecutive.hh
//
//   By Jixie: This class helps automatic instantiation of user session
//   according to your environment variable like G4UI_USE_XXX. 
//   For QT and XM, it will try to load gui.mac 
//
//   Usage in main():
//
//   ...
//   #include "HRSUIExecutive.hh"
//
//   int main(int argc, char** argv)
//   {
//     ...
//     HRSUIExecutive* myapp = new HRSUIExecutive(argc, argv, usegui);
//	   G4UIsession* session = myapp->GetSession();
//     if (session->IsGUI())
//       // Do any extra for a GUI session
//
//     myapp-> SessionStart();
//     ...
//     delete myapp;
//     ...
//
// ====================================================================
#ifndef HRS_UI_EXECUTIVE_H
#define HRS_UI_EXECUTIVE_H 1

#include "G4VUIshell.hh"
#include "G4UIsession.hh"
#include "G4UImanager.hh"

#include "G4UIterminal.hh"
#include "G4UIcsh.hh"

#if defined(G4UI_USE_QT)
#include "G4UIQt.hh"
#include "G4Qt.hh"
#endif

#if defined(G4UI_USE_XM)
#include "G4UIXm.hh"
#endif

#if defined(G4UI_USE_WIN32)
#include "G4UIWin32.hh"
#endif

#if defined(G4UI_USE_TCSH)
#include "G4UItcsh.hh"
#endif

#define DISCARD_PARAMETER(p) (void)p

class HRSUIExecutive {
private:
	G4UIsession* session;
	G4VUIshell* shell;
	G4bool isGUI;

public:

	HRSUIExecutive(G4int argc, char** argv, int usegui=0)
	{
		session=0; shell=0; isGUI=false;
		//G4UImanager* pUImanager =0;
		//pUImanager = G4UImanager::GetUIpointer();

		if(usegui)
		{
#if defined(G4UI_USE_QT)
			session = new G4UIQt(argc, argv);
			isGUI = true;
			//pUImanager->ApplyCommand("/control/execute gui.mac");
#elif defined(G4UI_USE_XM)
			session = new G4UIXm(argc, argv);
			isGUI = true;
			//pUImanager->ApplyCommand("/control/execute gui.mac");			
#elif defined(G4UI_USE_WIN32)
		session = new G4UIWin32();
#elif defined(G4UI_USE_TCSH)
		shell = new G4UItcsh;
		session = new G4UIterminal(shell);
#else
		shell = new G4UIcsh;
		session = new G4UIterminal(shell);
#endif
		}
		else
		{			
			DISCARD_PARAMETER(argc);
			DISCARD_PARAMETER(argv);

#if defined(G4UI_USE_WIN32)
		session = new G4UIWin32();
#elif defined(G4UI_USE_TCSH)
		shell = new G4UItcsh;
		session = new G4UIterminal(shell);
#else
		shell = new G4UIcsh;
		session = new G4UIterminal(shell);
#endif
		}
	};

	~HRSUIExecutive(){ if(shell) delete shell;};

	G4bool IsGUI(){ return isGUI;};

	G4UIsession* GetSession() const{ return session; };

	void SetPrompt(const G4String& prompt){ if(shell) shell-> SetPrompt(prompt); };

	void SetLsColor(TermColorIndex dirColor, TermColorIndex cmdColor)
	{
		if(shell) shell-> SetLsColor(dirColor, cmdColor);
	};

	void SessionStart() { session->SessionStart(); delete session;};

};

#endif
