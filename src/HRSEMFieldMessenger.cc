//HRS Exp, v 1.0, 2010/12/26 , By Jixie Zhang
// ********************************************************************
//
// $Id: HRSElectricFieldMessenger.cc,v 1.0, 2010/12/26 HRS Exp $
// GEANT4 tag $Name: geant4-09-04 $
//
// 

#include "HRSEMFieldMessenger.hh"
#include "G4UnitsTable.hh"
//////////////////////////////////////////////////////////////////////////////
//
HRSEMFieldMessenger::HRSEMFieldMessenger(HRSEMField* pEMfield)
:target(pEMfield)
{ 
	new G4UnitDefinition("kilovolt/cm", "kV/cm", "ElectricField", kilovolt/cm);
	ErDCCmd = new G4UIcmdWithADoubleAndUnit("/field/ErDC",this);  
	ErDCCmd->SetGuidance("Define Electric Field in Drift Region");
	ErDCCmd->SetGuidance("Electric field will be in -R direction.");
	ErDCCmd->SetParameterName("ErDC",false,false);
	ErDCCmd->SetDefaultUnit("kV/cm");  //kilovolt/cm
	ErDCCmd->AvailableForStates(G4State_Idle); 

	ErInnerCmd = new G4UIcmdWithADoubleAndUnit("/field/ErInner",this);  
	ErInnerCmd->SetGuidance("Define Electric Field in Inner Drift Region");
	ErInnerCmd->SetGuidance("Electric field will be in +R direction.");
	ErInnerCmd->SetParameterName("ErInner",false,false);
	ErInnerCmd->SetDefaultUnit("kV/cm");  
	ErInnerCmd->AvailableForStates(G4State_Idle);  

	EField3VCmd = new G4UIcmdWith3VectorAndUnit("/field/setEField3V",this);  
	EField3VCmd->SetGuidance("Set to use an uniform electric field");
	EField3VCmd->SetParameterName("Ex","Ey","Ez",false);
	EField3VCmd->SetDefaultUnit("kilovolt/cm");  
	EField3VCmd->AvailableForStates(G4State_Idle);  

	BField3VCmd = new G4UIcmdWith3VectorAndUnit("/field/setBField3V",this);  
	BField3VCmd->SetGuidance("Set to use an uniform magnetic field");
	BField3VCmd->SetParameterName("Bx","By","Bz",false);
	BField3VCmd->SetDefaultUnit("tesla");  
	BField3VCmd->SetUnitCandidates("tesla kilogauss gause");
	BField3VCmd->AvailableForStates(G4State_Idle);  
}


///////////////////////////////////////////////////////////////////////////////

HRSEMFieldMessenger::~HRSEMFieldMessenger()
{
	delete ErDCCmd;
	delete ErInnerCmd;
	delete EField3VCmd;
	delete BField3VCmd;
}

////////////////////////////////////////////////////////////////////////////
//
void HRSEMFieldMessenger::SetNewValue( G4UIcommand* command, G4String newValue)
{ 
	if( command == ErDCCmd )
	{ 
		target->SetErDC(ErDCCmd->GetNewDoubleValue(newValue));
	}  
	if( command == ErInnerCmd )
	{ 
		target->SetErInner(ErInnerCmd->GetNewDoubleValue(newValue));
	}

	if( command==EField3VCmd )
	{ 
		target->SetEField3V(EField3VCmd->GetNew3VectorValue(newValue)); 
	}

	if( command==BField3VCmd )
	{ 
		target->SetBField3V(BField3VCmd->GetNew3VectorValue(newValue)); 
	}
}


/////////////////////////////////////////////////////////////////////////
//
G4String HRSEMFieldMessenger::GetCurrentValue( G4UIcommand* command)
{ 
	G4String cv;
	if( command == ErDCCmd )
	{ 
		cv=command->ConvertToString(target->GetErDC());
	}  
	if( command == ErInnerCmd )
	{ 
		cv=command->ConvertToString(target->GetErInner());
	}

	if( command==EField3VCmd )
	{  
		cv=command->ConvertToString(target->GetEField3V());
	}

	if( command==BField3VCmd )
	{ 
		cv=command->ConvertToString(target->GetBField3V()); 
	}
	return cv;
}

