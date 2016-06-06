// ********************************************************************
// $Id: HRSMaterial.hh,v 1.0, 2010/12/26 HRS Exp $
// --------------------------------------------------------------
//

#ifndef HRSMaterial_h
#define HRSMaterial_h 1

#include <map>
#include "globals.hh"	//for units and g4io 

#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4MaterialTable.hh"

class HRSMaterial
{
public:
	HRSMaterial();
	virtual ~HRSMaterial();
	G4Material* FindMaterial(G4String name);

	//Static method which returns the singleton pointer of this class.
	static HRSMaterial* GetHRSMaterialManager();

private:
	static HRSMaterial* fHRSMaterialManager;


private:
	void ConstructMaterials();
	void DestroyMaterials();
	void BuildMaterialMap();

public:
	std::map<G4String, G4Material*>  mMaterialMap;

	G4Material* galaxy;
	G4Material* vacuum;
	G4Material* air;
	G4Material* argonGas;
	G4Material* heliumGas;
	G4Material* heliumGas4k;
	G4Material* aluminum;
	G4Material* calcium;
	G4Material* iron;
	G4Material* carbon;
	G4Material* diamond;

	G4Material* tantalum;
	G4Material* tungsten;
	G4Material* lead;
	G4Material* lead208;
	G4Material* solidNH3;
	G4Material* liquidHe;
	G4Material* liquidHe3;
	G4Material* liquidH2;
	G4Material* liquidD2;

	G4Material* scintillator;	//EJ-204 Plastic Scintillator
	G4Material* mylar;			//C10_H8
	G4Material* kapton;			//C22_H10_N2_O5
	G4Material* copper;
	G4Material* siliconsteel;
	G4Material* stainlesssteel;
	G4Material* stainlesssteel304;
	G4Material* NH3He;			//55% SolidNH3 + 45% LiquidHe

	G4Material* ultem;			//C37_H24_O6_N2
	G4Material* SiO2;
	G4Material* epoxy;			//C11_H12_O3
	G4Material* G10FR4;			//60% SiO2 and 40% epoxy
	G4Material* CH2;			//H(CH2)nH
	G4Material* ethane;			//C2H6
	G4Material* EthaneHe;		//ethane and He4 gas mixture
	G4Material* pcbNchip;
	G4Material* cable;
	
	G4Material* beryllium;
	G4Material* BeO;

	G4Material* boron;
	G4Material* plastic;
	G4Material* boratedpoly05;
	G4Material* boratedpoly30;

	//PCTFE (Polychlorotrifluoroethylene)  [CF2-CFCl]n
	//synthetic resin formed by the polymerization of chlorotrifluoroethylene
	//PCTFE is a homopolymer of chlorotrifluoroethylene. Features of PCTFE 
	//include high compressive strength and low deformation under load. In 
	//particular, its cold flow characteristic is lower than other fluoropolymers 
	//and it does not deform under load at room temperature. PCTFE has extremely 
	//low gas permeability and essentially does not absorb moisture.
	//PCTFE resin is produced by Daikin, under the trade name Neoflon.
	G4Material* PCTFE;   
	
	G4Material* MuMetal;
	G4Material* absorber;  //something with huge density to absorb particles

private:
	//the following will be read from configuration file detector.ini
	double mNH3WeightRatio;
	double mSolidNH3D,mLiquidHeD;
};
typedef HRSMaterial HRSMaterialManager;

#endif

