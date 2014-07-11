// ********************************************************************
// $Id: HRSVisAttribute.hh,v 1.0, 2010/12/26 HRS Exp $
// --------------------------------------------------------------
//

#ifndef HRSVisAttribute_h
#define HRSVisAttribute_h 1

class G4VisAttributes;

class HRSVisAttribute
{
public:
	HRSVisAttribute();
	virtual ~HRSVisAttribute();

private:
	void ConstructVisAttribute();

public:
	G4VisAttributes* HallVisAtt;
	G4VisAttributes* MagneticVisAtt;
	G4VisAttributes* WhiteVisAtt;
	G4VisAttributes* DarkRedVisAtt;
	G4VisAttributes* OrangeVisAtt;
	G4VisAttributes* GrayVisAtt;
	G4VisAttributes* YellowVisAtt;
	G4VisAttributes* PcbGreenVisAtt;
	G4VisAttributes* CuBrownVisAtt;
	G4VisAttributes* DarkBlueVisAtt;		//36;24;130
	G4VisAttributes* VioletVisAtt;			//238;130;238
	G4VisAttributes* LightYellowVisAtt;

	G4VisAttributes* BlackVisAtt; 
	G4VisAttributes* CableVisAtt;

	G4VisAttributes* RedVisAtt;				//255;204;50
	G4VisAttributes* YellowGreenVisAtt;		//153;204;50
	G4VisAttributes* LightPurpleVisAtt;		//155;48;255
	G4VisAttributes* PurpleVisAtt;			//128,0,128
	G4VisAttributes* DarkOrangeVisAtt;		//255,140,0
	G4VisAttributes* SkyBlueVisAtt;			//0;127;255
	G4VisAttributes* LightGreenVisAtt;
	G4VisAttributes* LightBlueVisAtt;

	G4VisAttributes* SteelVisAtt;
	G4VisAttributes* IronVisAtt;
	G4VisAttributes* SilverVisAtt;
	G4VisAttributes* LeadVisAtt;

};

#endif

