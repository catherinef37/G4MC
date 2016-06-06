// ********************************************************************
//
// $Id: HRSPrimaryRootEvent.cc,v 1.0, 2011/7/14  HRS Exp $
// GEANT4 tag $Name: geant4-09-04 $
//
// --------------------------------------------------------------
//
#include "HRSPrimaryRootEvent.hh"


#include "GlobalDebuger.hh"
//#define HRSPRIMARYROOTEVENT_DEBUG 3

HRSPrimaryRootEvent::HRSPrimaryRootEvent()
{
#ifdef HRSPRIMARYROOTEVENT_DEBUG
    if(Global_Debug_Level < HRSPRIMARYROOTEVENT_DEBUG) 
    {
        Global_Debug_Level = HRSPRIMARYROOTEVENT_DEBUG;
        cout<<"HRSPrimaryRootEvent(): Set Global_Debug_Level to "<<Global_Debug_Level<<endl;
    }
#endif

	for(int i=0;i<MaxTrackNum;i++)
	{
		mNtuple[i]=0;mCurrentEntry[i]=0;
	}

	//cout<<"HRSPrimaryRootEvent() construction done!"<<endl;
}

HRSPrimaryRootEvent::~HRSPrimaryRootEvent()
{
	for(int i=0;i<MaxTrackNum;i++)
	{
		if(mNtuple[i]) delete mNtuple[i];
	}
	//cout<<"delete HRSPrimaryRootEvent ... done!"<<endl;
}


void HRSPrimaryRootEvent::LoadNtuple(int trackid, const char* filename, const char* treename,
									 int trignum, int skipnum)
{ 
    mNtuple[trackid]=new HRSNtReader(filename,treename);
    if (mNtuple[trackid]->mEntries < 1) 
    {
        cout<<"\nWarning!!! Tree "<<treename<<" in ntuple "<<filename<<"is empty!"<<endl;
        return;
    }
	if(skipnum>0) mCurrentEntry[trackid]=skipnum-1;
	if(trignum>0) mNtuple[trackid]->mEntries=trignum;
}


bool HRSPrimaryRootEvent::EventCut(int trackid)
{	
	if( mNtuple[trackid]->P0<0 || mNtuple[trackid]->PdgId==-1 || 
		fabs(mNtuple[trackid]->Theta0+10.0)<0.1) return false;
	if(trackid==0)
	{
		if(mNtuple[trackid]->TrackClass >= 0 ) return true;
	}
	else
	{
		if(mNtuple[trackid]->TrackClass >= 0 ) return true;
	}
	return false;
}

//read vertex and initial momentum from the tree
//return   -1: reach end of file
//       <=-2: error
//          0: get one event correctly
int HRSPrimaryRootEvent::GetParticle(int trackid,int &pdgid,TVector3& v3x,TVector3& v3p)
{
   
	if(trackid<0 || trackid>=MaxTrackNum) return -2;
	if(!mNtuple[trackid]) return -3;

	//make it short
	HRSNtReader* pNt=mNtuple[trackid];

	//return -1 if already reach the end of file
	if(pNt->mEntries <= mCurrentEntry[trackid]) return -1;

	bool bFoundOneEvent=false;
	while(!bFoundOneEvent)
	{
		pNt->GetEntry(mCurrentEntry[trackid]++);
		if(!EventCut(trackid))
		{			
			if(pNt->mEntries <= mCurrentEntry[trackid]) return -1;
			else continue;
		}
		else
		{
			bFoundOneEvent=true;
		}
	}
	
	//get one valid events, load the vertex and momentum
	pdgid=pNt->PdgId;
	v3x.SetXYZ(pNt->X0,pNt->Y0,pNt->Z0);
	v3p.SetMagThetaPhi(pNt->P0,pNt->Theta0,pNt->Phi0); 

#ifdef HRSPRIMARYROOTEVENT_DEBUG
    if(Global_Debug_Level>=2)
    {
		cout<<"HRSPrimaryRootEvent::GetParticle(): pdgid="<<pdgid
			<<"  v3x="<<v3x<<"  v3p="<<v3p<<endl;
    }
#endif

    return 0;
}





