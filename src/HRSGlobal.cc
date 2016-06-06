// ********************************************************************
//
// $Id: HRSGlobal.cc,v 1.0, 2010/12/26   HRS Exp $
// GEANT4 tag $Name: geant4-09-04 $
//
//..............................................................................

#include "GlobalDebuger.hh"
#include "HRSGlobal.hh"
#include <iostream>
#include <fstream>

using namespace std;

//#define HRSCLOBAL_DEBUG 1

//////////////////////////////////////////////////////////////////////
//if unit==0, write into file only, otherwise write to both the file and the screen
void WriteLog(const char *str,const char *filename, int unit)
{
	FILE *gLog=fopen(filename,"a+");
	time_t TimeNow;
	time( &TimeNow );
	//if(unit!=0) printf("%s\t%s\n",ctime(&TimeNow),str);
	//fprintf(gLog,"%s\t%s\n",ctime(&TimeNow),str);
	fclose(gLog);
}


void InitLogStream(ofstream &logstream,const char *filename)
{ 
	if(!logstream.is_open()) logstream.open(filename,ios_base::app);
	time_t TimeNow;
	time( &TimeNow );
	//logstream<<ctime(&TimeNow);  //ctime already have a change line appended to the end
}


//////////////////////////////////////////////////////////////////////

bool CheckFile(const char *filename)
{
	string file(filename);
	return CheckFile(file);
}

bool CheckFile(string &val)
{
	ifstream pFile(val.c_str());
	if(pFile.is_open()) 
	{
		pFile.close();
		//cout<< ">>>>>>>>file " <<val.c_str()<<" is found!<<<<<<<<<<\n"; 
		return true;
	}
	else 
	{
		cout<< ">>>>>>>>file " <<val.c_str()<<" not exist!<<<<<<<<<<\n"; 
		return false;
	} 
}

//in c++, since string can resize automaticlly, no need to control the buffer
void ReplaceAll(string& strT,const char *oldsubstr, const char *newsubstr)
{

	int lenOld=strlen(oldsubstr);
	int lenNew=strlen(newsubstr);
#ifdef  HRSCLOBAL_DEBUG
	if( HRSCLOBAL_DEBUG >=1 )
	{
		cout<<"Before:"<<strT<<endl;
	}
#endif

	size_t pos;
	pos=strT.find(oldsubstr);
	while(pos!=string::npos) 
	{//found this substr already, replace it
		strT.replace(pos,lenOld,newsubstr);

#ifdef  HRSCLOBAL_DEBUG
		if( HRSCLOBAL_DEBUG >=2 )
		{ 
			printf("Replace one \"%s\" as \"%s\", the temp result is:\n%s\n",
				oldsubstr,newsubstr,strT.c_str());	  
		}
#endif
		//I should find the oldsubstr within the rest only
		//otherwise it will go into a daed loop if oldsubstr is part of newsubstr
		pos=strT.find(oldsubstr,pos+lenNew);
	}

#ifdef  HRSCLOBAL_DEBUG
	if( HRSCLOBAL_DEBUG >=1 )
	{
		cout<<" After:"<<strT<<endl;
	}
#endif 

}

void ReplaceAll(char *targetstr,const char *oldsubstr, const char *newsubstr)
{
	string strT=targetstr;
#ifdef  HRSCLOBAL_DEBUG
	if( HRSCLOBAL_DEBUG >=3 )
	{
		puts("====================================================================");
		puts("It is your responsibility to keep the size of the array pointed by\n\
			 targetstr shall be long enough to contain the string after replacing");
		puts("====================================================================");
	} 

#endif

	ReplaceAll(strT,oldsubstr,newsubstr);
	strcpy(targetstr,strT.c_str());

}
//////////////////////////////////////////////////////////////////////
/*
This subroutine will check the existance of the input string then decide the output file name.
One can specify these type of the output: root|bos|txt
for exsamle, CreateFileName("HRS_G4Sim_nt_10.root",outfilename,mytype) will have the
following outputs according to mytype: 
mytype=0 ==>HRS_G4Sim_nt_10_00.root 
mytype=1 ==>HRS_G4Sim_bos_10_00.A00 
mytype=2 ==>HRS_G4Sim_txt_10_00.txt 

arguments: instr: raw file name to be check, please use root file name;
outstr:  output file name;								      
type:    file type, 0=>root; 1=>bos; 2=>txt
*/

//c++ vertion
bool CreateFileName(char *instr,char *outstr,int outtype)
{
	const char *strI[]={"_nt_","_bos_","_txt_"};
	const char *strD[]={".root",".A00",".txt"};  


	size_t iPosDot=0;
	string strA,strB, In;
	In=instr;

	iPosDot=In.find_last_of('.');
	if(iPosDot==string::npos )
	{
		cout<<"there is no '.' in this input string "<<instr<<", \".root\" appended"<<endl;
		In.append(".root"); 
	}
	if(iPosDot+1==In.size() )
	{
		cout<<"there is nothing after '.' in this input string "<<instr<<", \"root\" appended"<<endl;
		In.append("root"); 
	}

	iPosDot=In.find_last_of('.');
	//cout<<"Last occurence of '.' found at "<<iPosDot<<endl;

	strA=In.substr(0,iPosDot);
	strB=In.substr(iPosDot,In.size()-iPosDot+1);

	int intype=0;
	for(int i=0;i<3;i++)
	{
		if(strB==strD[i]) {intype=i;break;}
	}

	FILE *f;
	int counter=0;
	bool bFoundFile=true;
	do{
		sprintf(outstr,"%s_%02d%s",strA.c_str(),counter,strB.c_str());
		//cout<<"checking file "<<outstr<<endl;
		if((f=fopen(outstr,"r"))==NULL) bFoundFile=false;
		else
		{
			bFoundFile=true;
			fclose(f);
			counter++;
		}

	}while(bFoundFile);

	//cout<<"intype="<<intype<< "  outtype="<<outtype<<endl;
	if(intype!=outtype && outtype>=0 && outtype<=2)
	{
		//replace "_nt_" with correct strI
		ReplaceAll(outstr,strI[intype],strI[outtype]);
		ReplaceAll(outstr,strD[intype],strD[outtype]);
	}


	char strLog[255];
	sprintf(strLog,"Create file %s\n",outstr);
	//puts(strLog);
	WriteLog(strLog);

	return true;

}
//////////////////////////////////////////////////////////////////////

