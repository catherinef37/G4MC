/////////////////////////////////////////////////////////////////////////////
// This is the implement of calss UsageManager
// Usage.ini will be used to build the structure which is stored at vector VBranch
// In order to acces the value quickly, I stored the values in the 3 maps according
// to its type: 
//	map<string,int> MapBr_i;
//	map<string,float> MapBr_f;
//	map<string,string> MapBr_s;
// For detail of how to built your own ini file, please refer to Usage.ini
// The class can also be used to read a regular parameter file, all parameter will 
// be stored as string in the map
//	map<string,string> MapParam_s;
/////////////////////////////////////////////////////////////////////////////
#include "GlobalDebuger.hh"
#include "UsageManager.hh"

//#define  USAGEMANAGER_DEBUG 1


UsageManager* UsageManager::fUsageManager=0;
UsageManager* UsageManager::GetUsageManager()
{ 
	if(!fUsageManager) 
	{
		std::cout<<"No instance of UsageManager found! "
			<<"Use UsageManager::UsageManager(filename,int argc, char **argv) to build one in your main()!"<<std::endl;
	}
	return fUsageManager; 
}


UsageManager::UsageManager(const char *inifile,int argc, char **argv,const char *logfilename)
{
#ifdef USAGEMANAGER_DEBUG
	if(Global_Debug_Level < (int)USAGEMANAGER_DEBUG) 
	{
		SetGlobalDebugLevel("UsageManager::UsageManager",(int)USAGEMANAGER_DEBUG);
	}
#endif

	if(fUsageManager)
	{ 
		cout<<"Error: UsageManager constructed twice.\n"
			<<"Use UsageManager::GetUsageManager() to get the static pointer UsageManager::fUsageManager\n";
		exit(-1);
	}
	fUsageManager=this;
	Executable=(argc>0 && argv)?argv[0]:"exe";
	LogFileName=logfilename;

	//store option structure and default values into the branch vector
	ReadIni(inifile);
	//process the comment line
	ProcessArgv(argc,argv);
}

UsageManager::~UsageManager()
{
	MapBr_i.clear();
	MapBr_f.clear();
	MapBr_s.clear();
	VBranch.clear();
	MapParam_s.clear();

	cout << "deleting UsageManager ......" << endl;
	fUsageManager = 0;
}

//////////////////////////////////////////////////////////////////////
//if unit==0, write into file only, otherwise write to both the file and the screen
void UsageManager::WriteLog(const char *str,const char *filename, int unit)
{
	FILE *gLog=fopen(filename,"a+");
	time_t TimeNow;
	time( &TimeNow );
	if(unit!=0) printf("%s\t%s\n",ctime(&TimeNow),str);
	fprintf(gLog,"%s\t%s\n",ctime(&TimeNow),str);
	fclose(gLog);
}

//if unit==0, write into file only, otherwise write to both the file and the screen
void UsageManager::WriteLog(const char *str)
{
	WriteLog(str,LogFileName.c_str(),1);
}
//////////////////////////////////////////////////////////////////////

//in c++, since string can resize automaticlly, no need to control the buffer
void UsageManager::ReplaceAll(string& strT,const char *oldsubstr, const char *newsubstr)
{
	size_t pos=strT.find(oldsubstr);
	if(pos==string::npos) return;

	int lenOld=int(strlen(oldsubstr));
	int lenNew=int(strlen(newsubstr));
#ifdef  USAGEMANAGER_DEBUG
	if( USAGEMANAGER_DEBUG >=4 )
	{
		cout<<"Before:"<<strT<<endl;
	}
#endif

	while(pos!=string::npos) 
	{//found this substr already, replace it
		strT.replace(pos,lenOld,newsubstr);

#ifdef  USAGEMANAGER_DEBUG
		if( USAGEMANAGER_DEBUG >=4 )
		{ 
			printf("Replace one \"%s\" as \"%s\", the temp result is:\n%s\n",
				oldsubstr,newsubstr,strT.c_str());	  
		}
#endif
		//I should find the oldsubstr within the rest only
		//otherwise it will go into a daed loop if oldsubstr is part of newsubstr
		pos=strT.find(oldsubstr,pos+lenNew);
	}

#ifdef  USAGEMANAGER_DEBUG
	if( USAGEMANAGER_DEBUG >=4 )
	{
		cout<<" After:"<<strT<<endl;
	}
#endif 

}

void UsageManager::ReplaceAll(char *targetstr,const char *oldsubstr, const char *newsubstr)
{
	string strT=targetstr;
#ifdef  USAGEMANAGER_DEBUG
	if( USAGEMANAGER_DEBUG >=3 )
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
//Create a unique file name using the given key stored in infilename
//if overwritten == true, will not care if the file exist or not, just generate key_#start#.root
//otherwise, it will check if  key_#start#.root exist or not, if exist, increase
//start by 1 and keep searching till find a name which is not exist.
//then return this name. for example: 
//CreateFileName("key.root",outfilename, 0, true) will return key_00.root
//This routine will find the last dot '.' and insert an index '_##' in fornt of the dot
void UsageManager::CreateFileName(const char *infilename,char outfilename[], int start, bool overwritten)
{
	FILE *f;
	char strA[1024],strB[1024];
	size_t counter=start, iPosDot=0;	  
	const char *pch=strrchr(infilename,'.');    //find the last char '.'
	iPosDot=size_t(pch-infilename+1);
	//printf ("Last occurence of '.' found at %d \n",pch-filename+1);
	strncpy(strA,infilename,iPosDot-1); //copy str before '.'
	strA[iPosDot-1]='\0';
	strcpy(strB,infilename+iPosDot-1);  //copy str after '.'

	bool bFoundFile=true;
	do{   
		sprintf(outfilename,"%s_%02d%s",strA,(int)counter,strB);
		if(overwritten) {bFoundFile=true; break;}
		else
		{
			if((f=fopen(outfilename,"r"))==NULL) bFoundFile=false;
			else
			{
				bFoundFile=true;
				fclose(f);
				counter++;
			}
		}

	}while(bFoundFile);

	char strLog[255];
	sprintf(strLog,"Create file %s\n",outfilename);
	WriteLog(strLog);
}

//////////////////////////////////////////////////////////////////////
//This routine was used in BoNuS Geant4 simulation. 
//It will check the existance of the input string then determine the output file name.
//Only the following types of the output is supported: root|bos|txt|mac
//If outtype is to 0|1|2|3, all key words, (i.e. '_nt') in the given string 'instr' will 
//be replaced according to the given outtype, otherwise no replacement happen.
//For exsamle, CreateFileName("Bonus_G4Sim_nt.root",outfilename,mytype) will have the
//following outputs according to mytype: 
//mytype=0 ==>Bonus_G4Sim_nt_00.root 
//mytype=1 ==>Bonus_G4Sim_bos_00.A00 
//mytype=2 ==>Bonus_G4Sim_txt_00.txt 
//mytype=3 ==>Bonus_G4Sim_txt_00.mac 
//If files of index 00 exist and overwritten==false, the index will increase 1 till that filename 
//is not exist.
//
//input: 
//instr:  the given file name to be check, please use root file name;
//outstr: output file name;								      
//outtype:output file type, 0=>root; 1=>bos; 2=>txt 3=>mac 	
void UsageManager::CreateFileName(char *instr,char *outstr,int outtype)
{
	const char *strI[]={"_nt","_bos","_txt","_mac"};
	const char *strD[]={".root",".A00",".txt",".mac"};  


	size_t iPosDot=0;
	string strA,strB, In;
	In=instr;

	iPosDot=In.find_last_of('.');
	if(iPosDot==string::npos )
	{
#ifdef  USAGEMANAGER_DEBUG
	if( USAGEMANAGER_DEBUG >=3 )
	{
		cout<<"there is no '.' in this input string "<<instr<<", \".root\" appended"<<endl;
	}
#endif
		In.append(".root"); 
	}
	if(iPosDot+1==In.size() )
	{
#ifdef  USAGEMANAGER_DEBUG
	if( USAGEMANAGER_DEBUG >=3 )
	{
		cout<<"there is nothing after '.' in this input string "<<instr<<", \"root\" appended"<<endl;;
	}
#endif
		In.append("root"); 
	}

	iPosDot=In.find_last_of('.');
	//cout<<"Last occurence of '.' found at "<<iPosDot<<endl;

	strA=In.substr(0,iPosDot);
	strB=In.substr(iPosDot,In.size()-iPosDot+1);

	int intype=0;
	for(int i=0;i<4;i++)
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

	//cout<<"intype="<<intype<<"  outtype="<<outtype<<endl;
	if(intype!=outtype && outtype>=0 && outtype<=3)
	{
		//replace "_nt_" with correct strI
		ReplaceAll(outstr,strI[intype],strI[outtype]);
		ReplaceAll(outstr,strD[intype],strD[outtype]);
	}


	char strLog[255];
	sprintf(strLog,"Create file %s\n",outstr);
	//puts(strLog);
	WriteLog(strLog);
}

//////////////////////////////////////////////////////////////////////
//remove space or tab \t at both end of the string
void  UsageManager::Trim(std::string &astr)
{
	if(astr.length()<1) return;
	size_t i=astr.length()-1;
	while(astr[i]==' ' || astr[i]=='\t' || astr[i]=='\r' )
	{
		astr.erase(astr.length()-1,1);
		i=astr.length()-1;
	}
	i=0;
	while(astr[i]==' ' || astr[i]=='\t' || astr[i]=='\r' )
	{
		astr.erase(0,1);
	}
}

void UsageManager::Trim(char *astr)
{
	if(strlen(astr)<1) return;

	size_t i=strlen(astr)-1;
	bool search=true;
	do{
		if(astr[i]==' ' || astr[i]=='\t' || astr[i]=='\r' ) 
		{
			astr[i]='\0';
			i--;
		}
		else search=false;
	}while(search);

	i=0;
	search=true;
	do{
		if(astr[i]==' ' || astr[i]=='\t' || astr[i]=='\r' ) 
		{
			i++;
		}
		else search=false;
	}while(search);

	if(i>0)
	{
		std::string cppstr(astr+i);
		strcpy(astr,cppstr.c_str());
	}
}

std::string UsageManager::TrimString(std::string &astr)
{
	Trim(astr);
	return astr;
}

std::string UsageManager::TrimString(char *astr)
{
	std::string cppstr(astr);
	Trim(cppstr);
	strcpy(astr,cppstr.c_str());
	return cppstr;
}


//Keep reading lines from stream ifs till find one valid content, treat ererything after the first
//'#' as comment, ignoring space ' ' and tab key '\t' at both end of the content
std::string UsageManager::ReadOneValidLine(ifstream &ifs)
{
	char line[1024];
	string wholeline, value;
	while (ifs.good() && value.size()<1) 
	{
		ifs.getline(line,1024);
		wholeline=line;
		size_t found=wholeline.find("#");
		if (found!=string::npos)
		{
			//chucate the line from this position
			value=wholeline.substr(0,found);
		}
		else
		{
			value=wholeline;
		}
	}
	Trim(value);
	return value;
}

//split a sentence (string) into words using the following characters " \t,;"
//allow no more than 10 words, return the number of words
size_t UsageManager::SplitList(const string str, string list[],size_t allownum)
{
	size_t i=0;
	char *cstr=0, *p=0;	
	cstr = new char [str.size()+1];
	strcpy (cstr, str.c_str());
	// cstr now contains a c-string copy of str
	const char *spliter=" \t,;";
	p=strtok(cstr,spliter);
	while (p!=NULL)
	{
		list[i]=p;
		Trim(list[i]);
		i++;
		//cout<<"list["<<i-1<<"]="<<p<<endl;
		if(i>=allownum) break;
		p=strtok(NULL,spliter);
	}

	delete cstr;
	return i;
}

//print one argument branch, if quit
void UsageManager::PrintBranch(branch_t &br,bool quit)
{
	if(br.argc_require<0) return;
				
	cout.flush(); 
	if(quit) cout<<"The detail usage for option '-"<<br.key1<<"' or '-"<<br.key2<<"' is: "<<endl;
	cout<<endl<<"  -"<<br.key1<<" or -"<<br.key2;
	for(size_t j=0;j<br.argc_given;j++) 
	{
		if(j<size_t(br.argc_require))
		{
			cout<<" <"<<br.namelist[j]<<"("<<br.valuelist[j]<<")>";
		}
		else
		{
			cout<<" ["<<br.namelist[j]<<"("<<br.valuelist[j]<<")]";
		}
		if(j!=br.argc_given-1  && !((j+1)%4) ) cout<<endl<<"\t  ";
	}
	//replace $t with \t, $n with \n 
	cout<<":"<<endl;
	string astr(br.description);
	//check if there is a "description" in the string, if it exist, no need to print it 
	size_t found=astr.find("description");
	if(found==string::npos || found>15) cout<<"\t description: ";
	if(astr.find("\\t")!=string::npos) ReplaceAll(astr,"\\t","\t");
	if(astr.find("\\n")!=string::npos) ReplaceAll(astr,"\\n","\n");	
	if(astr.find("$t")!=string::npos) ReplaceAll(astr,"$t","\t");	
	if(astr.find("$n")!=string::npos) ReplaceAll(astr,"$n","\n");	
	cout<<astr<<endl;
	if(quit) exit(-999);
}

//Record all the arguments into the map
void UsageManager::BuildMap()
{
	for(size_t i=0;i<VBranch.size();i++)
	{		
		size_t argc=(VBranch[i].argc_given==0)?1:VBranch[i].argc_given;
		for(size_t j=0;j<argc;j++) 
		{		
			if(VBranch[i].typelist[j].compare("int")==0 || VBranch[i].typelist[j].compare("size_t")==0) 
				MapBr_i[VBranch[i].namelist[j]]=atoi(VBranch[i].valuelist[j].c_str());
			else if(VBranch[i].typelist[j].compare("float")==0 || VBranch[i].typelist[j].compare("double")==0) 
				MapBr_f[VBranch[i].namelist[j]]=float(atof(VBranch[i].valuelist[j].c_str()));
			else if(VBranch[i].typelist[j].compare("string")==0 || VBranch[i].typelist[j].compare("char")==0) 
				MapBr_s[VBranch[i].namelist[j]]=VBranch[i].valuelist[j];
		}
	}
}
	

//print the argument map, which is the name and value for all command line arguments 
void UsageManager::PrintMap(ostream &pOut)
{
	if(MapBr_i.size() + MapBr_f.size() + MapBr_s.size() <= 0) return;
	pOut<<"*************************MAP**VALUE**START***************************"<<endl;
	pOut<<"PrintMap(): The content of the arguments in the map are:\n";

	map<string,int>::iterator it_i;
	for ( it_i=MapBr_i.begin() ; it_i != MapBr_i.end(); it_i++ )
	{
		pOut<<"MapBr_i["<<setw(22)<<(*it_i).first<<"] = "<<(*it_i).second<<endl;
	}

	map<string,float>::iterator it_f;
	for ( it_f=MapBr_f.begin() ; it_f != MapBr_f.end(); it_f++ )
	{
		pOut<<"MapBr_f["<<setw(22)<<(*it_f).first<<"] = "<<(*it_f).second<<endl;
	}

	map<string,string>::iterator it_s;
	for ( it_s=MapBr_s.begin() ; it_s != MapBr_s.end(); it_s++ )
	{
		pOut<<"MapBr_s["<<setw(22)<<(*it_s).first<<"] = \""<<(*it_s).second<<"\""<<endl;
	}
	pOut<<"*************************MAP**VALUE****END***************************"<<endl;

}
	
//print the parameter map, which is the name and values written in the parameter files
void UsageManager::PrintParamMap(ostream &pOut)
{
	if(MapParam_s.size() <= 0) return;
	pOut<<"*************************MAP**VALUE**START***************************"<<endl;
	pOut<<"PrintParamMap(): The content of the parameters in the map are:\n";
	map<string,string>::iterator it_s;
	for ( it_s=MapParam_s.begin() ; it_s != MapParam_s.end(); it_s++ )
	{
		pOut<<"MapParam_s["<<setw(22)<<(*it_s).first<<"] = \""<<(*it_s).second<<"\""<<endl;
	}
	pOut<<"*************************MAP**VALUE****END***************************"<<endl;
}


//print all the arguments
void UsageManager::PrintOpt(ostream &pOut)
{
	pOut<<"****************************VALUE**START*****************************"<<endl;
	pOut<<"PrintOpt(): The content of the arguments in the vector are:\n";
	for(size_t i=0;i<VBranch.size();i++)
	{		
		size_t argc=(VBranch[i].argc_given<=0)?1:VBranch[i].argc_given;
		pOut<<"-"<<VBranch[i].key1<<" or -"<<VBranch[i].key2<<":"<<endl;
		for(size_t j=0;j<argc;j++) 
		{		
			pOut<<setw(12)<<VBranch[i].namelist[j]<<"="<<VBranch[i].valuelist[j]<<"  ";
			if(j==argc-1 || !((j+1)%4) ) pOut<<endl;
		}
	}	
	pOut<<"****************************VALUE****END*****************************"<<endl;
}

void UsageManager::PrintUsage(bool bExitAfterPrint)
{		
	cout<<"*********************************MENU**START******************************"<<endl;
	cout<<"This is the usage generated from class UsageManager. \n"
		<<"Usage: "<<Executable<<" [option1 argument_list] [option2 argument_list] [...]\n\n"
		<< " To get the detail usage of each option, just type 'help' after it \n"
		<<"  Available options and default values(value inside the parenthesis):"<<endl;

	for(size_t i=0;i<VBranch.size();i++)
	{
		PrintBranch(VBranch[i]);
	}
	cout<<"*********************************MENU****END******************************"<<endl;
	if(bExitAfterPrint) exit(-999);
}

void UsageManager::ReadIni(const char *inifile)
{
	ifstream ifs;
	ifs.open(inifile,ios_base::in);
	
	char strLog[255];
	if(!ifs.is_open()) 
	{
		cout<<"*******************************************************"<<endl;
		sprintf(strLog,"File \"%s\" is not found. Don't know how to precess arguments.",inifile);
		WriteLog(strLog);
		cout<<"*******************************************************"<<endl;
	}
	else
	{
		sprintf(strLog,"UsageManager::ReadIni() is trying to read file \"%s\" ",inifile);
		WriteLog(strLog);
	}

	string value;
	while (ifs.good())     // loop while extraction from file is possible
	{
		value=ReadOneValidLine(ifs);
		if(value.size()<1) continue;
		if(value.compare(0,8,"[branch]")==0 || value.compare(0,8,"[BRANCH]")==0)
		{
			branch_t aBranch;
			aBranch.id=atoi(ReadOneValidLine(ifs).c_str());
			aBranch.key1=ReadOneValidLine(ifs);
			aBranch.key2=ReadOneValidLine(ifs);
			aBranch.argc_require=atoi(ReadOneValidLine(ifs).c_str());
			int tmp=atoi(ReadOneValidLine(ifs).c_str());
			if(tmp<0) 
			{ 
				cout<<"***Error: According to file "<<inifile<<", option -"<<aBranch.key1
					<<": argc_given should be non-negative! exit..."<<endl; 
				exit(-99);
			}
			aBranch.argc_given=size_t(tmp);
			SplitList(ReadOneValidLine(ifs),aBranch.typelist,aBranch.argc_given);
			SplitList(ReadOneValidLine(ifs),aBranch.namelist,aBranch.argc_given);
			SplitList(ReadOneValidLine(ifs),aBranch.valuelist,aBranch.argc_given);
			//description part
			aBranch.description=ReadOneValidLine(ifs);
			char b=aBranch.description[aBranch.description.length()-2];
			char c=aBranch.description[aBranch.description.length()-1];
#ifdef USAGEMANAGER_DEBUG
			c=aBranch.description[aBranch.description.length()-1];
			//this part is used to debug the difference between between windows and linux
			// if string="and exit.", The following lines will print: 
			//linux:   char b=. (int)b=46 char c=^M (int)c=13   change line \=\  (int)\=92
			//windows: char b=t (int)b=116 char c=. (int)c=46   change line \=\  (int)\=92
			// if string=".root ", The following lines will print: 
			//linux:   char b=  (int)b=32 char c=^M (int)c=13   change line \=\  (int)\=92
			//windows: char b=t (int)b=116 char c=  (int)c=32   change line \=\  (int)\=92
			// if "events. \", The following lines will print: 
			//linux:   char b=\ (int)b=92 char c=^M (int)c=13   change line \=\  (int)\=92
			//windows: char b=  (int)b=32 char c=\ (int)c=92   change line \=\  (int)\=92

			//conclusion:  The linux treat the last
			if(Global_Debug_Level>=5)
			{
				cout<<"char b="<<b<<" (int)b="<<(int)b <<" char c="<<c<<" (int)c="<<(int)c 
					<<"   change line \\="<<'\\'<<"  (int)\\="<<int('\\')<<endl;
			}
#endif
			//take care of the linux difference
			//It looks like have something relative to clrf and rf end line
			while( (int)c==92 || ((int)c==13 && (int)b==92) )
			{//read the next line and append it to the description
#ifdef USAGEMANAGER_DEBUG
				if(Global_Debug_Level>=4)
				{
					cout<<"A change line sysble '\\' found, will be replace by '\\n':"<<endl;
					cout<<aBranch.description<<endl;
				}
#endif		
				size_t found=aBranch.description.find_last_of('\\');
				aBranch.description=aBranch.description.substr(0,found);
				aBranch.description.append("\n");
				c='0';

				//need to check if the next line is the key word [branch]
				//if it is, then roll the position back
				//this will happen if user put a '\' at the end of the description but missing the continuous line
				int pos_marker1=ifs.tellg();
				string nextline=ReadOneValidLine(ifs);
				if( nextline.find("[branch]")!=string::npos ||
					nextline.find("[BRANCH]")!=string::npos )
				{
					int pos_marker2=ifs.tellg();
					for(int i=0;i<pos_marker2-pos_marker1-1;i++) ifs.unget();

#ifdef USAGEMANAGER_DEBUG
					if(Global_Debug_Level>=5)
					{
						cout<<" pos_marker1="<<pos_marker1<<" pos_marker2="<<pos_marker2<<endl;	
						cout<<"pos_marker2:"<<nextline<<endl;	
						char *buffer = new char [pos_marker2-pos_marker1+1];
						ifs.read (buffer,pos_marker2-pos_marker1-1);
						buffer[pos_marker2-pos_marker1-1]='\0';
						cout<<"pos_marker1:"<<buffer<<endl;	
						cout.write(buffer,pos_marker2-pos_marker1-1);
						delete buffer;
					}
#endif
				}
				else
				{
					aBranch.description.append(nextline);
					c=aBranch.description[aBranch.description.length()-1];
				}
			}

			VBranch.push_back(aBranch);
#ifdef USAGEMANAGER_DEBUG
			if(Global_Debug_Level>=3)
			{
				PrintBranch(aBranch);
			}
#endif
		}
	}

	ifs.close();           // close file
	return;

}

void UsageManager::ProcessArgv(int argc, char **argv)
{  
#ifdef USAGEMANAGER_DEBUG
	if(Global_Debug_Level>=3)
	{
		cout<<"argc="<<argc<<endl;
		for(int i=0;i<argc;i++) cout<<argv[i]<<"  ";
		cout<<endl;
	}
	if(Global_Debug_Level>=2)
	{
		cout<<"Before processing the input argument:"<<endl;
		PrintOpt();
	}
#endif

	int marker;
	branch *br=0;
	char *argptr;
	for (int i=1; i<argc; i++)
	{
		argptr = argv[i]; 
		if (*argptr == '-')
		{
			argptr++;

			if(strcmp(argptr,"help")==0  || strcmp(argptr,"h")==0 ||
				strcmp(argptr,"-help")==0 || strcmp(argptr,"-h")==0 ) 
			{
				PrintUsage(true);
			}
			//////////////////////////////////////////////////////////////////////
			//reset the branch
			br=0;
			marker=i;
			size_t b,j;
			for(b=0;b<VBranch.size();b++)
			{//locate the branch
				if(VBranch[b].key1.compare(argptr)==0  || VBranch[b].key2.compare(argptr)==0 )
				{
					br=&VBranch[b]; 
					break;
				}
			}
			if(!br) 
			{
				if(strcmp(argptr,"help")==0  || strcmp(argptr,"HELP")==0) PrintUsage(true);
				cout<<endl<<"Warning: !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ..."<<endl;
				cout<<"Warning: unrecognized option: \""<<argv[i]<<"\". Ignored ..."<<endl;
				cout<<"Warning: !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ..."<<endl<<endl;
				continue;
			}

			//not allow '-' exist, except it is a number
			if(argc<=i+br->argc_require) 
			{			  
				cout<<"***Error: in option: \""<<argv[i]<<"\", argument missing (require "
					<<br->argc_require<<" arguments). Exit ..."<<endl;
				PrintBranch(*br,true);
			}
			//required arguments and optional arguments
			for(j=0;j<br->argc_given;j++)
			{	
				//if this is the last input argument, break		 
				if(i>=argc-1) break;
				//check if this argument is "help" or not, if it is, print the help for this option
				if(strcmp(argv[i+1],"help")==0) 
				{
					PrintBranch(*br,true);
				}

				if(argv[i+1][0]!='-' || (argv[i+1][0]=='-' && isdigit(argv[i+1][1])) )
				{
					br->valuelist[j]=string(argv[++i]);	
				}
				else
				{
					if(j<size_t(br->argc_require))  
					{
						cout<<"***Error: option '-"<<br->key1<<"' or '-"<<br->key2<<"' requires "<<br->argc_require
							<<"argument(s), but only "<<j<<" argument(s) found! exit..."<<endl;
						PrintBranch(*br,true);
					}
					else  break;
				}		
			}

			//I add this routine here to take care of reserve options
			ProcessReservedOption(br);

			//reset the argc_require
			//argc_given is used in BuildMap(), should not be updated
			br->argc_require=(int)j;
			//
#ifdef USAGEMANAGER_DEBUG
			if(Global_Debug_Level>=3)
			{
				for(int jj=marker;jj<=i;jj++) cout<<"  argv["<<jj<<"]="<<argv[jj];
				cout<<endl;
				PrintBranch(*br);
				STOP4DEBUG;
			}
#endif

			//////////////////////////////////////////////////////////////////////
		}//end if (*argptr == '-')
		else
		{				
			if(strcmp(argptr,"help")==0  || strcmp(argptr,"HELP")==0) PrintUsage(true);
			cout<<endl<<"Warning: !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ..."<<endl;
			cout<<"Warning: unrecognized option: \""<<argv[i]<<"\". Ignored ..."<<endl;
			cout<<"Warning: !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ..."<<endl<<endl;
			continue;
		}

	}//end for

#ifdef USAGEMANAGER_DEBUG
	if(Global_Debug_Level>=2)
	{
		cout<<"After processing the input argument:"<<endl;
		PrintOpt();
		PrintMap();
	}
	if(Global_Debug_Level>=1)	  PrintMap();
#endif

	BuildMap();


	////////////////////////////////////////////////////////////////////
	//write all arguments into the log file
	////////////////////////////////////////////////////////////////////
	ofstream fLog;
	if(!fLog.is_open()) fLog.open(LogFileName.c_str(),ios_base::app);
	time_t TimeNow;
	time( &TimeNow );
	fLog<<ctime(&TimeNow); 
	fLog<<"The input command is:"<<endl;
	for(int i=0;i<argc;i++)	fLog<<argv[i]<<" ";
	fLog<<endl;
	this->PrintOpt(fLog);
	fLog.close();
	////////////////////////////////////////////////////////////////////

}

//////////////////////////////////////////////////////////////////////
void UsageManager::ProcessReservedOption(branch *br)
{
	//take care of a "switch" type of option, which has argc_require=argc_given=0 
	//and it is data type should be treated as size_t or int
	//if this option is invoked, just set the default value to its inverse
	if( br->argc_require==0 && br->argc_given==0 && 
		(br->typelist[0].compare("int")==0 || br->typelist[0].compare("size_t")==0) )
	{
		br->valuelist[0]=(br->valuelist[0].compare("0")==0)?"1":"0";
	}

	if( (br->key1.compare("INPUTFILE")==0 || br->key1.compare("inputfile")==0) )
	{
		for(size_t j=0;j<br->argc_given;j++)
		{	
			//take care of the reserved "-inputfile" or "-INPUTFILE" option, 
			//if this option is invoked, all giving config file will be read in. 
			// All parameters wirtten in these configuration files cab be achieved by 
			//calling UsageManager::Getparameter(string name, T &para);
			if(CheckPath(br->valuelist[j])) ReadFile(br->valuelist[0].c_str());
		}
	}
}
//////////////////////////////////////////////////////////////////////

bool UsageManager::GetArgument(const std::string name,int &argument)
{
	bool found=false;
	map<std::string,int>::iterator it_i;

	it_i=MapBr_i.find(name);
	if(it_i != MapBr_i.end()) 
	{
		argument=(*it_i).second; 
		found=true;
	}
	else 
	{
		std::cout<<"***Error: GetArgument(Name="<<name<<",int): wrong argument name!!! "<<std::endl;
		exit(-999);
	}
	return found;
}

bool UsageManager::GetArgument(const std::string name,size_t &argument)
{
	int arg;
	bool found=GetArgument(name,arg);
	argument=(size_t)arg;
	return found;
}

bool UsageManager::GetArgument(const std::string name,float &argument)
{
	bool found=false;
	map<std::string,float>::iterator it_f;

	it_f=MapBr_f.find(name);
	if(it_f != MapBr_f.end()) 
	{
		argument=(*it_f).second; 
		found=true;
	}
	else 
	{
		std::cout<<"***Error: GetArgument(Name="<<name<<",float): wrong argument name!!! "<<std::endl;
		exit(-999);
	}
	return found;
}

bool UsageManager::GetArgument(const std::string name,double &argument)
{
	float arg;
	bool found=GetArgument(name,arg);
	argument=(double)arg;
	return found;
}

bool UsageManager::GetArgument(const std::string name,std::string &argument)
{
	bool found=false;
	map<std::string,std::string>::iterator it_s;

	it_s=MapBr_s.find(name);
	if(it_s != MapBr_s.end()) 
	{
		argument=(*it_s).second; 
		found=true;
	}
	else 
	{
		std::cout<<"***Error: GetArgument(Name="<<name<<",string): wrong argument name!!! "<<std::endl;
		exit(-999);
	}
	return found;
}

bool UsageManager::GetArgument(const std::string name,char* argument)
{
	string arg;
	bool found=GetArgument(name,arg);
	if(found)
	{
		size_t length = strlen(argument);
		if(length<arg.length()+1) argument=new char [arg.length()+1];
		strcpy(argument,arg.c_str());
	}
	return found;
}

//using the branch vector to find the variables
std::string UsageManager::GetArgument(std::string name)
{	
	bool found=false;
	std::string argument("");
	for(size_t i=0;i<VBranch.size();i++)
	{		
		size_t argc=(VBranch[i].argc_given==0)?1:VBranch[i].argc_given;
		for(size_t j=0;j<argc;j++) 
		{		
			if(VBranch[i].namelist[j].compare(name)==0) 
			{
				argument=VBranch[i].valuelist[j];
				found=true;	
				break;
			}		
		}
		if(found) break;
	}
	if(!found) 
	{
		std::cout<<"***Error: GetArgument(Name="<<name<<"): wrong argument name!!! "<<std::endl;
		exit(-999);
	}
	return argument;
}

std::string UsageManager::GetArgument(const char *name)
{ 	
	std::string argname(name);
	return GetArgument(argname);
}


bool UsageManager::CheckPath(const char *filename)
{
	string file(filename);
	return CheckPath(file);
}

bool UsageManager::CheckPath(string &val)
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

//c++ version of function to read a parameter file
//Input: filename c string
//Output: return false if file can not open succesfully, otherwise return true 
//        all parameters written in the given file will be stored and get be
//        achieved by calling bool GetParameter<typename>(char *name, T &value);
//        Please NOTE that the typename could be size_t,int,float,double,bool,char*,string only
//feature: 
//1)The parameter files must in a format as "parameter=value", allow space ' ' beside '='. 
//2)The value is ended by comma ',', semi-comma ';' or space ' ', care should be taken for string parameters
//2)Use '#' to comment out lines or part of lines. Everything after # will be ignored.
//3)A line without '=' will not be understood. All non-understood lines will be printed out and ignored.
//4)The maximum length of each line is 1024.

bool UsageManager::ReadFile(const char *filename)
{
	ifstream ini(filename,ios_base::in);
	if(!ini.good())
	{
		cout<<"***Error: Can not open parameter file \""<<filename<<"\"!"<<endl;
		return false;
	}
	else  cout<<"Reading parameter file \""<<filename<<"\""<<endl;

	int nline=0;
	size_t found=0;
	static const int MaxLineSize=1024;
	char line[MaxLineSize];
	std::string strLine;
	std::string strName;
	std::string strValue;

	while(!ini.eof())
	{
		ini.getline(line,MaxLineSize);
		nline++;
		//Search forward. for the '#' to skip the afterwards
		strLine=line;
		found=strLine.find('#');
		if(found!=string::npos)	strLine=strLine.substr(0,found);
	
		//trim the string
		Trim(strLine);
		if(strLine.length()==0) continue; //this is a comment line whose first char is '#'

		//try to find '=', if not find, ignore this line		
		found=strLine.find('=');
		if (found==string::npos || strLine.length()<3) 
		{
			cout<<"Warning! line "<<nline<<" is not recognized: \""<<line<<"\" \n";
			continue;
		}
		strName=strLine.substr(0,found);
		Trim(strName);	

		strValue=strLine.substr(found+1);
		found=strValue.find(';');
		if (found!=string::npos) strValue=strValue.substr(0,found);
		Trim(strValue);	
		cout<<setw(25)<<strName<<" = "<<strValue<<endl;
		//store the value into the map
		MapParam_s[strName]=strValue;

	}
	ini.close();

	////////////////////////////////////////////////////////////////////
	//write all parameters into log file
	////////////////////////////////////////////////////////////////////
	ofstream fLog;
	if(!fLog.is_open()) fLog.open(LogFileName.c_str(),ios_base::app);
	time_t TimeNow;
	time( &TimeNow );
	fLog<<ctime(&TimeNow); 
	this->PrintParamMap(fLog);
	fLog.close();
	////////////////////////////////////////////////////////////////////

	return true;

}

bool UsageManager::ReadFile_c(const char *filename)
{
	FILE *ini;
	if((ini=fopen(filename,"r"))==NULL)
	{
		printf("***Error: Can not open parameter file \"%s\"!\n",filename);
		return false;
	}
	else  printf("Reading parameter file \"%s\"\n",filename);

	char ch[]="=;\n";
	int nline=0;
	static const int MaxLineSize=1024;
	char *name,*value,line[MaxLineSize];
	char *pDest;
	std::string strName,strValue;

	while(!feof(ini))
	{
		//empty the buffer
		sprintf(line,"%s","");
		fgets(line,MaxLineSize,ini);
		nline++;

		//Search forward. for the '#' to skip the fterwards
		pDest = strchr( line, '#' );
		if(pDest!=NULL) *pDest='\0'; 
		//trim the string
		Trim(line);
		if(strlen(line)==0) continue; //this is a comment line

		//try to find '=', if not find, ignore this line		
		pDest = strchr( line, '=' );
		if(pDest==NULL || strlen(line)<3) 
		{
			printf("Warning! line %d is not recognized: \"%s\" \n",nline,line);
			continue;
		}

		name=strtok(line,ch);
		value=strtok(0,ch);

		//store the value into the map
		Trim(name);
		Trim(value);
		strName=name;
		strValue=value;
		printf("%15s = %s\n",strName.c_str(),strValue.c_str());
		MapParam_s[strName]=strValue;
		
	}
	fclose(ini);

	////////////////////////////////////////////////////////////////////
	//write all parameters into log file
	////////////////////////////////////////////////////////////////////
	ofstream fLog;
	if(!fLog.is_open()) fLog.open(LogFileName.c_str(),ios_base::app);
	time_t TimeNow;
	time( &TimeNow );
	fLog<<ctime(&TimeNow); 
	this->PrintParamMap(fLog);
	fLog.close();
	////////////////////////////////////////////////////////////////////

	return true;
}


//using the branch vector to find the variables
std::string UsageManager::GetParameter(std::string name)
{	
	bool found=false;
	std::string parameter("");
	map<string,string>::iterator it;

	it=MapParam_s.find(name);
	if(it != MapParam_s.end()) 
	{
		parameter=(*it).second; 
		found=true;
	}

	if(!found) 
	{
		std::cout<<"***Error: GetParameter(Name="<<name<<"): wrong parameter name!!!"<<std::endl; 
		std::cout<<"***Error: One of your ini files is out of date or written with a wrong parameter name!"<<std::endl;
		exit(-999);
	}
	return parameter;
}

std::string UsageManager::GetParameter(const char *name)
{ 	
	std::string para(name);
	return GetParameter(para);
}

/*
bool UsageManager::GetParameter(const std::string name,size_t &parameter)
{
std::string para=GetParameter(name);
bool found=(para.length()==0)?false:true;
parameter=(size_t) atoi(para.c_str());
return found;
}
bool UsageManager::GetParameter(const std::string name,int    &parameter)
{
std::string para=GetParameter(name);
bool found=(para.length()==0)?false:true;
parameter=(int) atoi(para.c_str());
return found;
}
bool UsageManager::GetParameter(const std::string name,float  &parameter)
{
std::string para=GetParameter(name);
bool found=(para.length()==0)?false:true;
parameter=(float) atof(para.c_str());
return found;
}
bool UsageManager::GetParameter(const std::string name,double &parameter)
{
std::string para=GetParameter(name);
bool found=(para.length()==0)?false:true;
parameter=(double) atof(para.c_str());
return found;
}

*/
bool UsageManager::GetParameter(const std::string name,std::string &parameter)
{
	parameter=GetParameter(name);
	bool found=(parameter.length()==0)?false:true;
	return found;
}
bool UsageManager::GetParameter(const std::string name,char *parameter)
{
	std::string para=GetParameter(name);
	bool found=(para.length()==0)?false:true;
	sprintf(parameter,"%s",para.c_str());
	return found;
}

//////////////////////////////////////////////////////////////////

bool UsageManager::SetParameter(std::string name,std::string value)
{
	bool found=false;
	map<string,string>::iterator it;

	it=MapParam_s.find(name);
	if(it != MapParam_s.end()) 
	{
		it->second=value; 
		found=true;
		std::cout<<"###SetParameter(name="<<name<<", value="<<value
			<<")  called!###"<<std::endl;
	}

	if(!found) 
	{
		std::cout<<"***Error: SetParameter(name="<<name<<", value="<<value
			<<"): wrong parameter name!!! exit ... "<<std::endl;
		exit(-999);
	}
	return found;
}

bool UsageManager::SetParameter(std::string name,const char* value)
{
	std::string theStr(value);
	return SetParameter(name,theStr);
}

bool UsageManager::SetParameter(std::string name,char value[])
{
	std::string theStr(value);
	return SetParameter(name,theStr);
}
//////////////////////////////////////////////////////////////////

bool UsageManager::SetArgument(std::string name,std::string value)
{
	bool found=false;
	map<string,string>::iterator it;

	it=MapBr_s.find(name);
	if(it != MapBr_s.end()) 
	{
		it->second=value; 
		found=true;
		std::cout<<"###SetArgument(name="<<name<<", value="<<value
			<<")  called!###"<<std::endl;
	}

	if(!found) 
	{
		std::cout<<"***Error: SetArgument(name="<<name<<", value="<<value
			<<"): wrong parameter name!!! exit ... "<<std::endl;
		exit(-999);
	}
	return found;
}

bool UsageManager::SetArgument(std::string name,const char* value)
{
	std::string theStr(value);
	return SetArgument(name,theStr);
}

bool UsageManager::SetArgument(std::string name,char value[])
{
	std::string theStr(value);
	return SetArgument(name,theStr);
}

