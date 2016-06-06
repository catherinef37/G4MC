/*////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#By Jixie Zhang, 11/11/2008  
#This is the config file for the program menu. This file can be used as an parameter configuration.
#The following feature is included:
#(1) The length of each line should be less than 1024.
#(2) Use ($n ot \n) and ($t or \t) to change line and indent 
#(3) Everthing after the first # in the line will be ignored,therefore you can place your comment any where.
#(4) Use [branch] or [BRANCH] to start one option. In each option, you can config to have either 
#    required arguments or optional arguments or both. Allow no more than 10 arguments in each option. Argument
#    types, names and default values should be provided for all allowed arguments.
#(5) There is no limit on the number of branchs. The order of the branchs is not important.   
#(6) "-h" or "-help" has been built. An usage (help menu) is built automaticly using all the description. 
#    You can still add you help menu "--h or --help" like the branch id 0. Please note that once the usage 
#    menu printed, the program will be terminated. 
#(7) Please NOTE that the branch and argument should be input like the following structures and
#    in the exact order. But the description part can have more than one lines, use character \ at 
#    the end of the line to tell that the next line is the continuous part of this one.
#
#typedef struct branch{
#	int     id;	            //menu 
#	string	key1;
#	string	key2;
#	int     argc_require;   //if less than 0, this branch is not an option, it is just a list of config variabls
#                           //which will not be printed out in help menu 
#	size_t  argc_given;		//allow no more than 10 arguments in each branch, 
#	                        //this number will set the number of arguments of this branch
#	string  typelist[10];	//support [size_t,] int, float, [double,] string only
#	string	namelist[10];	//name list
#	string	valuelist[10];	//value list
#	string  description;	//description of this option, will be printed only if argc_require NOT less than 0
#	                        //Please use "$t" or "\t" for tab key, "\" "\n" or "$n" to change line
#	                        //if key word "description" is missing at the first 15 characters, "$t description:" will be printed 
#}branch_t;
#
#    In the command line, your input can be: "exe -key1 argc_required argc_optional" or "exe -key2 argc_required argc_optional"
#    you don't need to add a minus sign - in front of the keys.
#(8) In the above list, you can use any of these 3 characters( tab key '\t' or space ' ' or comma ',') to split the members in the list
#(9) Please note every branch is optional. you don't have to specify all or any of them. However, if you specify
#    one and the required arguments is not provided, the help menu will be printed and the program will exit. If you
#    provide too many arguments, the program take what it need and ignore the rest and meanwhile the help menu will
#    be print but program will still run.
#(10)To create config variabls which do not need to be input in the command line, just set argc_require less than 0 !
#    Its description will not be printed in the help menu. But user can use these vaiables in the code. 
#(11)There may be some special options which do not followed by any arguments. I called this type of option as a "switch" type.
#    For switch type, the argc_require and argc_allow must be zero and the data type must be size_t or int. Once a switch is invoked,
#    its value will be set to its inverse. Please note that if you invoke this option twice then it is equal to not be invoked at all. 
#(12)There is a reserve option, whose key1 and key2 are not allowed to change. When this arument is invoked the giving config files
#    will be read in. All parameters wirtten in these configuration files cab be achieved by calling 
#    UsageManager::Getparameter(string name, T &para);
#Enjoy! 
#########################################################################
*/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifndef _USAGET_H_
#define _USAGET_H_

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <time.h>
#include <string>
#include <vector>
#include <map>
#include <typeinfo>

#include <ctype.h>
using namespace std;


typedef struct branch{
	int     id;	            //menu 
	string	key1;
	string	key2;
	int     argc_require;   //if less than 0, this branch is not an option, it is just a list of config variabls
                            //which will not be printed out in help menu 
	size_t  argc_given;		//allow no more than 10 arguments in each branch, 
	                        //this number will set the number of arguments of this branch
	string  typelist[10];	//support [size_t,] int, float, [double,] string only
	string	namelist[10];	//name list
	string	valuelist[10];	//value list
	string  description;	//description of this option, will be printed only if argc_require NOT less than 0
	                        //Please use "\t" or "$t" for tab key, "\", "\n" or "$n" to change line
	                        //if key word "description" is missing at the first 15 characters, "$t description:" will be printed 
}branch_t;


static std::string LogFileName;

class UsageManager
{
public: 
	//Static method which returns the singleton pointer of UsageManager or its derived class.
	static UsageManager* GetUsageManager();
private:
	static UsageManager* fUsageManager;

public: 
	//thie construction of this class
	UsageManager(const char *inifile="Usage.ini", int argc=0, char **argv=0,const char *logfilename="Log.txt~");
	virtual ~UsageManager();

	//process the command line
	void ProcessArgv(int argc, char **argv);

	//print the usage|menual, will exit the program if bExitAfterPrint==true 
	void PrintUsage(bool bExitAfterPrint=false);

	//print all the arguments
	void PrintOpt(ostream &pOut=cout);

	//print the argument map, which is the name and value for all command line arguments 
	void PrintMap(ostream &pOut=cout);

	//print the parameter map, which is the name and values written in the parameter files
	void PrintParamMap(ostream &pOut=cout);

	//The major GET function for the arguments which were written in the usage ini file 
	//Get the argument by name, which can be found in the branch, return an empty string if not found
	std::string GetArgument(std::string name);
	std::string GetArgument(const char *name);
	//Get argument which has the same type as the given one, return true if found otherwise return false 
	bool GetArgument(const std::string name,size_t &argument);
	bool GetArgument(const std::string name,int    &argument);
	bool GetArgument(const std::string name,float  &argument);
	bool GetArgument(const std::string name,double &argument);
	bool GetArgument(const std::string name,std::string &argument);
	bool GetArgument(const std::string name,char *argument);

	//changing the parameter in the map, both argument map and parameter map
	bool SetParameter(std::string name,std::string value);
	bool SetParameter(std::string name,const char * value);
	bool SetParameter(std::string name,char value[]);
	template <typename T>
	bool SetParameter(const std::string name,T value)
	{	
	  //bool found=false;
		char tmpstr[255];
		//std::cout<<"\n$$$$$$$$$this type is "<<typeid(value).name()<<"$$$$$$$$$"<<std::endl;
		if(typeid(value)==typeid(float)    || typeid(value)==typeid(double)) 
			sprintf(tmpstr,"%f",value);
		else if(typeid(value)==typeid(int) || typeid(value)==typeid(unsigned int) ||
			typeid(value)==typeid(long)    || typeid(value)==typeid(unsigned long) ||
			typeid(value)==typeid(short)   || typeid(value)==typeid(unsigned short) ) 
			sprintf(tmpstr,"%d",(int)value);
		else if(typeid(value)==typeid(bool)) 
			sprintf(tmpstr,"%d",(value)?1:0);
		//this type doesnot support very well
		//else if(typeid(value)==typeid(char*)) 
		//	sprintf(tmpstr,"%s",value);
		else 
		{
			std::cout<<"Unknown input type \""<<typeid(value).name()<<std::endl;
			//found=false;
		}
		std::string valuestr(tmpstr);
		return SetParameter(name,valuestr);
	}

	//////////////////////////////////////////////////////////////
	bool SetArgument(std::string name,std::string value);
	bool SetArgument(std::string name,const char* value);
	bool SetArgument(std::string name,char value[]);
	template <typename T>
	bool SetArgument(const std::string name,T value)
	{	
	  bool found=false;

		if(typeid(value)==typeid(float)    || typeid(value)==typeid(double)) 
		{
			float theValue=(float)value;
			map<string,float>::iterator it=MapBr_f.find(name);
			if(it != MapBr_f.end()) 
			{
				it->second=theValue; 
				found=true;
			}
		}
		else if(typeid(value)==typeid(int) || typeid(value)==typeid(unsigned int) ||
			typeid(value)==typeid(long)    || typeid(value)==typeid(unsigned long) ||
			typeid(value)==typeid(short)   || typeid(value)==typeid(unsigned short) ||
			typeid(value)==typeid(bool) ) 	
		{
			int theValue;
			if(typeid(value)==typeid(bool)) theValue=(value)?1:0;
			else theValue=(int)value;
			map<string,int>::iterator it=MapBr_i.find(name);
			if(it != MapBr_i.end()) 
			{
				it->second=theValue; 
				found=true;
			}
		}
		//this type doesnot support very well
		//else if(typeid(value)==typeid(char*))
		//{
		//	std::string theValue(value);
		//	map<string,string>::iterator it=MapBr_s.find(name);
		//	if(it != MapBr_s.end()) 
		//	{
		//		it->second=theValue; 
		//		found=true;
		//	}
		//}
		else 
		{
			std::cout<<"Unknown input type \""<<typeid(value).name()<<std::endl;
			found=false;
		}

		if(found)
		{
			std::cout<<"###SetArgument(name="<<name<<", value="<<value<<") called! ###"<<std::endl;
		}
		else
		{
			std::cout<<"***Error: SetArgument(name="<<name<<", value="<<value
			<<"): wrong argument name!!! exit ... "<<std::endl;
			exit(-999);
		}
		return found;
	}


	//The majoy GET function for the parameters which were written in the parameter file
	//Get the parameter by name, which can be found in the parameter file, return an empty string if not found
	std::string GetParameter(std::string name);
	std::string GetParameter(const char *name);
	
	//Get argument which has the same type as the given one, return true if found otherwise return false 
	bool GetParameter(const std::string name,std::string &argument);
	bool GetParameter(const std::string name,char *argument);
	//the following five have been included in the template
	//bool GetParameter(const std::string name,size_t &argument);
	//bool GetParameter(const std::string name,int    &argument);
	//bool GetParameter(const std::string name,float  &argument);
	//bool GetParameter(const std::string name,double &argument);
	//bool GetParameter(const std::string name,bool  &argument);
	
	//Get argument which has the same type as the given one, return true if found otherwise return false 
	//Basiclly all regular 'value' type, i.e signed and unsigned short, int, long, float and double
	//will be included in this template
	template <typename T>
	bool GetParameter(const std::string name,T &argument)
	{
		bool found=true;
		std::string arg=GetParameter(name);
		if(arg.length()<1) return false;
 
		if(typeid(argument)==typeid(float)    || typeid(argument)==typeid(double)) 
			argument=(T)atof(arg.c_str());
		else if(typeid(argument)==typeid(int) || typeid(argument)==typeid(unsigned int) ||
			typeid(argument)==typeid(long)    || typeid(argument)==typeid(unsigned long) ||
			typeid(argument)==typeid(short)   || typeid(argument)==typeid(unsigned short) ) 
			argument=(T)atoi(arg.c_str());
		else if(typeid(argument)==typeid(bool)) 
			argument=(arg.compare("true")==0  || arg.compare("TRUE")==0) ? true : false;
		else 
		{
			std::cout<<"Unknown input type \""<<typeid(argument).name()<<std::endl;
			found=false;
		}
		return found;
	}


	//c++ version of function to read a parameter file
	//Input: filename c string
	//Output: return false if file can not open succesfully, otherwise return true 
	//        all parameters written in the given file will be stored and get be
	//        achieved by calling bool GetParameter<typename>(char *name, T &value);
	//        Please NOTE that the typename could be size_t,int,float,double,bool,char*,string only
	//feature: 
	//1)The parameter files must in a format as "parameter=value", allow space ' ' beside '='. 
	//2)The value is ended by semi-comma ';'. Care should be taken for string parameters
	//2)Use '#' to comment out lines or part of lines. Everything after # will be ignored.
	//3)A line without '=' will not be understood. All non-understood lines will be printed out and ignored.
	//4)The maximum length of each line is 1024.
	bool ReadFile_c(const char *filename);

	//c version of function to read a parameter file
	//Input: filename c string
	//Output: return false if file can not open succesfully, otherwise return true 
	//        all parameters written in the given file will be stored and get be
	//        achieved by calling bool GetParameter<typename>(char *name, T &value);
	//        Please NOTE that the typename could be size_t,int,float,double,bool,char*,string only
	//feature: 
	//1)The parameter files must in a format as "parameter=value", allow space ' ' beside '='. 
	//2)The value is ended by semi-comma ';'. Care should be taken for string parameters
	//2)Use '#' to comment out lines or part of lines. Everything after # will be ignored.
	//3)A line without '=' will not be understood. All non-understood lines will be printed out and ignored.
	//4)The maximum length of each line is 1024.
	bool ReadFile(const char *filename);

	//Append current time and the given string into the log file 
	//if unit==0, write into file only, otherwise write to both the file and the screen
	static void WriteLog(const char *str,const char *filename, int unit=1);
	static void WriteLog(const char *str);

	//Create a unique file name using the given key stored in infilename
	//if overwritten == true, will not care if the file exist or not, just generate key_#start#.root
	//otherwise, it will check if  key_#start#.root exist or not, if exist, increase
	//start by 1 and keep searching till find a name which is not exist.
	//then return this name. for example: 
	//CreateFileName("key.root",outfilename, 0, true) will return key_00.root
	//In short, this routine will find the last dot '.' and insert an index '_##' in fornt of the dot.
	static void CreateFileName(const char *infilename,char outfilename[], int start, bool overwritten);

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
	static void CreateFileName(char *instr,char *outstr,int outtype);

	//check the path for the given string to see if it exist or not
	static bool CheckPath(const char *filename);
	static bool CheckPath(std::string &filename);

	//remove space or tab \t at both ends of the string
	static void Trim(std::string &astr);
	static void Trim(char *astr);

	//replace the oldsubstr with newsubstr in everywhere  
	static void ReplaceAll(char *targetstr,const char *oldsubstr, const char *newsubstr);
	static void ReplaceAll(string& targetstr,const char *oldsubstr, const char *newsubstr);

private:
	//print the help menu for the given brach, exit from program if needed
	void PrintBranch(branch_t &br,bool quit=false);

	//Build maps using the branch vector 
	void BuildMap();


	//split a sentence (string) into words by the following characters " \t,;"
	//allow no more than 10 words, return the number of words
	size_t SplitList(const string str, string list[],size_t allownum=10);

	//remove space or tab \t at both ends of the string
	std::string TrimString(char *astr);
	std::string TrimString(std::string &astr);

	//Keep reading lines from stream ifs till find one valid content, treat ererything after the first
	//'#' as comment, ignoring space ' ' and tab key '\t' at both end of the content
	std::string ReadOneValidLine(ifstream &ifs);

	//Read the Usage configuration file, it must have the format described in the head of this file
	void ReadIni(const char *inifile="Usage.ini");

	//process all reserved options, such as switch or "inputfile"
	//this routine will be invoked for each '-' found in the command line  
	void ProcessReservedOption(branch *br);

private:
	vector <branch_t> VBranch;
	map<string,int> MapBr_i;
	map<string,float> MapBr_f;
	map<string,string> MapBr_s;
	map<string,string> MapParam_s;
	std::string Executable;
};

typedef UsageManager UsageT;
typedef UsageManager BonusUsage;

#endif //_USAGET_H_

