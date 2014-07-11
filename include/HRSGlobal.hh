#ifndef _HRSGLOBAL_H_
#define _HRSGLOBAL_H_

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <time.h>
#include <string>
#include<math.h>
using namespace std;

//return true if file  exist
bool CheckFile(const char *filename);
bool CheckFile(std::string &val);

//use this one to open log file and put the time stamp
void InitLogStream(ofstream &logstream,const char *filename="_G4Sim_Log.txt");

//if unit==0, write into file only, otherwise write to both the file and the screen
void WriteLog(const char *str,const char *filename="_G4Sim_Log.txt", int unit=1);
bool CreateFileName(char *instr,char *outstr,int type=-1);
void ReplaceAll(char *targetstr,const char *oldsubstr, const char *newsubstr);
void ReplaceAll(string& targetstr,const char *oldsubstr, const char *newsubstr);


#endif //_HRSGLOBAL_H_

