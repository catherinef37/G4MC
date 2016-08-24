#include <fstream>
#include <iostream>
#include <sstream>

#include "string.h"
#include "math.h"

#include "TChain.h"
#include "TH1F.h"
#include "TMinuit.h"
#include "TSystem.h"
#include "TPad.h"

using std::cout;
using std::ifstream;
using std::ofstream;
using std::ostringstream;

//Headers
void fcn(Int_t &nparam, Double_t *gin, Double_t &outValue, Double_t *inValues, Int_t flag);
Double_t minFunction(Double_t X, Double_t Y);
void printarray(string *arg, int length);
void changeFieldCmds(Double_t BL1Field, Double_t BL2Field, Double_t BL4Field);
void insertCommands(string inFile, string outFile, string *cmd, string *value, const int nElem);

//This function minimizes the minFunction value by adjusting the three quad fields.
void fit(){
	//Create minimizer and set its function.
	TMinuit *minimizer = new TMinuit(5);
	minimizer -> SetFCN(fcn);

	//Set initial parameter values, step sizes, and parameter limits.
	Double_t snakemagnumber = 0.83762;
	Double_t startValues[3] = {0.260387/snakemagnumber, 
				   0.93528 /snakemagnumber,
				   1.15762 /snakemagnumber};
	Double_t stepSizes[3] = {1E-4, 1E-4, 1E-4};
	Int_t ierflg = 0;
	
	minimizer -> mnparm(0, "BL1Field", startValues[0], stepSizes[0],
			    0,0, ierflg);
	minimizer -> mnparm(1, "BL2Field", startValues[1], stepSizes[1],
			    0,0, ierflg);
	minimizer -> mnparm(2, "BL4Field", startValues[2], stepSizes[2],
			    0,0, ierflg);

	//Minimizing process
	Double_t arglist[10];
	arglist[0] = 500;
	arglist[1] = 1.;
	minimizer -> mnexcm("MIGRAD", arglist, 2, ierflg);

	//Print results so far
	minimizer -> SetPrintLevel(1);
	Double_t amin, edm, errdef;
	Int_t nvpar, nparx, icstat;
	minimizer -> mnstat(amin, edm, errdef, nvpar, nparx, icstat);
	minimizer -> mnprin(1, amin);
}

void fcn(Int_t &nparam, Double_t *gin, Double_t &outValue, Double_t *inValues, Int_t flag){
	//Change magnetic fields
	changeFieldCmds(inValues[0], inValues[1], inValues[2]);
	gSystem -> CopyFile("test1.mac", "beam.mac", kTRUE);
	
	//Run simulation
	gSystem -> Exec("./G4MC -m 1 ./beam.mac");	

	//Read in tree
	TChain *T = new TChain("track0");
	T -> Add("g4mc_out_prex_minimize.root");
  
	//Find value of sqrt(X_RMS^2 + Y_RMS^2)
	T -> Draw("X_fp_tr >> X_hist", "X_fp_tr > -1 & X_fp_tr < 1 & abs(P_sen_tr - 1.063) < 0.005");
	TH1F *X_hist = (TH1F*)gPad -> GetPrimitive("X_hist");
	T -> Draw("Y_fp_tr >> Y_hist", "Y_fp_tr > -1 & Y_fp_tr < 1 & abs(P_sen_tr - 1.063) < 0.005");
	TH1F *Y_hist = (TH1F*)gPad -> GetPrimitive("Y_hist");

	outValue = minFunction(X_hist->GetRMS(), Y_hist->GetRMS()); 
	cout << "X_hist = " << X_hist->GetRMS() << "\n";
	cout << "Y_hist = " << Y_hist->GetRMS() << "\n";
	cout << "outValue = " << outValue << "\n";

	//Delete histograms and TChain
	delete Y_hist;
	delete X_hist;
	delete T;
}

//Function to minimize
Double_t minFunction(Double_t X, Double_t Y){
	return sqrt(X*X + Y*Y);
}

//Prints out an array of strings with 'length' elements.
void printarray(string *arg, int length){
	for(int n = 0; n<length; ++n)
		cout << *(arg+n) << ' ';
}

void changeFieldCmds(Double_t BL1Field, Double_t BL2Field, Double_t BL4Field){
  
  	const int SIZE = 3;
  	string commands[SIZE] = {"/field/setFZBL1Field", 
				 "/field/setFZBL2Field",
				 "/field/setFZBL4Field"};
  	string values[SIZE];

	//Store input double field values as strings in the array
  	ostringstream s;
	s << BL1Field << "";
	values[0] = s.str() + " tesla";
	s.str("");
	s.clear();
	s << BL2Field << "";
	values[1] = s.str() + " tesla";
	s.str("");
	s.clear();
	s << BL4Field << "";
	values[2] = s.str() + " tesla";
	
	//Change macro to have new commands
 	cout << "Commands are: ";
	printarray(commands, SIZE); 
	cout << "\n";
  	cout << "Values are: ";
	printarray(values, SIZE);
	cout << "\n";
	
	insertCommands("beam.mac", "test1.mac", commands, values, SIZE);
}

//Replaces the specified commands in inFile and writes them to outFile.
//cmd is an array of size nElem whose elements will be compared on a line by line basis with those of inFile.
//value is an array of size nElem whose elements will not be compared, but will be used for replacement.
//cmd and value are parallel arrays.
//Lines that contain a substring of any element in cmd will be replaced by "cmd[i] value[i]" for some index i.
//Strings in cmd that are not found in inFile will be appended to the end of outFile as described above.
void insertCommands(string inFile, string outFile, string *cmd, string *value, const int nElem)
{ 
  //The line read in
  string str; 
  //The index counter
  int i;
  //Which strings have been found and replaced -- initial value is false.
  bool found[nElem];
	for(int j = 0; j < nElem; j++){
		found[j] = false;
	}

  //Create streams for input and output files -- names specified by user.
  ifstream inStream(inFile.c_str());
  ofstream outStream(outFile.c_str());
  
  if(inStream.is_open()){
    while(std::getline(inStream, str)){

	//If the line is a comment, copy it to the new file
	if(str.find("#")!=std::string::npos){
		//cout << "found comment" << str << "\n";
		outStream << str << "\n";
	}
	//Otherwise, check if the line should be replaced
	else {
	  i = 0;
	  //Find the index of the replacement line in cmd
	  while(i < nElem){
		//cout << i << " is the counter and is less than " << nElem <<  "\n";
		//cout << "Comparing " << str << " from file with " << *(cmd+i) << " at index " << i << "\n";
		if(str.find(*(cmd+i))== std::string::npos)
	  		i++;
		else 
			break;
	  }
	  //If the line is found in cmd, write the new line to outFile
	  if(i < nElem){
		//cout << "found " << str << " at array position " << i << "\n";
		outStream << *(cmd+i) << ' ' << *(value+i) << "\n";
		found[i] = true;
	  }
	  //Otherwise write the original line to outFile
	  else{
		outStream << str << "\n";
	  }
	}	  
    }
    //Append any new lines to outFile
    for(int n = 0; n < nElem; n++){
	if(found[n] == false){
	  outStream << *(cmd+n) << ' ' << *(value+n) << "\n";
	 //cout << "Appending" << *(cmd+n) << ' ' << *(value+n) << "\n";
    	}
    }
  }

  inStream.close();
  outStream.close();
}
