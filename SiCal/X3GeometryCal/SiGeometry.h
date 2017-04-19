//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// 
// 
//
// Developed by : Jon Lighthall, 2017.04
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//C++
#include <fstream>
#include <exception>
#include <time.h>
//ROOT
#include <TMath.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TGraph.h>
#include <TF1.h>
#include <TSpectrum.h>
#include <TCutG.h>
#include <TVector.h>
#include <TLegend.h>
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
using namespace std;
const Int_t ndets=24;
const Int_t nchan=12;
const Int_t range=4096*4;
ofstream outfile;

class Gains {
 public:
  void Save(TString);
};

class Time {
 public:
  char stamp[80];
  void Get();
};

class BadDetectors {
 public:
  Int_t det[ndets*nchan];
  Int_t front[ndets*nchan];
  Int_t back[ndets*nchan];
  Int_t count;
  void Add(Int_t,Int_t,Int_t BackChNum=-1);
  void Print();
};

void Gains::Save(TString fname) {
  Time time;
  time.Get();
  outfile.open(Form("%s_%s.dat",fname.Data(),time.stamp));
  outfile << "DetNum\tFrontCh\tBackCh\t\n";
}

void Time::Get() {
  time_t rawtime;
  struct tm * timeinfo;
  time (&rawtime);
  timeinfo = localtime (&rawtime);
  strftime (stamp,80,"%y%m%d.%H%M%S",timeinfo);
  char date[80];
  strftime (date,80,"%a %b %d %Y at %H:%M:%S. ",timeinfo);
  printf("The date is %s",date);
  printf("Time stamp is %s\n",stamp);
}

void BadDetectors::Add(Int_t DetNum, Int_t FrontChNum, Int_t BackChNum) {
  det[count] = DetNum;
  front[count] = FrontChNum;
  back[count] = BackChNum;
  count++;
}

void BadDetectors::Print() {
  printf("List of bad detectors:\n");
  printf(" DetNum\tFrontCh\tBackCh\n");
  for (Int_t i=0; i<count; i++){
    cout << " " << det[i] << "\t" << front[i] << "\t" << back[i] << endl;
  }
  for(int i=0;i<3;i++) {//print beeps at end of program
    printf(" beep!\a\n");
    sleep(1);
  }
}
