//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Fit methods for relative calibration of Si gains for SX3
// See readme.md for general instructions.
//
// Developed by : Jon Lighthall, 2017.04
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//C++
#include <fstream>
#include <exception>
#include <time.h>
#include <iomanip>
//ROOT
#include <TMath.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH2F.h>
#include <TGraph.h>
#include <TF1.h>
#include <TCutG.h>
#include <TVector.h>
#include <TLegend.h>
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
using namespace std;
const Int_t ndets=28;
Int_t range=4096*4;
ofstream outfile;
Int_t counter;

class Gains {
 public:
  Double_t old[ndets];
  void Load(TString);
  void Print();
  void Save(TString);
  void Add(Int_t,Double_t,Double_t);
  ~Gains() {
    outfile.close();
  };
};

class Time {
 public:
  char stamp[80];
  void Get();
};

class BadDetectors {
 public:
  Int_t det[ndets];
  Int_t count;
  BadDetectors() {
    count=0;
  };
  void Add(Int_t);
  void Print();
};

void Gains::Load(TString fname) {
  ifstream infile;
  printf("Loading file %s...",fname.Data());
  infile.open(fname.Data());
  Int_t det=0,ch=0;
  Double_t dummy = 0;
  if (infile.is_open()) {
    cout << "Read OK"<<endl;
    infile.ignore(100,'\n');//read in dummy line
    while (!infile.eof()){
      infile >> det >> ch >> dummy;
      old[det] = dummy;
    }
  }else{
    cout << "Error: Dat file " << fname.Data() << " does not exist\n";
    exit(EXIT_FAILURE);
  }
  infile.close();
}

void Gains::Print() {
  printf("DetNum\tGain\n");
  for (Int_t i=0; i<ndets; i++){
      printf("%d\t%f\n",i,old[i]);
  }
}

void Gains::Save(TString fname) {
  Time time;
  time.Get();
  outfile.open(Form("%s_%s.dat",fname.Data(),time.stamp));
  outfile << "DetNum\tOffset\tSlope\n";
}

void Gains::Add(Int_t DetNum,Double_t slope,Double_t new_gain) {
  if(new_gain)
    printf(" Previous gain = %f \t Slope = %f \t New gain = %f\n",old[DetNum],slope,old[DetNum]*new_gain);
  old[DetNum]*=new_gain;
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

void BadDetectors::Add(Int_t DetNum) {
  det[count] = DetNum;
  count++;
}

void BadDetectors::Print() {
  printf("List of bad detectors:\n");
  printf(" DetNum\n");
  for (Int_t i=0; i<count; i++){
    cout << " " << det[i] << endl;
  }
  for(int i=0;i<3;i++) {//print beeps at end of program
    printf(" beep!\a\n");
    sleep(1);
  }
}
