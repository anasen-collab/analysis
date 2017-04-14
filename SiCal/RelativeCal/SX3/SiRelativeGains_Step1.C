//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Relative calibration of Si gains for SX3 Step 1
// See readme.md for general instructions.
// Usage: root -l SiRelativeGains_Step1.C+
//
// Edited by : John Parker , 2016Jan22
// Edited by : Maria Anastasiou, 2016Sept20
// Developed by : Jon Lighthall, 2017.02
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//C++
#include <fstream>
#include <exception>
#include <iomanip>
//ROOT
#include <TMath.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
//Methods
#include "SiRelativeGains.h"
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void SiRelativeGains_Step1(void)
{
  using namespace std;

  //TFile *f1 = new TFile("run236out_nocal.root");//front
  //TFile *f1 = new TFile("/data0/nabin/ANASEN/ANASEN_NKJ/New/evt2root/run251_NSCL11_Pulser.root");//back
  //TFile *f1 = new TFile("../../../OrganizeRaw_root/run567_051116.root");//front
  //TFile *f1 = new TFile("/data0/manasta/OrganizeRaw_files/run930_931_nospacer_X3slope1_divideback.root"); 
  //TFile *f1 = new TFile("/home/lighthall/anasen/root/run1226-9m.root");
  TFile *f1 = new TFile("/home/lighthall/anasen/root/run1255-61m.root");//all proton scattering
  //TFile *f1 = new TFile("/home/lighthall/anasen/root/run1255-7m.root");
  if ( !f1->IsOpen() ){
    cout << "Error: Root file does not exist\n";
    exit(EXIT_FAILURE);
  }
  
  //Input the .dat file used by Main.cpp to generate the .root file given above
  ifstream infile;
  infile.open("saves/X3RelativeGains_Slope1.dat");
  Int_t det=0,ch=0;
  Double_t slope[24][12];
  Double_t dummy = 0;
  if (infile.is_open()) {
    infile.ignore(100,'\n');//read in dummy line
    while (!infile.eof()){
      infile >> det >> ch >> dummy;
      slope[det-4][ch] = dummy;
    }
  }else{
    cout << "Error: Dat file does not exist\n";
    exit(EXIT_FAILURE);
  }
  infile.close();

  time_t rawtime;
  struct tm * timeinfo;
  char filename [80];
  time (&rawtime);
  timeinfo = localtime (&rawtime);

  ofstream outfile;
  strftime (filename,80,"saves/X3RelativeGains_Step1_%y%m%d.%H%M%S.dat",timeinfo);
  outfile.open(filename);
  outfile << "DetNum\tFrontCh\tGain\n";

   ofstream outfile2;
   strftime (filename,80,"saves/X3RelativeGains_Step1_%y%m%d.%H%M%S_back.dat",timeinfo);  // file2 may be used for diagnostics
   outfile2.open(filename);
   outfile2 << "DetNum\tFrontCh\tOld\t\tSlope\t\tNew\n";

  TCanvas *can = new TCanvas("can","can",800,600);

  Int_t bad_det[288];
  Int_t bad_front[288];
  Int_t bad_back[288];
  Int_t count_bad = 0;

  GainMatch gainmatch;

  for (Int_t DetNum=4; DetNum<28; DetNum++) {
    //if(DetNum!=21) continue;
    for (Int_t FrontChNum=0; FrontChNum<4; FrontChNum++) {
      TH2F *hist = NULL;
      //TString hname=Form("down_vs_updivideBack%i_f%i",DetNum,FrontChNum); //normalized
      TString hname=Form("down_vs_up%i_f%i",DetNum,FrontChNum);           //unnormalized
      hist = (TH2F*)f1->Get(hname.Data());
      if (hist==NULL) {
	cout << hname << " histogram does not exist\n";
	bad_det[count_bad] = DetNum;
	bad_front[count_bad] = FrontChNum;
	//bad_back[count_bad] = BackChNum;
	count_bad++;
	continue;
      }
      
      Double_t gain = gainmatch.Fit4(hist,can);
      Int_t deti=DetNum-4;
      Int_t frti=FrontChNum+8;
      printf("Previous gain = %f \t Slope = %f \t New gain = %f\n",slope[deti][frti],gain, -slope[deti][frti]/gain);
      outfile2 << DetNum << "\t" << frti << "\t"
       	       << left << fixed << setw(8) <<slope[deti][frti] << "\t"
       	       << left << fixed << setw(8) << gain << "\t"
	       << left << fixed << setw(8) << -slope[deti][frti]/gain << endl;
      slope[deti][frti] = -slope[deti][frti]/gain;
    }
    for (Int_t i=0; i<12; i++){
      outfile << DetNum << "\t" << i << "\t" << slope[DetNum-4][i] << endl;
    }
  }
  outfile.close();
  outfile2.close();
  cout << "List of bad detectors:\n";
  for (Int_t i=0; i<count_bad; i++){
    cout << bad_det[i] << "  " << bad_front[i] << "  " << bad_back[i] << endl;
  }
}
