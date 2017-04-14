//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Relative calibration of Si gains for SX3 Step 2
// See readme.md for general instructions.
// Usage: root -l SiRelativeGains_Step2.C+
//
// Edited by : John Parker , 2016Jan22
// Edited by : Maria Anastasiou, 2016Sept20
// Developed by : Jon Lighthall, 2017.04
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//C++
#include <fstream>
#include <exception>
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
void SiRelativeGains_Step2(void)
{
  using namespace std;

  //TFile *f1 = new TFile("/data0/manasta/OrganizeRaw_files/run930_931_nospacer_X3step1_divideback.root");
  TFile *f1 = new TFile("/home/lighthall/anasen/root/run1226-9mQ2S1.root");
  if ( !f1->IsOpen() ){
    cout << "Error: Root file does not exist\n";
    exit(EXIT_FAILURE);
  }
  
  //Input the .dat file used by Main.cpp to generate the .root file given above
  ifstream infile;
  infile.open("saves/X3RelativeGains_Step1.dat");
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

  ofstream outfile;
  outfile.open("saves/X3RelativeGains_Step2.dat"); 
  
  TCanvas *can = new TCanvas("can","can",800,600);

  Int_t bad_det[288];
  Int_t bad_front[288];
  Int_t bad_back[288];
  Int_t count_bad = 0;

  GainMatch gainmatch;

  for (Int_t DetNum=4; DetNum<28; DetNum++) {
    for (Int_t FrontChNum=0; FrontChNum<4; FrontChNum++){
      Int_t BackChNum = 0;   // some if-statements that differ between each data set
      if ( DetNum==8 ){
      	BackChNum = 3;
      }
      TH2F *hist = NULL;
      TString hname=Form("back_vs_front%i_%i_%i",DetNum,FrontChNum,BackChNum);
      hist = (TH2F*)f1->Get(hname.Data());
      if (hist==NULL) {
	cout << hname << " histogram does not exist\n";
	bad_det[count_bad] = DetNum;
	bad_front[count_bad] = FrontChNum;
	bad_back[count_bad] = BackChNum;
	count_bad++;
	continue;
      }

      Double_t gain = gainmatch.Fit1(hist,can,kFALSE);
      slope[DetNum-4][FrontChNum+4] = slope[DetNum-4][FrontChNum+4]*gain;
      slope[DetNum-4][FrontChNum+8] = slope[DetNum-4][FrontChNum+8]*gain;
    }
    for (Int_t i=0; i<12; i++){
      outfile << DetNum << "\t" << i << "\t" << slope[DetNum-4][i] << endl;
    }
  }
  outfile.close();
  cout << "List of bad detectors:\n";
  for (Int_t i=0; i<count_bad; i++){
    cout << bad_det[i] << "  " << bad_front[i] << "  " << bad_back[i] << endl;
  }
}
