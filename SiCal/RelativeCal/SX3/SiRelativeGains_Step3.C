//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Relative calibration of Si gains for SX3 Step 3
// See readme.md for general instructions.
// Usage: root -l SiRelativeGains_Step3.C+
//
// Edited by : John Parker , 2016Jan22
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
void SiRelativeGains_Step3(void)
{
  using namespace std;

  TFile *f1 = new TFile("/data0/manasta/OrganizeRaw_files/run930_931_nospacer_X3step2_divideback.root");
  if ( !f1->IsOpen() ){
    cout << "Error: Root file does not exist\n";
    exit(EXIT_FAILURE);
  }
  
  //Input the .dat file used by Main.cpp to generate the .root file given above
  Gains gains;
  gains.Load("saves/X3RelativeGains_Step2.dat");
  
  ofstream outfile;
  outfile.open("saves/X3RelativeGains_Step3.dat");
  
  TCanvas *can = new TCanvas("can","can",800,600);

  BadDetectors bad;
  bad.count=0;
  GainMatch gainmatch;

  for (Int_t DetNum=4; DetNum<28; DetNum++) {
    for (Int_t BackChNum=0; BackChNum<4; BackChNum++){
      Int_t FrontChNum = 1;
      if ( DetNum==8 ){
      	BackChNum = 3;
      }
      TH2F *hist = NULL;
      TString hname=Form("back_vs_front%i_back%i",DetNum,BackChNum);
      hist = (TH2F*)f1->Get(hname.Data());
      if (hist==NULL) {
	cout << hname << " histogram does not exist\n";
	bad.Add(DetNum,FrontChNum,BackChNum);
	continue;
      }

      Double_t gain = gainmatch.MyFit1(hist,can,kFALSE);
      gains.old[DetNum-4][BackChNum] = gains.old[DetNum-4][BackChNum]/gain;
    }
    for (Int_t i=0; i<12; i++){
      outfile << DetNum << "\t" << i << "\t" << gains.old[DetNum-4][i] << endl;
    }
  }
  outfile.close();
  bad.Print();
}
