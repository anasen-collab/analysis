//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Relative calibration of Si gains for SX3 Step 3
// See readme.md for general instructions.
// Usage: root -l SiRelativeGains_Step3.C+
//
// Edited by : John Parker , 2016Jan22
// Developed by : Jon Lighthall, 2017.04
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//C++
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

  TFile *f1 = new TFile("/home/lighthall/anasen/root/run1226-9mQ2S3.root");
  if ( !f1->IsOpen() ){
    cout << "Error: Root file does not exist\n";
    exit(EXIT_FAILURE);
  }
  
  //Input the .dat file used by Main.cpp to generate the .root file given above
  Gains gains;
  gains.Load("saves/X3RelativeGains_Step3_170418.dat");
  gains.Save("saves/X3RelativeGains_Step3");
  
  TCanvas *can = new TCanvas("can","can",800,600);

  BadDetectors bad;
  GainMatch gainmatch;

  for (Int_t DetNum=4; DetNum<28; DetNum++) {
    for (Int_t BackChNum=0; BackChNum<4; BackChNum++){
      TH2F *hist = NULL;
      TString hname=Form("back_vs_front%i_b%i",DetNum,BackChNum);
      hist = (TH2F*)f1->Get(hname.Data());
      if (hist==NULL) {
	cout << "****" << hname << " histogram does not exist\n";
	bad.Add(DetNum,-1,BackChNum);
	gains.Add(DetNum-4,BackChNum,0,0);
	continue;
      }

      Double_t gain = gainmatch.Fit4(hist,can,1);
      gains.Add(DetNum-4,BackChNum,gain,1.0/gain);
    }
    for (Int_t i=0; i<12; i++){
      outfile << DetNum << "\t" << i << "\t" << gains.old[DetNum-4][i] << endl;
    }
  }
  bad.Print();
}
