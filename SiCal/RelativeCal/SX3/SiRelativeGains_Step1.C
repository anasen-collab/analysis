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
//ROOT
#include <TMath.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
//Methods
#include "SiRelativeGains.h"
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void SiRelativeGains_Step1(void) {
  using namespace std;

  f1 = new TFile("/home/lighthall/anasen/root/main/10MeVmQ2S3.root");//10MeV only
  if ( !f1->IsOpen() ){
    cout << "Error: Root file does not exist\n";
    exit(EXIT_FAILURE);
  }
  
  //Input the .dat file used by Main.cpp to generate the .root file given above
  Gains gains;
  gains.Load("saves/X3RelativeGains_Step3_170525.dat");
  gains.Save("saves/X3RelativeGains_Step1");
  
  BadDetectors bad;
  GainMatch gainmatch;

  for (Int_t DetNum=4; DetNum<ndets+4; DetNum++) {
    //if(DetNum!=12) continue;
    for (Int_t FrontChNum=0; FrontChNum<4; FrontChNum++) {
      TH2F *hist = NULL;
      //TString hname=Form("down_vs_up_divideBack%i_f%i",DetNum,FrontChNum); //normalized
      TString hname=Form("down_vs_up%i_f%i",DetNum,FrontChNum);           //unnormalized
      hist = (TH2F*)f1->Get(hname.Data());
      if (hist==NULL) {
	cout << hname << " histogram does not exist\n";
	bad.Add(DetNum,FrontChNum);
	gains.Add(DetNum-4,FrontChNum+8,0,0);
	continue;
      }
      Double_t gain = gainmatch.Fit4(hist,-gains.old[DetNum-4][FrontChNum+8]);
      //gain = gainmatch.Fit9(DetNum,FrontChNum);//may only work if using cut from Fit4
      gains.Add(DetNum-4,FrontChNum+8,gain,-1.0/gain);
    }
    for (Int_t i=0; i<12; i++){
      outfile << DetNum << "\t" << i << "\t" << gains.old[DetNum-4][i] << endl;
    }
  }
  bad.Print();
}
