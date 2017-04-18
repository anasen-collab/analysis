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

  //TFile *f1 = new TFile("/home/lighthall/anasen/root/run1226-9m.root");//in gas
  //TFile *f1 = new TFile("/home/lighthall/anasen/root/run1255-61mQ2S1f.root");//all proton scattering
  TFile *f1 = new TFile("/home/lighthall/anasen/root/run1255-7mQ2S1.root");//10MeV only
  if ( !f1->IsOpen() ){
    cout << "Error: Root file does not exist\n";
    exit(EXIT_FAILURE);
  }
  
  //Input the .dat file used by Main.cpp to generate the .root file given above
  Gains gains;
  gains.Load("saves/X3RelativeGains_Step1_10MeV_rerun.dat");
  gains.Save("saves/X3RelativeGains_Step1");
  
  TCanvas *can = new TCanvas("can","can",800,600);

  BadDetectors bad;
  bad.count=0;
  GainMatch gainmatch;

  for (Int_t DetNum=4; DetNum<28; DetNum++) {
    //if(DetNum!=11) continue;
    for (Int_t FrontChNum=0; FrontChNum<4; FrontChNum++) {
      TH2F *hist = NULL;
      //TString hname=Form("down_vs_updivideBack%i_f%i",DetNum,FrontChNum); //normalized
      TString hname=Form("down_vs_up%i_f%i",DetNum,FrontChNum);           //unnormalized
      hist = (TH2F*)f1->Get(hname.Data());
      Int_t deti=DetNum-4;
      Int_t frti=FrontChNum+8;
      if (hist==NULL) {
	cout << hname << " histogram does not exist\n";
	bad.Add(DetNum,FrontChNum);
	outfile2 << DetNum << "\t" << frti << "\t"
		 << left << fixed << setw(8) <<gains.old[deti][frti] << "\t"
		 << left << fixed << setw(8) << "N/A\t\t"
		 << left << fixed << setw(8) << 0 << endl;
	gains.old[deti][frti]=0;
	continue;
      }
      
      Double_t gain = gainmatch.Fit4(hist,can,-1);
      
      printf("Previous gain = %f \t Slope = %f \t New gain = %f\n",gains.old[deti][frti],gain, -gains.old[deti][frti]/gain);
      outfile2 << DetNum << "\t" << frti << "\t"
       	       << left << fixed << setw(8) <<gains.old[deti][frti] << "\t"
       	       << left << fixed << setw(8) << gain << "\t"
	       << left << fixed << setw(8) << -gains.old[deti][frti]/gain << endl;
      gains.old[deti][frti] = -gains.old[deti][frti]/gain;
    }
    for (Int_t i=0; i<12; i++){
      outfile << DetNum << "\t" << i << "\t" << gains.old[DetNum-4][i] << endl;
    }
  }
  outfile.close();
  outfile2.close();
  bad.Print();
}
