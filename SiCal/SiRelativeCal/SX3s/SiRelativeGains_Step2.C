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
void SiRelativeGains_Step2(void)
{
  using namespace std;

  TFile *f1 = new TFile("/home/lighthall/anasen/root/run1226-9mQ2S2.root");
  //TFile *f1 = new TFile("/home/lighthall/anasen/root/run1255-7mQ2.root");//10MeV only
  if ( !f1->IsOpen() ){
    cout << "Error: Root file does not exist\n";
    exit(EXIT_FAILURE);
  }
  
  //Input the .dat file used by Main.cpp to generate the .root file given above
  Gains gains;
  gains.Load("saves/X3RelativeGains_Step2_no1.dat");
  gains.Open("saves/X3RelativeGains_Step2");  
 
  TCanvas *can = new TCanvas("can","can",800,600);

  BadDetectors bad;
  bad.count=0;
  GainMatch gainmatch;

  for (Int_t DetNum=4; DetNum<28; DetNum++) {
    //if(DetNum>8) continue;
    for (Int_t FrontChNum=0; FrontChNum<4; FrontChNum++) {

      Int_t BackChNum = 0;   // some if-statements that differ between each data set
      //if(DetNum==8) {
      //BackChNum = 3;
      //}
      TH2F *hist = NULL;
      TString hname=Form("back_vs_front%i_%i_%i",DetNum,FrontChNum,BackChNum);
      hist = (TH2F*)f1->Get(hname.Data());
      if (hist==NULL) {
	cout << "****" << hname << " histogram does not exist\n";
	bad.Add(DetNum,FrontChNum,BackChNum);
	outfile2 << DetNum << "\t" << FrontChNum+4 << "\t"
		 << left << fixed << setw(8) <<gains.old[DetNum-4][FrontChNum+4] << "\t"
		 << "N/A     \t" << 0.0 << "\t" << 0.0 << endl;
	gains.old[DetNum-4][FrontChNum+4] = 0;
	gains.old[DetNum-4][FrontChNum+8] = 0;
	continue;
      }

      Double_t gain = gainmatch.Fit4(hist,can,1.0);
      printf("Previous gain = %f \t Slope = %f \t New gain = %f\n",gains.old[DetNum-4][FrontChNum+4],gain, -gains.old[DetNum-4][FrontChNum+4]/gain);
      outfile2 << DetNum << "\t" << FrontChNum+4 << "\t"
       	       << left << fixed << setw(8) <<gains.old[DetNum-4][FrontChNum+4] << "\t"
       	       << left << fixed << setw(8) << gain << "\t"
	       << left << fixed << setw(8) << -gains.old[DetNum-4][FrontChNum+4]*gain << "\t"
      	       << left << fixed << setw(8) << -gains.old[DetNum-4][FrontChNum+8]*gain << endl;
      gains.old[DetNum-4][FrontChNum+4] = gains.old[DetNum-4][FrontChNum+4]*gain;
      gains.old[DetNum-4][FrontChNum+8] = gains.old[DetNum-4][FrontChNum+8]*gain;
    }
    for (Int_t i=0; i<12; i++){
      outfile << DetNum << "\t" << i << "\t" << gains.old[DetNum-4][i] << endl;
    }
  }
  outfile.close();
  bad.Print();
}
