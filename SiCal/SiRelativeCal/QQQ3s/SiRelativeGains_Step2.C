//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Relative calibration of Si gains for QQQ Step 2
// See readme.md for general instructions.
// Usage: root -l SiRelativeGains_Step2.C+
//
// Edited by : John Parker , 2016Jan22
// Developed by : Jon Lighthall, 2017.01
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
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
//Methods
#include "SiRelativeGains.h"
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void SiRelativeGains_Step2(void)
{
  using namespace std;
  TFile *f1 = new TFile("/home/lighthall/anasen/root/run1226-9mQ2z.root");
  if ( !f1->IsOpen() ){
    cout << "Error: Root file does not exist\n";
    exit(EXIT_FAILURE);
  }

  //Input the .dat file used by Main.cpp to generate the .root file given above
  Gains gains;
  gains.Load("saves/QQQRelativeGains_Step1.dat");
    
  gains.Open("saves/QQQRelativeGains_Step2");  
  
  TCanvas *can = new TCanvas("can","can",1362,656);
  can->SetWindowPosition(0,63);
  
  BadDetectors bad;
  bad.count=0;
  GainMatch gainmatch;
  
  for (Int_t DetNum=0; DetNum<ndets; DetNum++) {
    for (Int_t BackChNum=1; BackChNum<16; BackChNum++) {
      Int_t FrontChNum = 0;

      TH2F *hist = NULL;
      TString hname=Form("Q3_back_vs_front%i_%i_%i",DetNum,FrontChNum,BackChNum);
      hist = (TH2F*)f1->Get(hname.Data());
      if (hist==NULL) {
	cout << hname << " histogram does not exist\n";
	bad.Add(DetNum,FrontChNum,BackChNum);
	outfile2 << DetNum << "\t" << FrontChNum << "\t" << BackChNum << "\t"
		 << left << fixed << setw(8) << gains.old[DetNum][BackChNum] << "\t"
		 << left << fixed << setw(8) << "N/A\t\t"
		 << left << fixed << setw(8) << 0 << endl;
	gains.old[DetNum][BackChNum] = 0;
	continue;
      }

      Double_t gain = gainmatch.Fit6(hist,can); //set fit method here
      printf(" Previous gain = %f \t Slope = %f \t New gain = %f\n",gains.old[DetNum][BackChNum],gain, gains.old[DetNum][BackChNum]/gain);
      outfile2 << DetNum << "\t" << FrontChNum << "\t" << BackChNum << "\t"
      	       << left << fixed << setw(8) << gains.old[DetNum][BackChNum] << "\t"
      	       << left << fixed << setw(8) << gain << "\t"
      	       << left << fixed << setw(8) << gains.old[DetNum][BackChNum]/gain << endl;
      gains.old[DetNum][BackChNum] = gains.old[DetNum][BackChNum]/gain;
    }
    for (Int_t i=0; i<32; i++){
      outfile << DetNum << "\t" << i << "\t" << gains.old[DetNum][i] << endl;
    }
  }
  outfile.close();
  outfile2.close();
  bad.Print();
}
