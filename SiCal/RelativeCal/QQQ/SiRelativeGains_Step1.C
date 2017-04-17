//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Relative calibration of Si gains for QQQ Step 1
// See readme.md for general instructions.
// Usage: root -l SiRelativeGains_Step1.C+
//
// Edited by : John Parker , 2016Jan22
// Developed by : Jon Lighthall, 2016.12
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
void SiRelativeGains_Step1(void)
{
  using namespace std;

  TFile *f1 = new TFile("/home/lighthall/anasen/root/run1226-9m.root");
  if ( !f1->IsOpen() ){
    cout << "Error: Root file does not exist\n";
    exit(EXIT_FAILURE);
  }

  //Input the .dat file used by Main.cpp to generate the .root file given above
  Gains gains;
  //gains.Load("saves/QQQRelativeGains_Step1.dat");
  gains.Load("saves/QQQRelativeGains_Slope1.dat");

  gains.Print();
    
  time_t rawtime;
  struct tm * timeinfo;
  char filename [80];
  time (&rawtime);
  timeinfo = localtime (&rawtime);

  ofstream outfile;
  strftime (filename,80,"saves/QQQRelativeGains_Step1_%y%m%d.%H%M%S.dat",timeinfo);
  outfile.open(filename);
  outfile << "DetNum\tFrontCh\tGain\n";
  
  ofstream outfile2;
  strftime (filename,80,"saves/QQQRelativeGains_Step1_%y%m%d.%H%M%S_back.dat",timeinfo);  // file2 may be used for diagnostics
  outfile2.open(filename);
  outfile2 << "DetNum\tFrontCh\tBackCh\tOld\t\tSlope\t\tNew\n";
  
  TCanvas *can = new TCanvas("can","can",1362,656);
  can->SetWindowPosition(0,63);
  
  Int_t bad_det[128];
  Int_t bad_front[128];
  Int_t bad_back[128];
  Int_t count_bad = 0;

  GainMatch gainmatch;
  
  for (Int_t DetNum=0; DetNum<ndets; DetNum++) {
    for (Int_t FrontChNum=0; FrontChNum<16; FrontChNum++) {
      Int_t BackChNum = 0;
      //for (Int_t BackChNum=0; BackChNum<16; BackChNum++) {//loop over back (diagnostic)
      if(DetNum==2)
	BackChNum = 4;
      // if(DetNum==1 && (FrontChNum==4 || FrontChNum==14)){continue;}
      //	 if(DetNum==2 && FrontChNum==11){continue;}
      //if(!((FrontChNum==0)||(FrontChNum==13)))
      //if(!((FrontChNum==13)))
      //continue;

      TH2F *hist = NULL;
      TString hname=Form("Q3_back_vs_front%i_%i_%i",DetNum,FrontChNum,BackChNum);
      hist = (TH2F*)f1->Get(hname.Data());
      if (hist==NULL) {
	cout << hname << " histogram does not exist\n";
	bad_det[count_bad] = DetNum;
	bad_front[count_bad] = FrontChNum;
	bad_back[count_bad] = BackChNum;
	count_bad++;
	outfile2 << DetNum << "\t" << FrontChNum+16 << "\t" <<BackChNum << "\t"
		 << left << fixed << setw(8) <<gains.old[DetNum][FrontChNum+16] << "\t"
		 << left << fixed << setw(8) << "N/A\t\t"
		 << left << fixed << setw(8) << 0 << endl;
	gains.old[DetNum][FrontChNum+16] = 0;
	continue;
      }

      Double_t gain = gainmatch.Fit6(hist,can); //set fit method here
      printf("Previous gain = %f \t Slope = %f \t New gain = %f\n",gains.old[DetNum][FrontChNum+16],gain, gains.old[DetNum][FrontChNum+16]*gain);
      outfile2 << DetNum << "\t" << FrontChNum+16 << "\t" <<BackChNum << "\t"
	       << left << fixed << setw(8) <<gains.old[DetNum][FrontChNum+16] << "\t"
	       << left << fixed << setw(8) << gain << "\t"
	       << left << fixed << setw(8) << gains.old[DetNum][FrontChNum+16]*gain << endl;
      gains.old[DetNum][FrontChNum+16] *= gain;
      //}//back loop
    }
    for (Int_t i=0; i<32; i++){
      outfile << DetNum << "\t" << i << "\t" << gains.old[DetNum][i] << endl;
    }
  }
  outfile.close();
  outfile2.close();
  gains.Print();
  cout << "List of bad detectors:\n";
  for (Int_t i=0; i<count_bad; i++){
    cout << bad_det[i] << "  " << bad_front[i] << "  " << bad_back[i] << endl;
  }
}
