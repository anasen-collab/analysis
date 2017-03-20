//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//// Relative calibration of Si gains
////
//// Output file (e.g."X3RelativeGains_Slope1.dat") has the following columns:
//// Detector number, Front channel, Slope
////
////General Usage of all three steps
////
////If you loop over a subset of detectors, only the paramaters for those detectors will be written into the new .dat file
////First, run Organize.cpp with a given X3cal.dat file
////   Input organize.root and the .dat file into this code
////   This code will output another .dat file that you should define
////   root -l SiRelativeGains_Step1.C++
////Second, run Organize.cpp again with the new .dat file
////   Input the new organize.root and new .dat file into Step2.C
////   root -l SiRelativeGains_Step2.C++
////Finally, run Organize.cpp again with the new .dat file
////   Input the new organize.root and new .dat file into Step3.C
////   root -l SiRelativeGains_Step3.C++
////
//// Edited by : John Parker , 2016Jan22
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#include <TMath.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TGraph.h>
#include <TF1.h>
#include <TSpectrum.h>
#include <fstream>
#include <exception>

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Double_t MyFit(TH2F* hist, TCanvas* can){
  hist->Draw("colz");

  Int_t counter = 0;
  for (int i=1; i<hist->GetNbinsX(); i++){
    for (int j=1; j<hist->GetNbinsY(); j++){
	for (int k=0; k<hist->GetBinContent(i,j); k++){
	  counter++;
	}
      }
    }
  
  Double_t *x = new Double_t[counter];
  Double_t *y = new Double_t[counter];

  counter = 0;
  for (int i=1; i<hist->GetNbinsX(); i++){
    for (int j=1; j<hist->GetNbinsY(); j++){
      for (int k=0; k<hist->GetBinContent(i,j); k++){
	x[counter] = i*0.01;
	y[counter] = j*0.01;
	counter++;
      }
    }
  }

  TGraph *graph = new TGraph(counter,x,y);
  graph->Draw("*same");

  TF1 *fun2 = new TF1("fun2","[0]+[1]*x",0,2);
  fun2->SetParameter(0,1);
  fun2->SetParameter(1,-1);
  graph->Fit("fun2");
  can->Update();

  //cout << fun2->GetChisquare() << "  " << fun2->GetChisquare()/counter;

  //can->WaitPrimitive();

  Double_t gain = fun2->GetParameter(1);

  delete x;
  delete y;
  delete graph;
  delete fun2;

  if (gain < -3 || gain > -0.1 ){//basically if the slope is very far from one, I don't trust the fit
    gain = -1;
  }
  return gain;
}

void SiRelativeGains_Step1(void)
{
  using namespace std;

  //TFile *f1 = new TFile("run236out_nocal.root");//front
  //TFile *f1 = new TFile("/data0/nabin/ANASEN/ANASEN_NKJ/New/evt2root/run251_NSCL11_Pulser.root");//back
  TFile *f1 = new TFile("/home/lighthall/anasen/root/run1226-9m.root");

  TCanvas *can = new TCanvas("can","can",800,1100);
  
  ofstream outfile;

  Double_t average_slope = 0;
  
  outfile.open("X3RelativeGains_Step1.dat"); //output file name

  ifstream infile;
  infile.open("X3RelativeGains_Slope1.dat"); //input file name
  Int_t det=0,ch=0;
  Double_t dummy_slope = 0;
  Double_t slope[24][12];
  if (infile.is_open()){
    while (!infile.eof()){
      infile >> det >> ch >> dummy_slope;
      slope[det-4][ch] = dummy_slope;
    }
  }else{
    cout << "Infile not opened\n";
    exit(EXIT_FAILURE);
  }
  infile.close();
  
  for (Int_t DetNum=4; DetNum<28; DetNum++){
    for (Int_t FrontChNum=0; FrontChNum<4; FrontChNum++){
      TH2F *hist = NULL;
      TString hname=Form("down_vs_up%i_f%i",DetNum,FrontChNum);
      hist = (TH2F*)f1->Get(hname.Data());
      if (hist==NULL) {
	cout << hname << " histogram does not exist\n";
	continue;
      }
      
      average_slope += MyFit(hist,can);
      slope[DetNum-4][FrontChNum+4] = -slope[DetNum-4][FrontChNum+4]/average_slope;
    }
  }    
  for (Int_t i=0; i<24; i++){
    for (Int_t j=0; j<12; j++){
      outfile << i+4 << "\t" << j << "\t" << slope[i][j] << endl;
    }
  }
  outfile.close();
  delete can;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////		 

	    
	     
    
