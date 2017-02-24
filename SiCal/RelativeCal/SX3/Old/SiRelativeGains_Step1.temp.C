qy//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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
#include <TCutG.h>
#include <TVector.h>

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Double_t MyFit(TH2F* hist, TCanvas* can){
  hist->Draw("colz");
  hist->GetXaxis()->SetRange(0,180);
  hist->GetYaxis()->SetRange(0,180);
  
  Double_t x1[5] = { 450, 5000, 6200, 530, 450 };
  Double_t y1[5] = { 4800, 10, 1100, 6150, 4800 };
  TCutG *cut = new TCutG("cut",5,x1,y1);
  cut->Draw("same");
  
  
  /* vector<double> x1;
  vector<double> y1;

  TCutG *cut;
  cut = (TCutG*)can->WaitPrimitive("CUTG");
  x1.resize(cut->GetN());
  y1.resize(cut->GetN());
  cut->Draw("same");
  */
  Int_t counter = 0;
  for (int i=1; i<hist->GetNbinsX(); i++){
    for (int j=1; j<hist->GetNbinsY(); j++){
      if ( !cut->IsInside((Double_t)i*0.1,(Double_t)j*0.1) ){
	continue;
      }
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
      if ( !cut->IsInside((Double_t)i*0.1,(Double_t)j*0.1) ){
	continue;
      }
      for (int k=0; k<hist->GetBinContent(i,j); k++){
	x[counter] = (Double_t)i*0.1;
	y[counter] = (Double_t)j*0.1;
	counter++;
      }
    }
  }

  TGraph *graph = new TGraph(counter,x,y);
  graph->Draw("*same");

  TF1 *fun2 = new TF1("fun2","[0]+[1]*x",0,6000);
  fun2->SetParameter(0,10);
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
  //TFile *f1 = new TFile("../../../OrganizeRaw_root/run567_051116.root");//front
  TFile *f1 = new TFile("/home/lighthall/anasen/root/run1226-9m.root");

  TCanvas *can = new TCanvas("can","can",800,600);
  
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

  Int_t bad_det[288];
  Int_t bad_front[288];
  Int_t count_bad = 0;

  for (Int_t DetNum=4; DetNum<28; DetNum++){
    for (Int_t FrontChNum=0; FrontChNum<4; FrontChNum++){
      TH2F *hist = NULL;
      TString hname=Form("down_vs_up%i_f%i",DetNum,FrontChNum);
      hist = (TH2F*)f1->Get(hname.Data());
      if (hist==NULL) {
	cout << hname << " histogram does not exist\n";
	bad_det[count_bad] = DetNum;
	bad_front[count_bad] = FrontChNum;
	count_bad++;
	continue;
      }
      
      average_slope += MyFit(hist,can);
      slope[DetNum-4][FrontChNum+8] = -slope[DetNum-4][FrontChNum+8]/average_slope;
    }
    for (Int_t i=0; i<24; i++){
      outfile << DetNum << "\t" << i << "\t" << slope[DetNum-4][i] << endl;
    }
  }
  
 
  outfile.close();
cout << "List of bad detectors:\n";
  for (int i=0; i<count_bad; i++){
    cout << bad_det[i] << "  " << bad_front[i] << endl;
 }
  delete can;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////		 

	    
	     
    
