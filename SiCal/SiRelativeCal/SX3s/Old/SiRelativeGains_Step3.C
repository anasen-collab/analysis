//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//// Relative calibration of Si gains
////
//// Output file (e.g."Sipulser_2015Dec13.dat") has the following columns:
//// MBID, CBID, ASICs_Channel, ZeroShift(offset), Voltage_per_Ch(slope)
////
//// Usage: root -l SiPulser_All.C++ (from the same directory).
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
Double_t MyFit(TH2F* hist, TCanvas *can){
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
	x[counter] = i*100;
	y[counter] = j*100;
	counter++;
      }
    }
  }

  TGraph *graph = new TGraph(counter,x,y);
  graph->Draw("*same");

  TF1 *fun2 = new TF1("fun2","[0]*x",0,10000);
  graph->Fit("fun2");
  can->WaitPrimitive();

  Double_t gain = fun2->GetParameter(0);
  delete x;
  delete y;
  delete graph;
  delete fun2;

  return gain;
}

void SiRelativeGains_Step3(void)
{
  using namespace std;

  TFile *f1 = new TFile("run236out_NoClickStep2.root");//front
  TCanvas *can = new TCanvas("can","can",800,1100);
  TH2F *hist;
  
  ofstream outfile;
  ofstream outfile2;

  outfile.open("X3RelativeGains012216_NoClickStep3.dat");

  ifstream infile;
  infile.open("X3RelativeGains012216_NoClickStep2.dat");
  Int_t det=0,ch=0;
  Double_t slope[24][12];
  if (infile.is_open()){
    while (!infile.eof()){
      infile >> det >> ch >> slope[det-4][ch];
    }
  }else{
    cout << "Infile not opened\n";
    exit(EXIT_FAILURE);
  }
  infile.close();

  Double_t gain = 0;

  for (Int_t DetNum=16; DetNum<17; DetNum++){
    for (Int_t FrontChNum=1; FrontChNum<4; FrontChNum++){
      Int_t BackChannelNumber = 0;
      if ( DetNum==11 ){
	BackChannelNumber = 1;
      }
      if (DetNum==22 && (FrontChNum==2 || FrontChNum==3) ){
	continue;
      }
      if (DetNum==25 && FrontChNum==2){
	continue;
      }
	
      hist = (TH2F*)f1->Get(Form("back_vs_front%i_%i_%i",DetNum,FrontChNum,BackChannelNumber));

      gain = MyFit(hist,can);
      slope[DetNum-4][FrontChNum+4] = slope[DetNum-4][FrontChNum+4]*gain;
      slope[DetNum-4][FrontChNum+8] = slope[DetNum-4][FrontChNum+8]*gain;

    }
  }
  for (Int_t i=0; i<24; i++){
    for (Int_t j=0; j<12; j++){
      outfile << i+4 << "\t" << j << "\t" << slope[i][j] << endl;
    }
  }
  delete can;
  outfile.close();
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////		 

	    
	     
    
