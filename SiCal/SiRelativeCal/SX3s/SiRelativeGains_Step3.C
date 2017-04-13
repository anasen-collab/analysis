//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Relative calibration of Si gains for SX3 Step 3
// See readme.md for general instructions.
// Usage: root -l SiRelativeGains_Step3.C+
//
// Edited by : John Parker , 2016Jan22
// Developed by : Jon Lighthall, 2017.04
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
Double_t MyFit1(TH2F* hist, TCanvas *can){
  hist->Draw("colz");

  Double_t x1[5] = { 2, 12, 9, 0.4, 2 };
  Double_t y1[5] = { 0.8, 9.8, 11, 1.1, 0.8 };
  TCutG *cut = new TCutG("cut",5,x1,y1);
  //cut->Draw("same");
  
  Double_t maxbinNumberX = hist->GetXaxis()->GetXmax();
  Double_t maxbinNumberY = hist->GetYaxis()->GetXmax();
  Double_t maxbinX = (maxbinNumberX/hist->GetNbinsX());
  Double_t maxbinY = (maxbinNumberY/hist->GetNbinsY());

  Int_t counter = 0;
  for (int i=1; i<hist->GetNbinsX(); i++){
    for (int j=1; j<hist->GetNbinsY(); j++){
      //if ( !cut->IsInside((Double_t)i*maxbinX,j*maxbinY) ){
      //	continue;
      //}
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
      //if ( !cut->IsInside((Double_t)i*maxbinX,j*maxbinY) ){
      //	continue;
      //}
      for (int k=0; k<hist->GetBinContent(i,j); k++){
	x[counter] = i*maxbinX;
	y[counter] = j*maxbinY;
	counter++;
      }
    }
  }

  TGraph *graph = new TGraph(counter,x,y);
  graph->Draw("*same");

  TF1 *fun2 = new TF1("fun2","[0]*x +[1]",0,10000);
  graph->Fit("fun2");
  can->Update();
  can->WaitPrimitive();

  Double_t gain = fun2->GetParameter(0);
  delete x;
  delete y;
  delete graph;
  delete fun2;

  return gain;
}

Double_t MyFit3(TH2F* hist, TCanvas *can){//cut-as-line fit; copied from Old/Clickable_Step3
  hist->Draw("colz");

  Double_t x[10];
  Double_t y[10];

  TCutG *cut;
  cut = (TCutG*)can->WaitPrimitive("CUTG");
      
  for(int n=0;n<cut->GetN()-2;n++){
    cut->GetPoint(n,x[n],y[n]);
    cout << x[n] << "\t" << y[n] << endl;
  }
	
  TGraph *graph = new TGraph(cut->GetN()-2,x,y);
  hist->Draw("colz");
  graph->Draw("*same");
	
  TF1 *fun = new TF1("fun","[0] + [1]*x",0,1);
  graph->Fit("fun");
	
  can->Update();
  can->WaitPrimitive();

  Double_t gain = fun->GetParameter(1);
      
  delete graph;
  delete fun;
      
  return gain;
}

void SiRelativeGains_Step3(void)
{
  using namespace std;

  TFile *f1 = new TFile("/data0/manasta/OrganizeRaw_files/run930_931_nospacer_X3step2_divideback.root");
  if ( !f1->IsOpen() ){
    cout << "Error: Root file does not exist\n";
    exit(EXIT_FAILURE);
  }

  //Input the .dat file used by Main.cpp to generate the .root file given above
  ifstream infile;
  infile.open("X3RelativeGains_Step2.dat");
  Int_t det=0,ch=0;
  Double_t slope[24][12];
  Double_t dummy;
  if (infile.is_open()){
    while (!infile.eof()){
      infile >> det >> ch >> dummy;
      slope[det-4][ch] = dummy;
    }
  }else{
    cout << "Infile not opened\n";
    exit(EXIT_FAILURE);
  }
  infile.close();

  ofstream outfile;
  outfile.open("X3RelativeGains_Step3.dat");
    
  TCanvas *can = new TCanvas("can","can",800,600);

  Int_t bad_det[288];
  Int_t bad_front[288];
  Int_t bad_back[288];
  Int_t count_bad = 0;

  for (Int_t DetNum=4; DetNum<28; DetNum++) {
    for (Int_t BackChNum=0; BackChNum<4; BackChNum++){
      Int_t FrontChNum = 1;
      if ( DetNum==8 ){
      	BackChNum = 3;
      }
      TH2F *hist = NULL;
      hist = (TH2F*)f1->Get(Form("back_vs_front%i_back%i",DetNum,BackChNum));
      if (hist==NULL){
	cout << "Histo does not exist\n";
	bad_det[count_bad] = DetNum;
	bad_front[count_bad] = FrontChNum;
	bad_back[count_bad] = BackChNum;
	count_bad++;
	continue;
      }

      Double_t gain = MyFit1(hist,can);
      slope[DetNum-4][BackChNum] = slope[DetNum-4][BackChNum]/gain;
    }
    for (Int_t i=0; i<12; i++){
      outfile << DetNum << "\t" << i << "\t" << slope[DetNum-4][i] << endl;
    }
  }
  outfile.close();
  cout << "List of bad detectors:\n";
  for (int i=0; i<count_bad; i++){
    cout << bad_det[i] << "  " << bad_front[i] << "  " << bad_back[i] << endl;
  }
}
