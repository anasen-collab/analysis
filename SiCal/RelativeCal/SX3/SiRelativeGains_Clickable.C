//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////root -l SiRelativeGains_Clickable.C++
//// Edited by : John Parker , 2016Jan24
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

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Double_t MyFit(TH2F* hist, TCanvas *can){
  hist->Draw("colz");
  hist->GetXaxis()->SetRange(0,120);
  hist->GetYaxis()->SetRange(0,120);

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

void SiRelativeGains_Clickable(void)
{
  using namespace std;
  //TFile *f1 = new TFile("/data0/manasta/OrganizeRaw_files/run930_931_nospacer_X3slope1_divideback.root"); // root file name created with Main.cpp
  TFile *f1 = new TFile("/home/lighthall/anasen/root/run1255-61m.root");//all proton scattering
  if ( !f1->IsOpen() ){
    cout << "Error: Root file does not exist\n";
    exit(EXIT_FAILURE);
  }
  
  ifstream infile;
  infile.open("saves/X3RelativeGains_09182016_Slope1.dat"); //input file name
  Int_t det=0,ch=0;
  Double_t slope[24][12];
  Double_t dummy_slope = 0;
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

  ofstream outfile;
  outfile.open("saves/X3RelativeGains09192016_Step1_click.dat"); //output file name
  
  TCanvas *can = new TCanvas("can","can",800,600);

  Int_t bad_det[288];
  Int_t bad_front[288];
  Int_t bad_back[288];
  Int_t count_bad = 0;

  for (Int_t DetNum=4; DetNum<28; DetNum++){
    for (Int_t FrontChNum=0; FrontChNum<4; FrontChNum++){
      TH2F *hist = NULL;
      hist = (TH2F*)f1->Get(Form("down_vs_up_divideBack%i_front%i",DetNum,FrontChNum));//normalized
      //hist = (TH2F*)f1->Get(Form("down_vs_up_%i_front%i",DetNum,FrontChNum));//unnormalized
      if (hist==NULL){
	cout << "Histo does not exist\n";
	bad_det[count_bad] = DetNum;
	bad_front[count_bad] = FrontChNum;
	count_bad++;
	continue;
      }
      
      Double_t average_slope = MyFit(hist,can);
      slope[DetNum-4][FrontChNum+8] = -slope[DetNum-4][FrontChNum+8]/average_slope;
    }
    for (Int_t i=0; i<12; i++){
      outfile << DetNum << "\t" << i << "\t" << slope[DetNum-4][i] << endl;
    }
  }
  outfile.close();
  cout << "List of bad detectors:\n";
  for (int i=0; i<count_bad; i++){
    cout << bad_det[i] << "  " << bad_front[i] << endl;
  }
}
