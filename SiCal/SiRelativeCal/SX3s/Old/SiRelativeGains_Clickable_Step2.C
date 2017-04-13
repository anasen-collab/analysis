//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////root -l SiRelativeGains_Clickable_Step2.C++
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

void SiRelativeGains_Clickable_Step2(void)
{
  using namespace std;

  //TFile *f1 = new TFile("run236out_Step1.root");
  TFile *f1 = new TFile("/home/lighthall/anasen/root/run1255-61m.root");//all proton scattering

  ifstream infile;
  infile.open("../saves/X3RelativeGains_09182016_Slope1.dat"); 
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
  outfile.open("saves/X3RelativeGains012216_Step2_click.dat");
  
  TCanvas *can = new TCanvas("can","can",800,600);

  Int_t bad_det[288];
  Int_t bad_front[288];
  Int_t bad_back[288];
  Int_t count_bad = 0;

  for (Int_t DetNum=4; DetNum<28; DetNum++) {
    for (Int_t BackChNum=0; BackChNum<4; BackChNum++){
      if ( DetNum==11 && (BackChNum==0 || BackChNum==2 || BackChNum==3) ){
	continue;
      }
      if (DetNum==18 && BackChNum==3){
	continue;
      }
      if (DetNum==24 && BackChNum==3){
	continue;
      }
      if (DetNum==26 && (BackChNum==1 || BackChNum==2 || BackChNum==3) ){
	continue;
      }
      if (DetNum==27 && (BackChNum==2 || BackChNum==3) ){
	continue;
      }
      TH2F *hist = NULL;
      hist = (TH2F*)f1->Get(Form("back_vs_front%i_0_%i",DetNum,BackChNum));
      if (hist==NULL){
	cout << "Histo does not exist\n";
	bad_det[count_bad] = DetNum;
	bad_front[count_bad] = FrontChNum;
	bad_back[count_bad] = BackChNum;
	count_bad++;
	continue;
      }  
      
      Double_t average_slope = MyFit(hist,can);
      slope[DetNum-4][BackChNum] = 1./average_slope*slope[DetNum-4][BackChNum];
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
