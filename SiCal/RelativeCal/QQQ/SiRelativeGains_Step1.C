//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////Relative calibration of Si gains for QQQ
////Essentially the same progam as that for the SX3
////root -l SiRelativeGains_Step1.C+


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
Double_t MyFit(TH2F* hist, TCanvas *can) {
  hist->Draw("colz");
   
  Double_t x1[12] = {1450, 630, 3150, 6340, 9200, 10540, 13200, 13600, 11670, 7550, 2400, 1450};
  Double_t y1[12] = {800, 2250, 4900, 8050, 10380, 11520, 13800, 12600, 8700, 5400, 1100, 800};
  TCutG *cut = new TCutG("cut",12,x1,y1);
  //cut = (TCutG*)can->WaitPrimitive("CUTG");
  cut->Draw("same");
  
  Double_t maxbinNumberX = hist->GetXaxis()->GetXmax();
  Double_t maxbinNumberY = hist->GetYaxis()->GetXmax();
  Double_t maxbinX = (maxbinNumberX/hist->GetNbinsX());
  Double_t maxbinY = (maxbinNumberY/hist->GetNbinsY());

  if(0) {
    printf("axis range:\n");
    cout << " maxbinNumberX " << maxbinNumberX << endl;
    cout << " maxbinNumberY " << maxbinNumberY << endl;
    printf("number of bins:\n");
    cout << " binX " << hist->GetNbinsX() << endl;
    cout << " binY " << hist->GetNbinsY() << endl;
    printf("bin width:\n");
    cout << " maxbinsX" << maxbinX << endl;
    cout << " maxbinsY" << maxbinY << endl;
  }
  
  Int_t counter = 0;
  for (int i=1; i<hist->GetNbinsX(); i++) {
    for (int j=1; j<hist->GetNbinsY(); j++) {
      if ( !cut->IsInside((Double_t)i*maxbinX,j*maxbinY)) {
      continue;
       }
      for (int k=0; k<hist->GetBinContent(i,j); k++) {
	counter++;
      }
    }
  }

  Double_t *x = new Double_t[counter];
  Double_t *y = new Double_t[counter];

  counter = 0;
  for (int i=1; i<hist->GetNbinsX(); i++){
    for (int j=1; j<hist->GetNbinsY(); j++){
      if ( !cut->IsInside((Double_t)i*maxbinX,j*maxbinY) ){
      continue;
      }
      for (int k=0; k<hist->GetBinContent(i,j); k++){
	x[counter] = i*maxbinX; 
	y[counter] = j*maxbinY; 
	counter++;
      }
    }
  }

  TGraph *graph = new TGraph(counter,x,y);
  graph->Draw("*same");

  TF1 *fun2 = new TF1("fun2","[0]*x +[1]",0,16000);
  graph->Fit("fun2","qROB");
  can->Update();
  //can->WaitPrimitive();

  Double_t gain = fun2->GetParameter(0);
  delete x;
  delete y;
  delete graph;
  delete fun2;

  return gain;
}

void SiRelativeGains_Step1(void)
{
  using namespace std;

  //TFile *f1 = new TFile("/data0/manasta/OrganizeRaw_files/run924_16O_sp7_slope1.root");
  //TFile *f1 = new TFile("/home/lighthall/repository/analysis/ANASEN/anasen_analysis_software/output.root");
  TFile *f1 = new TFile("/home/lighthall/repository/analysis/ANASEN/root/run1209m.root"); 
  if ( !f1->IsOpen() ){
    cout << "Error: Root File Does Not Exist\n";
    exit(EXIT_FAILURE);
  }
  TCanvas *can = new TCanvas("can","can",800,600);

  ofstream outfile;
  outfile.open("QQQRelativeGains_Step1.dat");

  ifstream infile;
  infile.open("QQQRelativeGains09122016_Slope1.dat");
  Int_t det=0,ch=0;
  Double_t slope[4][32];
  Double_t dummy;
  if (infile.is_open()){
    while (!infile.eof()){
      infile >> det >> ch >> dummy;
      slope[det][ch] = dummy;
    }
  }else{
    cout << "Infile not opened\n";
    exit(EXIT_FAILURE);
  }
  infile.close();

  Int_t bad_det[128];
  Int_t bad_front[128];
  Int_t bad_back[128];
  Int_t count_bad = 0;

  for (Int_t DetNum=0; DetNum<4; DetNum++){
    for (Int_t FrontChNum=0; FrontChNum<16; FrontChNum++){
      Int_t BackChNum = 0;
      
      // if(DetNum==1 && (FrontChNum==4 || FrontChNum==14)){continue;}
      //	 if(DetNum==2 && FrontChNum==11){continue;}

      TH2F *hist = NULL;
      TString hname=Form("Q3_back_vs_front%i_%i_%i",DetNum,FrontChNum,BackChNum);
      hist = (TH2F*)f1->Get(hname.Data());
      if (hist==NULL){
	cout << hname << "Histo does not exist\n";
	bad_det[count_bad] = DetNum;
	bad_front[count_bad] = FrontChNum;
	bad_back[count_bad] = BackChNum;
	count_bad++;
	continue;
      }

      Double_t gain = MyFit(hist,can);
      slope[DetNum][FrontChNum+16] = slope[DetNum][FrontChNum+16]*gain;
    }
    for (Int_t i=0; i<32; i++){
      outfile << DetNum << "\t" << i << "\t" << slope[DetNum][i] << endl;
    }
  }
  outfile.close();
  cout << "List of bad detectors:\n";
  for (int i=0; i<count_bad; i++){
    cout << bad_det[i] << "  " << bad_front[i] << "  " << bad_back[i] << endl;
  }
  
  delete can;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////		 

	    
	     
    
