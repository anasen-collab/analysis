//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////Relative calibration of Si gains for QQQ
////Essentially the same progam as that for the SX3
////root -l SiRelativeGains_Step1.C+


//// Edited by : John Parker , 2016Jan22
//   Developed by : Jon Lighthall, 2016.01
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
#include <TROOT.h>
#include <TProfile.h>

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Double_t MyFit1(TH2F* hist, TCanvas *can) {
  //Method 1 - calculates slope of points wihtin pre-defined cut using TGraph
  hist->Draw("colz");
  
  Double_t x1[12] = {1450, 630, 3150, 6340, 9200, 10540, 13200, 13600, 11670, 7550, 2400, 1450};
  Double_t y1[12] = {800, 2250, 4900, 8050, 10380, 11520, 13800, 12600, 8700, 5400, 1100, 800};
  TCutG *cut = new TCutG("cut",12,x1,y1);
  cut->Draw("same");
  
  Int_t counter = 0;
  for (int i=1; i<hist->GetNbinsX(); i++){
    for (int j=1; j<hist->GetNbinsY(); j++){
      if ( !cut->IsInside(hist->GetXaxis()->GetBinCenter(i),hist->GetYaxis()->GetBinCenter(j))) {
	continue;
      }
      counter+=(Int_t)hist->GetBinContent(i,j);
    }
  }

  Double_t *x = new Double_t[counter];
  Double_t *y = new Double_t[counter];
  
  counter = 0;
  for (int i=1; i<hist->GetNbinsX(); i++){
    for (int j=1; j<hist->GetNbinsY(); j++){
      if ( !cut->IsInside(hist->GetXaxis()->GetBinCenter(i),hist->GetYaxis()->GetBinCenter(j))) {
	continue;
      }
      for (int k=0; k<hist->GetBinContent(i,j); k++){
  	x[counter] = hist->GetXaxis()->GetBinCenter(i);
  	y[counter] = hist->GetYaxis()->GetBinCenter(j);
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

Double_t MyFit2(TH2F* hist, TCanvas *can) {
  //Method 2 - calculates slope of points wihtin user-defined cut using TGraph
  // This method works very similar to method 1, except that instead of a predefined cut, the user must
  // draw their own graphical cut. To do this, in the canvas: View->Toolbar and click on the scissors on the top
  // right hand corner. Then, select the region around your data by clicking. Double click to close the cut.
  // A best fit line should appear through your data
  hist->Draw("colz");
  
  TCutG *cut;
  cut = (TCutG*)can->WaitPrimitive("CUTG");
  x1.resize(cut->GetN());
  y1.resize(cut->GetN());

  for(int n=0;n<cut->GetN();n++){
    cut->GetPoint(n,x1[n],y1[n]);
    cout << x1[n] << "\t" << y1[n] << endl;
  }
  
  cut->Draw("same");
  
  Int_t counter = 0;
  for (int i=1; i<hist->GetNbinsX(); i++){
    for (int j=1; j<hist->GetNbinsY(); j++){
      if ( !cut->IsInside(hist->GetXaxis()->GetBinCenter(i),hist->GetYaxis()->GetBinCenter(j))) {
	continue;
      }
      counter+=(Int_t)hist->GetBinContent(i,j);
    }
  }

  Double_t *x = new Double_t[counter];
  Double_t *y = new Double_t[counter];
  
  counter = 0;
  for (int i=1; i<hist->GetNbinsX(); i++){
    for (int j=1; j<hist->GetNbinsY(); j++){
      if ( !cut->IsInside(hist->GetXaxis()->GetBinCenter(i),hist->GetYaxis()->GetBinCenter(j))) {
	continue;
      }
      for (int k=0; k<hist->GetBinContent(i,j); k++){
  	x[counter] = hist->GetXaxis()->GetBinCenter(i);
  	y[counter] = hist->GetYaxis()->GetBinCenter(j);
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

Double_t MyFit3(TH2F* hist, TCanvas *can) {
  //Expanded Method 1 - three slopes are calculated:
  // 1) slope from un-gated fit of 2D histogram
  // 2) slope from gated fit of TGraph (same as Method 1)
  // 3) slope from un-gated fit of profile of 2D histogram (uses the average y-position for each x-position)
  can->Clear();

  TF1 *fun1 = new TF1("fun1","[0]*x +[1]",0,16384);
  fun1->SetLineColor(3);
  fun1->SetLineStyle(2);
  fun1->SetLineWidth(1);
  hist->Fit("fun1","q");

  Double_t x1[12] = {1450, 630, 3150, 6340, 9200, 10540, 13200, 13600, 11670, 7550, 2400, 1450};
  Double_t y1[12] = {800, 2250, 4900, 8050, 10380, 11520, 13800, 12600, 8700, 5400, 1100, 800};
  TCutG *cut = new TCutG("cut",12,x1,y1);
  
  Int_t counter = 0;
  for (int i=1; i<hist->GetNbinsX(); i++) {
    for (int j=1; j<hist->GetNbinsY(); j++) {
      if ( !cut->IsInside(hist->GetXaxis()->GetBinCenter(i),hist->GetYaxis()->GetBinCenter(j))) {
	continue;
      }
      counter+=(Int_t)hist->GetBinContent(i,j);
    }
  }

  printf(" for %s counts = %d\n",hist->GetName(),counter);
  Double_t *x = new Double_t[counter];
  Double_t *y = new Double_t[counter];

  counter = 0;
  for (int i=1; i<hist->GetNbinsX(); i++){
    for (int j=1; j<hist->GetNbinsY(); j++){
      if ( !cut->IsInside(hist->GetXaxis()->GetBinCenter(i),hist->GetYaxis()->GetBinCenter(j))) {
	continue;
      }
      for (int k=0; k<hist->GetBinContent(i,j); k++){
  	x[counter] = hist->GetXaxis()->GetBinCenter(i);
  	y[counter] = hist->GetYaxis()->GetBinCenter(j);
  	counter++;
      }
    }
  }

  TGraph *graph = new TGraph(counter,x,y);
  TF1 *fun2 = new TF1("fun2","[0]*x +[1]",0,16000);
  fun2->SetLineWidth(1);
  fun2->SetLineColor(1);
  graph->Fit("fun2","qROB");
  
  hist->ProfileX();
  TString hname;
  hname=hist->GetName();
  hname+="_pfx";
  TProfile *xprof=(TProfile *)gROOT->FindObject(hname.Data());
  xprof->SetMarkerStyle(4);
  xprof->SetMarkerSize(0.5);
  TF1 *fun3 = new TF1("fun3","[0]*x +[1]",0,16384);
  fun3->SetLineColor(2);
  fun3->SetLineStyle(1);
  fun3->SetLineWidth(1);
  
  hist->Draw("colz");
  cut->Draw("same");
  graph->Draw("*same");
  fun2->Draw("same");
  xprof->Draw("same");
  xprof->Fit("fun3","q");

  can->Update();
  can->WaitPrimitive();

  Double_t gain = fun2->GetParameter(0);
  Double_t gain2 = fun3->GetParameter(0);
  printf(" gain ratio is %f\n",fabs(gain-gain2)/gain);
  delete x;
  delete y;
  delete graph;
  delete xprof;
  delete fun1;
  delete fun2;
  delete fun3;

  return gain2;
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
	cout << hname << " histogram does not exist\n";
	bad_det[count_bad] = DetNum;
	bad_front[count_bad] = FrontChNum;
	bad_back[count_bad] = BackChNum;
	count_bad++;
	continue;
      }

      Double_t gain = MyFit1b(hist,can);
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

	    
	     
    
