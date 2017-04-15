//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Fit methods for relative calibration of Si gains for SX3
// See readme.md for general instructions.
//
// Developed by : Jon Lighthall, 2017.04
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//C++
#include <fstream>
#include <exception>
//ROOT
#include <TMath.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TGraph.h>
#include <TF1.h>
#include <TSpectrum.h>
#include <TCutG.h>
#include <TVector.h>
#include <TLegend.h>
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
using namespace std;

class Gains {
 public:
  Double_t old[24][12];
  void Load(TString);
  void Print();
};

class GainMatch {
 public:
  Double_t Fit1(TH2F*,TCanvas*,Bool_t docut=kTRUE);
  Double_t Fit2(TH2F*,TCanvas*);
  Double_t Fit3(TH2F*,TCanvas*);
  Double_t Fit4(TH2F*,TCanvas*);
};

void Gains::Load(TString fname) {
  ifstream infile;
  printf("Loading file %s\n",fname.Data());
  infile.open(fname.Data());
  Int_t det=0,ch=0;
  Double_t dummy = 0;
  if (infile.is_open()) {
    cout << "Read OK"<<endl;
    infile.ignore(100,'\n');//read in dummy line
    while (!infile.eof()){
      infile >> det >> ch >> dummy;
      old[det-4][ch] = dummy;
    }
  }else{
    cout << "Error: Dat file " << fname.Data() << " does not exist\n";
    exit(EXIT_FAILURE);
  }
  infile.close();
}

void Gains::Print() {
  printf("DetNum\tFrontCh\tGain\n");
  for (Int_t i=0; i<24; i++){
    for (Int_t j=0; j<12; j++){
      printf("%d\t%d\t%f\n",i+4,j,old[i][j]);
    }
  }
}

Double_t GainMatch::Fit1(TH2F* hist, TCanvas* can,Bool_t docut) {
  hist->Draw("colz");
  Int_t up=6000;
  hist->GetXaxis()->SetRangeUser(0,up);
  hist->GetYaxis()->SetRangeUser(0,up);
  
  const Int_t nv = 5;
  // Step 1
  Double_t x1[nv] = { 450, 5000, 6200, 530, 450 };
  Double_t y1[nv] = { 4800, 10, 1100, 6150, 4800 };
  //Double_t x1[nv] = {190, 240, 240, 13900, 4230, 1700, 165, 190};
  //Double_t y1[nv] = {315, 3215, 13100, 780, 420, 170, 136, 315};
  // Step 2, 3
  //Double_t x1[nv] = { 2, 12, 9, 0.4, 2 };
  //Double_t y1[nv] = { 0.8, 9.8, 11, 1.1, 0.8 };
    
  TCutG *cut = new TCutG("cut",nv,x1,y1);
  if(docut)cut->Draw("same");

  Int_t counter = 0;
  for (int i=1; i<hist->GetNbinsX(); i++){
    for (int j=1; j<hist->GetNbinsY(); j++){
      if(docut) if ( !cut->IsInside(hist->GetXaxis()->GetBinCenter(i),hist->GetYaxis()->GetBinCenter(j))) {
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
      if(docut)if ( !cut->IsInside(hist->GetXaxis()->GetBinCenter(i),hist->GetYaxis()->GetBinCenter(j))) {
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

  TF1 *fun2 = new TF1("fun2","[0]+[1]*x",0,16000);
  fun2->SetParameter(0,10);
  fun2->SetParameter(1,-1);
  graph->Fit("fun2","qROB");
  can->Update();
  
  //cout << fun2->GetChisquare() << "  " << fun2->GetChisquare()/counter;
  Double_t gain = fun2->GetParameter(1);
  delete x;
  delete y;
  delete graph;
  delete fun2;
  
  return gain;
}

Double_t GainMatch::Fit2(TH2F* hist, TCanvas* can){//manual cut
  if(!(can->GetShowEventStatus()))can->ToggleEventStatus();
  if(!(can->GetShowToolBar()))can->ToggleToolBar();
  hist->Draw("colz");
  Int_t up=6000;
  hist->GetXaxis()->SetRangeUser(0,up);
  hist->GetYaxis()->SetRangeUser(0,up);
    
  vector<double> x1;
  vector<double> y1;
  
  TCutG *cut;
  cut = (TCutG*)can->WaitPrimitive("CUTG");
  x1.resize(cut->GetN());
  y1.resize(cut->GetN());

  printf("Verticies of cut are:\n");
    for(int n=0;n<cut->GetN();n++){
      cut->GetPoint(n,x1[n],y1[n]);
      cout << "\t" << x1[n] << "\t" << y1[n] << endl;
  }
  cut->SetLineColor(6);
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

  TF1 *fun2 = new TF1("fun2","[0]+[1]*x",0,6000);
  fun2->SetParameter(0,10);
  fun2->SetParameter(1,-1);
  graph->Fit("fun2","qROB");
  can->Update();
  can->WaitPrimitive();

  //cout << fun2->GetChisquare() << "  " << fun2->GetChisquare()/counter;
  Double_t gain = fun2->GetParameter(1);
  delete x;
  delete y;
  delete graph;
  delete fun2;

  return gain;
}

Double_t GainMatch::Fit3(TH2F* hist, TCanvas *can){//cut-as-line fit
  if(!(can->GetShowEventStatus()))can->ToggleEventStatus();
  if(!(can->GetShowToolBar()))can->ToggleToolBar();
  hist->Draw("colz");
  Int_t up=6000;
  hist->GetXaxis()->SetRangeUser(0,up);
  hist->GetYaxis()->SetRangeUser(0,up);

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

Double_t GainMatch::Fit4(TH2F* hist, TCanvas* can){//auto calc cut
  can->Clear();
  hist->Draw("colz");
  Int_t up=6000;
  hist->GetXaxis()->SetRangeUser(0,up);
  hist->GetYaxis()->SetRangeUser(0,up);

  Int_t maxbin=0;
  Int_t minbin=hist->GetNbinsX();
  for (int i=1; i<hist->GetNbinsX(); i++) {
    for (int j=1; j<hist->GetNbinsY(); j++){
      if(hist->GetBinContent(i,j)>0) {
	if(j>maxbin) {
	  maxbin=j;
	  minbin=i;
	}
      }
    }
  }
  printf(" for %s ",hist->GetName());
  printf("maxbin is %4d, at %6g ",maxbin,hist->GetYaxis()->GetBinCenter(maxbin));
  printf("minbin is %4d, at %6g\n",minbin,hist->GetXaxis()->GetBinCenter(minbin));

  Double_t slope = -1;
  Double_t offset = hist->GetYaxis()->GetBinCenter(maxbin)+hist->GetXaxis()->GetBinCenter(minbin);

  Int_t steps=2;

  TF1 *fun2;// = new TF1("fun2","[0]*x +[1]",x1,x2);
  for (int k=steps; k>-1; k--) {
    //Set cut shape here; assumes form y=mx+b
    Double_t x1=10; 
    Double_t x2=13000;
    Double_t width=400;
    width*=TMath::Power(2,k);
    printf(" Step %d: width=%5.0f slope=%9.5f offset=%7.2f",steps-k+1,width,slope,offset);
    //The corners of a parallelepiped cut window are then calculated
    Double_t y1=slope*x1+offset;
    Double_t y2=slope*x2+offset;
    Double_t dx=(width/2.0)*slope/sqrt(1+slope*slope);
    Double_t dy=(width/2.0)*1/sqrt(1+slope*slope);
    const Int_t nv = 5;//set number of verticies
    Double_t xc[nv] = {x1+dx,x2+dx,x2-dx,x1-dx,x1+dx};
    Double_t yc[nv] = {y1-dy,y2-dy,y2+dy,y1+dy,y1-dy};
    TCutG *cut = new TCutG("cut",nv,xc,yc);
    cut->SetLineColor(6);
    cut->Draw("same");
  
    Int_t counter = 0;
    for (int i=1; i<hist->GetNbinsX(); i++) {//determine number of counts inside window
      for (int j=1; j<hist->GetNbinsY(); j++) {
	if ( !cut->IsInside(hist->GetXaxis()->GetBinCenter(i),hist->GetYaxis()->GetBinCenter(j))) {
	  continue;
	}
	counter+=(Int_t)hist->GetBinContent(i,j);
      }
    }

    printf(" counts = %d in window\n",counter);
    Double_t *x = new Double_t[counter];
    Double_t *y = new Double_t[counter];

    counter = 0;
    for (int i=1; i<hist->GetNbinsX(); i++){//fill vectors with histogram entries
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
    //graph->SetMarkerSize();
    graph->Draw("same*");
    fun2 = new TF1("fun2","[0]*x +[1]",x1,x2);
    //fun2->SetLineWidth(1);
    fun2->SetLineColor(k+2);
    graph->Fit("fun2","qROB");
  
    fun2->Draw("same");
  
    TLegend *leg = new TLegend(0.1,0.75,0.2,0.9);
   
    leg->AddEntry(cut,Form("cut %.0f wide",width),"l");
    //leg->AddEntry(graph,"graph","p");
    leg->AddEntry(fun2,Form("TGraph fit #%d",steps-k+1),"l");
    leg->Draw();
  
    can->Update();
    //if(k==0) can->WaitPrimitive();
    slope=fun2->GetParameter(0);
    offset=fun2->GetParameter(1);
  
    delete x;
    delete y;
    delete graph;
  }
    
  Double_t gain = slope;
  delete fun2;
 
  return gain;
}
