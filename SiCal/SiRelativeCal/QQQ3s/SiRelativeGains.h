//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Fit methods for relative calibration of Si gains for QQQ
// See readme.md for general instructions.
//
// Developed by : Jon Lighthall, 2017.01
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//C++
#include <fstream>
#include <exception>
#include <time.h>
#include <iomanip>  //file formatting; std::setfill, std::setw
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
#include <TProfile.h>
#include <TROOT.h>
#include <TLegend.h>
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
using namespace std;
const Int_t ndets=4;
const Int_t nchan=32;
const Int_t range=4096*4;

class Gains {
 public:
  Double_t old[ndets][nchan];
  void Load(TString);
  void Print();
};

class GainMatch {
 public:
  Double_t Fit1(TH2F*,TCanvas*);
  Double_t Fit2(TH2F*,TCanvas*);
  Double_t Fit4(TH2F*,TCanvas*);
  Double_t Fit6(TH2F*,TCanvas*);
};

void Gains::Load(TString fname) {
  ifstream infile;
  printf("Loading file %s...",fname.Data());
  infile.open(fname.Data());
  Int_t det=0,ch=0;
  Double_t dummy = 0;
  if (infile.is_open()) {
    cout << "Read OK"<<endl;
    infile.ignore(100,'\n');//read in dummy line
    while (!infile.eof()){
      infile >> det >> ch >> dummy;
      old[det][ch] = dummy;
    }
  }else{
    cout << "Error: Dat file " << fname.Data() << " does not exist\n";
    exit(EXIT_FAILURE);
  }
  infile.close();
}

void Gains::Print() {
  printf("DetNum\tFrontCh\tGain\n");
  for (Int_t i=0; i<ndets; i++){
    for (Int_t j=0; j<nchan; j++){
      printf("%d\t%d\t%f\n",i,j,old[i][j]);
    }
  }
}

Double_t GainMatch::Fit1(TH2F* hist, TCanvas *can) {
  //Method 1 - calculates slope of points wihtin pre-defined cut using TGraph
  hist->Draw("colz");

  const Int_t nv = 12;
  
  //Double_t x1[nv] = { 2, 12, 9, 0.4, 2 };
  //Double_t y1[nv] = { 0.8, 9.8, 11, 1.1, 0.8 };
  //Double_t x1[nv] = {1168, 450, 4490, 13160, 13440, 7650, 3790, 1116, 1168};//tight cut (exclude random points which are problem for some plots)
  //Double_t y1[nv] = {700, 1600, 6400, 15320, 13450, 6750, 2965, 780,700};
  Double_t x1[nv] = {1450, 630, 3150, 6340, 9200, 10540, 13200, 13600, 11670, 7550, 2400, 1450};
  Double_t y1[nv] = {800, 2250, 4900, 8050, 10380, 11520, 13800, 12600, 8700, 5400, 1100, 800};
  //Double_t x1[nv] = {1500, 6360, 10690, 14210, 11228, 4670, 800, 730, 940, 1500 }; //a bit broader cut
  //Double_t y1[nv] = {315, 800, 5500, 11600, 15675, 10900, 4220, 1425, 570, 315};
  //Double_t x1[nv] = {1000, 10000, 14000, 14000, 11500, 1470, 500, 90, 90, 1000}; //broadest cut
  //Double_t y1[nv] = {100, 100, 100, 11800, 15170, 13400, 13130, 8260, 920, 100}; 
  TCutG *cut = new TCutG("cut",nv,x1,y1);
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

  TF1 *fun2 = new TF1("fun2","[0]*x +[1]",0,range);
  graph->Fit("fun2","qROB");

  TLegend *leg = new TLegend(0.1,0.75,0.2,0.9);
  leg->AddEntry(cut,"cut","l");
  leg->AddEntry(graph,"graph","p");
  leg->AddEntry(fun2,"TGraph fit","l");
  leg->Draw();
  
  can->Update();
 
  Double_t gain = fun2->GetParameter(0);
  delete x;
  delete y;
  delete graph;
  delete fun2;
  
  return gain;
}

Double_t GainMatch::Fit2(TH2F* hist, TCanvas *can) {
  //Method 2 - calculates slope of points wihtin user-defined cut using TGraph
  if(!(can->GetShowEventStatus()))can->ToggleEventStatus();
  if(!(can->GetShowToolBar()))can->ToggleToolBar();
  hist->Draw("colz");

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

  TF1 *fun2 = new TF1("fun2","[0]*x +[1]",0,range);
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

Double_t GainMatch::Fit4(TH2F* hist, TCanvas *can) {
  //Method 4 - automated cut generation based on TProfile slope
  can->Clear();
  hist->Draw("colz");
  Int_t up=6000;
  hist->GetXaxis()->SetRangeUser(0,up);
  hist->GetYaxis()->SetRangeUser(0,up);
  hist->ProfileX();
  TString hname;
  hname=hist->GetName();
  hname+="_pfx";
  TProfile *xprof=(TProfile *)gROOT->FindObject(hname.Data());
  xprof->SetMarkerStyle(4);
  xprof->SetMarkerSize(0.125);
  xprof->Draw("same");
  TF1 *fun3 = new TF1("fun3","[0]*x +[1]",0,range);
  fun3->SetLineColor(4);
  fun3->SetLineStyle(2);
  fun3->SetLineWidth(2);
  xprof->Fit("fun3","q");
  Double_t slope = fun3->GetParameter(0);
  Double_t offset = fun3->GetParameter(1);


  printf(" for %s ",hist->GetName());
  Int_t steps=3;

  TF1 *fun2;// = new TF1("fun2","[0]*x +[1]",x1,x2);
  for (int k=steps; k>-1; k--) {
    //Set cut shape here; assumes form y=mx+b
    Double_t x1=500; 
    Double_t x2=13000;
    Double_t width=400;
    width*=TMath::Power(2,k);
    printf(" Step: %d width=%5.0f slope=%9.5f offset=%7.2f",steps-k+1,width,slope,offset);
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
    //graph->Draw("same*");
    fun2 = new TF1("fun2","[0]*x +[1]",x1,x2);
    //fun2->SetLineWidth(1);
    fun2->SetLineColor(k+2);
    graph->Fit("fun2","qROB");
  
    fun2->Draw("same");
    fun3->Draw("same");
  
    TLegend *leg = new TLegend(0.1,0.75,0.2,0.9);
    leg->AddEntry(xprof,"x-profile","pe");  
    leg->AddEntry(fun3,"TProfile fit","l");   
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
  delete xprof;
  delete fun2;
  delete fun3;

  return gain;
}

Double_t GainMatch::Fit6(TH2F* hist, TCanvas *can) {//used for Det 2; using fixed initial slope
  //Method 6 - automated cut generation based on TProfile slope; cone-shaped
  can->Clear();
  hist->Draw("colz");
  Int_t up=6000;
  hist->GetXaxis()->SetRangeUser(0,up);
  hist->GetYaxis()->SetRangeUser(0,up);
  hist->ProjectionX();
  TString hname;
  hname=hist->GetName();
  hname+="_px";
  TH1 *xproj=(TH1 *)gROOT->FindObject(hname.Data());
  TF1 *fun3 = new TF1("fun3","[0]*x +[1]",0,range);
  fun3->SetLineColor(4);
  fun3->SetLineStyle(2);
  fun3->SetLineWidth(2);
  fun3->FixParameter(0,1); //default slope
  fun3->FixParameter(1,0); //default offset
  Double_t slope = fun3->GetParameter(0);
  Double_t offset = fun3->GetParameter(1);

  Int_t highb=0;
  for(Int_t i=0; i<xproj->GetXaxis()->GetNbins();i++) {
    Double_t cont=xproj->GetBinContent(i);
    if(cont>0)
      if(i>highb)
	highb=i;
  }
  Double_t high=xproj->GetBinCenter(highb);
  printf("For %s: ",hist->GetName());
  printf("high bin is %d %.0f\n",highb,high);
					
  Int_t steps=4;
  TF1 *fun2;// = new TF1("fun2","[0]*x +[1]",x1,x2);
  for (int k=steps; k>-1; k--) {
    //Set cut shape here; assumes form y=mx+b
    //Set variable low-x position
    Double_t xa=400;
    Double_t xb=1500;
    Double_t dxs=(xb-xa)/steps;
    Double_t x1=xa+(dxs*k); 

    //Set high-x and width
    Double_t x2=1.5*high;
    Double_t width0=200;
    Double_t width=width0;
    
    //The corners of a truncated cone cut window are then calculated
    Double_t y1=slope*x1+offset;
    Double_t y2=slope*x2+offset;
    Double_t dx0=(width0/2.0)*slope/sqrt(1+slope*slope);
    Double_t dy0=(width0/2.0)*1/sqrt(1+slope*slope);
    
    width*=TMath::Power(2,k);
    printf(" Step %d: low=%4.0f width=%5.0f slope=%9.5f offset=%7.2f",steps-k+1,x1,width,slope,offset);
    Double_t dx=(width/2.0)*slope/sqrt(1+slope*slope);
    Double_t dy=(width/2.0)*1/sqrt(1+slope*slope);
    const Int_t nv = 5;//set number of verticies
    Double_t xc[nv] = {x1+dx0,x2+dx,x2-dx,x1-dx0,x1+dx0};
    Double_t yc[nv] = {y1-dy0,y2-dy,y2+dy,y1+dy0,y1-dy0};
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
    //graph->Draw("same*");
    fun2 = new TF1("fun2","[0]*x +[1]",x1,x2);
    //fun2->SetLineWidth(1);
    fun2->SetLineColor(k+2);
    graph->Fit("fun2","qROB");
  
    fun2->Draw("same");
    fun3->Draw("same");
  
    TLegend *leg = new TLegend(0.1,0.75,0.2,0.9);

    leg->AddEntry(fun3,"TProfile fit","l");   
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
  delete fun3;

  return gain;
}
