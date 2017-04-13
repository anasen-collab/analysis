//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Relative calibration of Si gains for SX3 Step 1
// See readme.md for general instructions.
// Usage: root -l SiRelativeGains_Step1.C+
//
// Edited by : John Parker , 2016Jan22
// Edited by : Maria Anastasiou, 2016Sept20
// Developed by : Jon Lighthall, 2017.02
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
#include <TLegend.h>

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Double_t MyFit1(TH2F* hist, TCanvas* can) {
  hist->Draw("colz");
  Int_t up=6000;
  hist->GetXaxis()->SetRangeUser(0,up);
  hist->GetYaxis()->SetRangeUser(0,up);

  Bool_t docut=kTRUE;//added to duplicate functionality of some previous versions
  const Int_t nv = 5;
  Double_t x1[nv] = { 450, 5000, 6200, 530, 450 };
  Double_t y1[nv] = { 4800, 10, 1100, 6150, 4800 };
  //Double_t x1[nv] = {190, 240, 240, 13900, 4230, 1700, 165, 190};
  //Double_t y1[nv] = {315, 3215, 13100, 780, 420, 170, 136, 315};
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
  can->WaitPrimitive();

  //cout << fun2->GetChisquare() << "  " << fun2->GetChisquare()/counter;
  Double_t gain = fun2->GetParameter(1);
  delete x;
  delete y;
  delete graph;
  delete fun2;
  
  return gain;
}

Double_t MyFit2(TH2F* hist, TCanvas* can){//manual cut, from FrontFirst directory
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

  for(int n=0;n<cut->GetN();n++){
    cut->GetPoint(n,x1[n],y1[n]);
    cout << x1[n] << "\t" << y1[n] << endl;
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

Double_t MyFit3(TH2F* hist, TCanvas *can){//cut-as-line fit
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

Double_t MyFit4(TH2F* hist, TCanvas* can){//auto calc cut
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
  printf("maxbin is %d, at %f\n",maxbin,hist->GetYaxis()->GetBinCenter(maxbin));
  printf("minbin is %d, at %f\n",minbin,hist->GetXaxis()->GetBinCenter(minbin));

  Double_t slope = -1;
  Double_t offset = hist->GetYaxis()->GetBinCenter(maxbin)+hist->GetXaxis()->GetBinCenter(minbin);

  Int_t steps=1;

  TF1 *fun2;// = new TF1("fun2","[0]*x +[1]",x1,x2);
  for (int k=steps; k>-1; k--) {
    //Set cut shape here; assumes form y=mx+b
    Double_t x1=10; 
    Double_t x2=13000;
    Double_t width=400;
    width*=TMath::Power(2,k);
    printf("width=%f slope=%f offset=%f",width,slope,offset);
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

    printf(" for %s counts = %d in window\n",hist->GetName(),counter);
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
    
  Double_t gain = fun2->GetParameter(1);
  delete fun2;
 
  return gain;
}

void SiRelativeGains_Step1(void)
{
  using namespace std;

  //TFile *f1 = new TFile("run236out_nocal.root");//front
  //TFile *f1 = new TFile("/data0/nabin/ANASEN/ANASEN_NKJ/New/evt2root/run251_NSCL11_Pulser.root");//back
  //TFile *f1 = new TFile("../../../OrganizeRaw_root/run567_051116.root");//front
  //TFile *f1 = new TFile("/data0/manasta/OrganizeRaw_files/run930_931_nospacer_X3slope1_divideback.root"); 
  //TFile *f1 = new TFile("/home/lighthall/anasen/root/run1226-9m.root");
  TFile *f1 = new TFile("/home/lighthall/anasen/root/run1255-61m.root");//all proton scattering
  //TFile *f1 = new TFile("/home/lighthall/anasen/root/run1255-7m.root");
  if ( !f1->IsOpen() ){
    cout << "Error: Root file does not exist\n";
    exit(EXIT_FAILURE);
  }

  //Input the .dat file used by Main.cpp to generate the .root file given above
  ifstream infile;
  infile.open("saves/X3RelativeGains_Slope1.dat");
  Int_t det=0,ch=0;
  Double_t slope[24][12];
  Double_t dummy = 0;
  if (infile.is_open()) {
    infile.ignore(100,'\n');//read in dummy line
    while (!infile.eof()){
      infile >> det >> ch >> dummy;
      slope[det-4][ch] = dummy;
    }
  }else{
    cout << "Error: Dat file does not exist\n";
    exit(EXIT_FAILURE);
  }
  infile.close();

  time_t rawtime;
  struct tm * timeinfo;
  char filename [80];
  time (&rawtime);
  timeinfo = localtime (&rawtime);

  ofstream outfile;
  strftime (filename,80,"saves/X3RelativeGains_Step1_%y%m%d.%H%M%S.dat",timeinfo);
  outfile.open(filename);
  outfile << "DetNum\tFrontCh\tGain\n";

  // ofstream outfile2;
  // strftime (filename,80,"saves/QQQRelativeGains_Step1_%y%m%d.%H%M%S_back.dat",timeinfo);  // file2 may be used for diagnostics
  // outfile2.open(filename);
  // outfile2 << "DetNum\tFrontCh\tBackCh\tOld\t\tSlope\t\tNew\n";

  TCanvas *can = new TCanvas("can","can",800,600);

  Int_t bad_det[288];
  Int_t bad_front[288];
  Int_t bad_back[288];
  Int_t count_bad = 0;

  for (Int_t DetNum=4; DetNum<28; DetNum++) {
    //if(DetNum!=21) continue;
    for (Int_t FrontChNum=0; FrontChNum<4; FrontChNum++) {
      TH2F *hist = NULL;
      //TString hname=Form("down_vs_updivideBack%i_f%i",DetNum,FrontChNum); //normalized
      TString hname=Form("down_vs_up%i_f%i",DetNum,FrontChNum);           //unnormalized
      hist = (TH2F*)f1->Get(hname.Data());
      if (hist==NULL) {
	cout << hname << " histogram does not exist\n";
	bad_det[count_bad] = DetNum;
	bad_front[count_bad] = FrontChNum;
	bad_back[count_bad] = BackChNum;
	count_bad++;
	continue;
      }
      
      Double_t gain = MyFit4(hist,can);
      //Int_t deti=DetNum-4;
      //Int_t frti=FrontChNum+8;
      // printf("Previous gain = %f \t Slope = %f \t New gain = %f\n",slope[DetNum-4][FrontChNum+8],gain, -slope[DetNum][FrontChNum+8]/gain);
      // outfile2 << DetNum << "\t" << FrontChNum+16 << "\t" <<BackChNum << "\t"
      // 	       << left << fixed << setw(8) <<slope[DetNum-4][FrontChNum+8] << "\t"
      // 	       << left << fixed << setw(8) << gain << "\t"
      // 	       << left << fixed << setw(8) << slope[DetNum-4][FrontChNum+8]*gain << endl;
      // slope[DetNum][FrontChNum+16] *= gain;
      // slope[DetNum-4][FrontChNum+8] = -slope[DetNum-4][FrontChNum+8]/average_slope;
    }
    for (Int_t i=0; i<12; i++){
      outfile << DetNum << "\t" << i << "\t" << slope[DetNum-4][i] << endl;
    }
  }
  outfile.close();
  //outfile2.close();
  cout << "List of bad detectors:\n";
  for (Int_t i=0; i<count_bad; i++){
    cout << bad_det[i] << "  " << bad_front[i] << "  " << bad_back[i] << endl;
  }
}
