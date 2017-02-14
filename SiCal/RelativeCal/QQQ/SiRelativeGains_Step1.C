//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////Relative calibration of Si gains for QQQ
////Essentially the same progam as that for the SX3
////root -l SiRelativeGains_Step1.C+


//// Edited by : John Parker , 2016Jan22
//   Developed by : Jon Lighthall, 2016.12
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
#include <TLegend.h>

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
TFile *f1;

Double_t MyFit1(TH2F* hist, TCanvas *can) {
  //Method 1 - calculates slope of points wihtin pre-defined cut using TGraph
  hist->Draw("colz");
  
  Double_t x1[12] = {1450, 630, 3150, 6340, 9200, 10540, 13200, 13600, 11670, 7550, 2400, 1450};
  Double_t y1[12] = {800, 2250, 4900, 8050, 10380, 11520, 13800, 12600, 8700, 5400, 1100, 800};
  TCutG *cut = new TCutG("cut",12,x1,y1);
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

  TF1 *fun2 = new TF1("fun2","[0]*x +[1]",0,16000);
  //fun2->SetLineWidth(1);
  //fun2->SetLineColor(2);
  graph->Fit("fun2","qROB");

  TLegend *leg = new TLegend(0.1,0.75,0.2,0.9);
  leg->AddEntry(cut,"cut","l");
  leg->AddEntry(graph,"graph","p");
  leg->AddEntry(fun2,"TGraph fit","l");
  leg->Draw();
  
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
  if(!(can->GetShowEventStatus()))can->ToggleEventStatus();
  if(!(can->GetShowToolBar()))can->ToggleToolBar();
  hist->Draw("colz");

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

Double_t MyFit4(TH2F* hist, TCanvas *can) {
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
  TF1 *fun3 = new TF1("fun3","[0]*x +[1]",0,16384);
  fun3->SetLineColor(4);
  fun3->SetLineStyle(2);
  fun3->SetLineWidth(2);
  xprof->Fit("fun3","q");
  Double_t slope = fun3->GetParameter(0);
  Double_t offset = fun3->GetParameter(1);

  Int_t steps=3;

  TF1 *fun2;// = new TF1("fun2","[0]*x +[1]",x1,x2);
  for (int k=steps; k>-1; k--) {
    //Set cut shape here; assumes form y=mx+b
    Double_t x1=500; 
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

Double_t MyFit6(TH2F* hist, TCanvas *can) {//used for Det 2; using fixed initial slope
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
  TF1 *fun3 = new TF1("fun3","[0]*x +[1]",0,16384);
  fun3->SetLineColor(4);
  fun3->SetLineStyle(2);
  fun3->SetLineWidth(2);
  fun3->FixParameter(0,1);
  fun3->FixParameter(1,0);
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
  printf("high bin is %d %.0f: ",highb,high);
					
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
    
    printf("width=%f slope=%f offset=%f",width,slope,offset);
    //The corners of a truncated cone cut window are then calculated
    Double_t y1=slope*x1+offset;
    Double_t y2=slope*x2+offset;
    Double_t dx0=(width0/2.0)*slope/sqrt(1+slope*slope);
    Double_t dy0=(width0/2.0)*1/sqrt(1+slope*slope);
    
    width*=TMath::Power(2,k);
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

Double_t MyFit7(TH2F* hist, TCanvas *can, Int_t DetNum, Int_t FrontChNum, Int_t BackChNum) {//used for Det 2; using fixed initial slope
  //Method 7 - automated cut based on TGraph
  can->Clear();
  hist->Draw("colz");
  Int_t up=6000;
  hist->GetXaxis()->SetRangeUser(0,up);
  hist->GetYaxis()->SetRangeUser(0,up);

  printf("det = %d, fr = %d bk = %d \n", DetNum,FrontChNum, BackChNum );

  TTree *intree = (TTree*) f1->Get("MainTree");
  
  intree->Draw("Si.Detector.EBack_Pulser:Si.Detector.EFront_Pulser","Si.Hit.HitType==11&&Si.Detector.DetID==0&&Si.Detector.FrontChNum==0&&Si.Detector.BackChNum==0","");
  
  /* TGraph *graph = new TGraph(intree->GetSelectedRows(),intree->GetV1(),intree->GetV2());
  TF1* fun2 = new TF1("fun2","[0]*x +[1]",0,up);
  graph->Fit("fun2","qROB");
  //fun2->Draw("same");
  Double_t slope=fun2->GetParameter(0);
  
  //  delete graph;
  Double_t gain = slope;
  //delete fun2;
  */
  Double_t gain = 1;
  return gain;
}

void SiRelativeGains_Step1(void)
{
  using namespace std;

  //TFile *f1 = new TFile("/data0/manasta/OrganizeRaw_files/run924_16O_sp7_slope1.root");
  f1 = new TFile("/home/lighthall/anasen/root/run1226-9m.root");
  if ( !f1->IsOpen() ){
    cout << "Error: Root File Does Not Exist\n";
    exit(EXIT_FAILURE);
  }
  TCanvas *can = new TCanvas("can","can",1362,656);
  can->SetWindowPosition(0,63);

  ofstream outfile;
  outfile.open("QQQRelativeGains_Step1_all.dat");

  ofstream outfile2;
  outfile2.open("QQQRelativeGains_Step1_list.dat");

  ifstream infile;
  infile.open("QQQRelativeGains_Slope1.dat");
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

  for (Int_t DetNum=0; DetNum<1; DetNum++) {
    for (Int_t FrontChNum=0; FrontChNum<1; FrontChNum++) {
      //for (Int_t BackChNum=0; BackChNum<16; BackChNum++) {//loop over back (diagnostic)
	Int_t BackChNum = 0;
	if(DetNum==2)
	  BackChNum = 4;
	//if(DetNum==0 && FrontChNum==0){continue;}
	//	 if(DetNum==2 && FrontChNum==11){continue;}
	//if(!((FrontChNum==0)||(FrontChNum==13)))
	//if(!((FrontChNum==13)))
	//continue;
	printf("det= %i front = %i  back = %i\n",DetNum,FrontChNum,BackChNum);
	TH2F *hist = NULL;
	TString hname=Form("Q3_back_vs_front%i_%i_%i",DetNum,FrontChNum,BackChNum);
	hist = (TH2F*)f1->Get(hname.Data());
	if (hist==NULL) {
	  cout << hname << " histogram does not exist\n";
	  bad_det[count_bad] = DetNum;
	  bad_front[count_bad] = FrontChNum;
	  bad_back[count_bad] = BackChNum;
	  count_bad++;
	  continue;
	}

	//Double_t gain = MyFit6(hist,can); //set fit method here
	Double_t gain = MyFit7(hist,can,DetNum,FrontChNum,BackChNum); //set fit method here
	slope[DetNum][FrontChNum+16] = slope[DetNum][FrontChNum+16]*gain;
	outfile2 << DetNum << "\t" << FrontChNum << "\t" <<BackChNum << "\t" << gain << endl;
	//}//back loop
    }
    for (Int_t i=0; i<32; i++){
      outfile << DetNum << "\t" << i << "\t" << slope[DetNum][i] << endl;
    }
  }
  outfile.close();
  cout << "List of bad detectors:\n";
  for (Int_t i=0; i<count_bad; i++){
    cout << bad_det[i] << "  " << bad_front[i] << "  " << bad_back[i] << endl;
  }
  
  //delete can;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////		     
