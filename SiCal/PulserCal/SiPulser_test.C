//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ASICS all channels Pulser calibration
// See readme.md for general instructions.
//
// Edited by : Nabin Rijal , 2015Dec13
//             Jon Lighthall Nov 2016
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//C++
#include <fstream>
//ROOT
#include <TMath.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1I.h>
#include <TGraph.h>
#include <TF1.h>
#include <TSpectrum.h>

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void SiPulser_test (void) {
  const Int_t run = 1262;

  if(run==1262||run==1264) {
    const Int_t npeaks = 5;
    Float_t Volts[npeaks] = { 0.5, 1.5, 3.0, 7.0, 9.0};
    if(run==1262) {
      TFile *f1 = new TFile("/data0/lighthall/root/main/run1262.root");
      TString lname="EBack_Pulser[0]";
        }
    if(run==1264) {
      TFile *f1 = new TFile("/data0/lighthall/root/main/run1264.root");
      TString lname="EFront_Pulser[0]";
    }
  }
  if ( !f1->IsOpen() ) {
    cout << "Error: Root file does not exist\n";
    exit(EXIT_FAILURE);
  }

  const int npar=npeaks*3;
  
  Float_t Volts1[npeaks-1];//(n-1) peaks, replaces Volts5
  for (Int_t i=0;i<(npeaks-1);i++) {
    Volts1[npeaks-1-i-1]=Volts[npeaks-i-1];
  }

  Float_t Volts2[npeaks-2];//(n-2 peaks) replaces Volts4
  for (Int_t i=0;i<(npeaks-2);i++) {
    Volts2[npeaks-2-i-1]=Volts[npeaks-i-1];
    }

  Float_t Volts3[npeaks-3];//(n-3 peaks) 
  for (Int_t i=0;i<(npeaks-3);i++) {
    Volts3[npeaks-3-i-1]=Volts[npeaks-i-1];
    }

  Float_t Volts4[npeaks-4];//(n-4 peaks) 
  for (Int_t i=0;i<(npeaks-4);i++) {
    Volts4[npeaks-4-i-1]=Volts[npeaks-i-1];
    }

  Int_t cwide=1362;
  TCanvas *c1 = new TCanvas();
  c1->SetWindowPosition(0,63);
  c1->SetWindowSize(cwide,656);
  c1->Divide(1,2);

  Bool_t dowait=kFALSE; //wait betweeen fits
  Bool_t dotest=kFALSE; //test results
  
  if(dowait) {
    TCanvas *c2 = new TCanvas("c2","double-click me",270,100);
    c2->cd();
    TButton *but2 = new TButton("quit ROOT",".q",.05,.1,.45,.45);
    but2->Draw();
  }

  if(dotest) {
    TCanvas *c3 = new TCanvas("c3","test fit");
    c3->SetWindowPosition(c1->GetWindowTopX()+c1->GetWindowWidth(),c1->GetWindowTopY()-22);      
    c3->SetWindowSize(cwide,656);
    c3->Divide(1,2);
  }

  Int_t  size = 4096*4;
  TH1I *h1 = new TH1I("h1","h1",size,0,size);//min was 300
  TTree *MainTree = NULL;
  MainTree = (TTree*)f1->Get("MainTree");
  if (MainTree==NULL) {
    cout << "Tree does not exist\n";
    exit(EXIT_FAILURE);
  }
  if(dotest) {
    TH1I *h3 = new TH1I("h3","h3",size,0,size);
  }
  
  TF1 *fit = new TF1("fit","pol1",0,size/2);
  
  Double_t zeroshift = 0;//zero voltage shift (fit offset)
  Double_t vperch = 0;//volts per channel (fit slope)
  Double_t q0 = 0;
  TSpectrum *s = 0;
  TGraph *FitGraph = 0;
  
  ofstream outfile; //offsets, ROB
  
  outfile.open("saves/Sipulser_offsets.dat");
  
  if(dotest)
    TGraph *FitGraph3 = 0;
    
  TPolyMarker *pm2; //calculated intercepts
  
  for (Int_t DetNum=0; DetNum<28; DetNum++) {
    c1->cd(1);
    MainTree->Draw(Form("%s >> h1",lname.Data()),Form("DetID==%i",DetNum),"");
    h1->SetTitle(Form("Det%i",DetNum));
    
    c1->Update();
    //c1->WaitPrimitive();
    if(dowait) {
      c2->Update();      
      //c2->WaitPrimitive();
    }
    if(h1->GetEntries()==0) continue;
    
    if(s!=0) {
      delete s;
    }
    
    //h1->GetXaxis()->UnZoom();
    h1->GetXaxis()->SetRangeUser(100,16000);
    TSpectrum *s = new TSpectrum();
    //s->SetResolution(3);
    Float_t ssigma=10;
    Float_t sthresh=0.15;
    Int_t nfound = s->Search(h1,ssigma," ",sthresh);//run 1264
    
    //sort peaks in order of channle number
    Float_t *xpeaks = s->GetPositionX();
    Float_t Temp=0;
    for(Int_t i=0;i<nfound;i++) {
      for(Int_t j=i;j<nfound;j++) {
	if (xpeaks[j] < xpeaks[i]) {
	  Temp = xpeaks[i];
	  xpeaks[i] = xpeaks[j];
	  xpeaks[j] = Temp;  
	      }
      }
    }
    
    c1->cd(2);
	  if(FitGraph!=0) {
	    delete FitGraph;
	  }
	  
	  if(nfound==npeaks) 
	    FitGraph = new TGraph(nfound,xpeaks, &(Volts[0]));
	  else if(nfound==(npeaks-1)) 
	    FitGraph = new TGraph(nfound,xpeaks, &(Volts1[0]));
	  else if(nfound==(npeaks-2)) 
	    FitGraph = new TGraph(nfound,xpeaks, &(Volts2[0]));
	  else if(nfound==(npeaks-3)) 
	    FitGraph = new TGraph(nfound,xpeaks, &(Volts3[0]));
	  else if(nfound==(npeaks-4)) 
	    FitGraph = new TGraph(nfound,xpeaks, &(Volts4[0]));
	  else {
	    printf("number of peaks found = %d of %d\n", nfound, npeaks);
	    
	    continue;
	  }
		
	  pm2=new TPolyMarker(4);
	  pm2->SetMarkerStyle(2);
	  pm2->SetMarkerSize(1.2);
	  pm2->SetMarkerColor(6);
	       
	  FitGraph->GetHistogram()->GetXaxis()->SetTitle("Peak positions (channel)");
	  FitGraph->GetHistogram()->GetYaxis()->SetTitle("Pulser setting (voltage)");
    	  FitGraph->Draw("AP*");
	  if(nfound>3)
	    FitGraph->Fit("fit","qROB=0.95");
	  else
	    FitGraph->Fit("fit","qROB=1");

	  leg = new TLegend(0.1,0.75,0.2,0.9);
	  leg->AddEntry(fit,"pol1, ROB=95","l");
	  leg->Draw();
		
	  zeroshift = fit->GetParameter(0);
	  vperch = fit->GetParameter(1);
	  q0 = -zeroshift/vperch;
	  pm2->SetPoint(0,q0,0);
	  printf("          ROB fit intercept = %9.4f, offset = %f\n",q0,zeroshift);
	  //calculate root of quadratic fit
	  Float_t a,b,c;
	  pm2->Draw("same");

	  c1->Update();
	  
	  if(dotest) {
	    c3->cd(1);	    
	    MainTree->Draw(Form("%s >> h3",lname.Data()),Form("DetID==%i",DetNum),"");
	    h3->SetTitle(Form("Det%i",DetNum));  
	    h3->GetXaxis()->SetRangeUser(100,16000);
	    nfound = s->Search(h3,ssigma," ",sthresh);//run 1264
	    c3->Update();      
	    //sort peaks in order of channle number
	    xpeaks = s->GetPositionX();
	    Temp=0;
	    for(Int_t i=0;i<nfound;i++) {
	      for(Int_t j=i;j<nfound;j++) {
		if (xpeaks[j] < xpeaks[i]) {
		  Temp = xpeaks[i];
		  xpeaks[i] = xpeaks[j];
		  xpeaks[j] = Temp;  
		}
	      }
	    }
	    if(nfound==npeaks) 
	      FitGraph3 = new TGraph(nfound,xpeaks, &(Volts[0]));
	    else if(nfound==(npeaks-1)) 
	      FitGraph3 = new TGraph(nfound,xpeaks, &(Volts1[0]));
	    else if(nfound==(npeaks-2)) 
	      FitGraph3 = new TGraph(nfound,xpeaks, &(Volts2[0]));
	    else if(nfound==(npeaks-3)) 
	      FitGraph3 = new TGraph(nfound,xpeaks, &(Volts3[0]));
	    else if(nfound==(npeaks-4)) 
	      FitGraph3 = new TGraph(nfound,xpeaks, &(Volts4[0]));
	    else {
	      printf("number of peaks found = %d of %d\n", nfound, npeaks);
	      continue;
	    }
	    c3->cd(2);

	    FitGraph3->GetHistogram()->GetXaxis()->SetTitle("Peak positions (channel)");
	    FitGraph3->GetHistogram()->GetYaxis()->SetTitle("Pulser setting (voltage)");
	    FitGraph3->Draw("AP*");
	    FitGraph3->Fit("fit","qROB=0.95");
	    zeroshift = fit->GetParameter(0);
	    vperch = fit->GetParameter(1);
	    q0 = -zeroshift/vperch;
	    printf("         test fit intercept = %9.4f, offset = %f\n",q0,zeroshift);
	    c3->Update();      
	  }
	  if(dowait) {
	    c2->Update();      
	    c2->WaitPrimitive();
	  }
	  
	  cout << zeroshift << " " << vperch << endl;;
	  outfile << DetNum << "\t" << zeroshift << "\t" << vperch <<endl;
	}//end channel loop

  //delete FitGraph;
  //if(docentroid) delete FitGraph2;
  //delete s;
}//end SiPulser_All
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
