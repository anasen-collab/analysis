//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//// ASICS all channels Pulser calibration.
////
//// Output file (e.g."Sipulser_2015Dec13.dat") has the following columns:
//// MBID, CBID, ASICs_Channel, ZeroShift(offset), Voltage_per_Ch(slope)
////
//// Usage: root -l SiPulser_All.C (from the same directory).
////
//// Edited by : Nabin Rijal , 2015Dec13
////             Jon Lighthall Nov 2016
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#include <TMath.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1I.h>
#include <TGraph.h>
#include <TF1.h>
#include <TSpectrum.h>
#include <fstream>

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void SiPulser_All (void)
{
  //MBID CID ASICs_Chn ZeroShift VperCh

  //run250//run251 //
  //const Int_t npeaks = 6;
  //Float_t Volts[npeaks] = { 0.5, 1.0, 3.0, 5.0, 7.0, 9.0 };
  //Float_t Volts5[5] = { 1.0, 3.0, 5.0, 7.0, 9.0 };
  //Float_t Volts4[4] = { 3.0, 5.0, 7.0, 9.0 };
  //TFile *f1 = new TFile("/data0/nabin/ANASEN/ANASEN_NKJ/New/evt2root/run251_NSCL11_Pulser.root");//back

  //run 585
  //const Int_t npeaks = 6;
  //Float_t Volts[npeaks] = { 0.5, 1.0, 3.0, 5.0, 7.0, 9.0 };
  //Float_t Volts5[5] = { 1.0, 3.0, 5.0, 7.0, 9.0 };
  //Float_t Volts4[4] = { 3.0, 5.0, 7.0, 9.0 };
  //TFile *f1 = new TFile("run585.root");//front

  //run643 - negative pol//run645 - positive pol//
  // const Int_t npeaks = 8;
  // Float_t Volts[npeaks] = { 0.2, 0.3, 0.6, 1.2, 1.5, 2.0, 3.0, 3.5 };
  // Float_t Volts5[5] = { 1.2, 1.5, 2.0, 3.0, 3.5 };
  // Float_t Volts4[4] = { 1.5, 2.0, 3.0, 3.5 };
  // TFile *f1 = new TFile("/home/manasta/Desktop/parker_codes/evt2root_files/run643.root");//front

  //run 1034-5
  const Int_t npeaks = 8;
  Float_t Volts[npeaks] = { 0.5, 0.8, 1.0, 1.5, 3.0, 5.5, 7.0, 9.0};
  TFile *f1 = new TFile("/data1/lighthall/root/run1034.root");
  //TFile *f2 = new TFile("/data1/lighthall/root/run1034.root"); 
  
  //run1262-4
  //const Int_t npeaks = 5;
  //Float_t Volts[npeaks] = { 0.5, 1.5, 3.0, 7.0, 9.0};
  //TFile *f1 = new TFile("/data1/lighthall/root/run1262.root");
  //TFile *f1 = new TFile("/data1/lighthall/root/run1263.root");

  const int npar=npeaks*3;
  
  Float_t Volts1[npeaks-1];//(n-1) peaks, replaces Volts5
  for (Int_t i=0;i<(npeaks-1);i++) {
    Volts1[npeaks-1-i-1]=Volts[npeaks-i-1];
    //printf("Volts1[%d]=%f\n",npeaks-1-i-1,Volts1[npeaks-1-i-1]);
  }

  Float_t Volts2[npeaks-2];//(n-2 peaks) replaces Volts4
  for (Int_t i=0;i<(npeaks-2);i++) {
    Volts2[npeaks-2-i-1]=Volts[npeaks-i-1];
    //printf("Volts2[%d]=%f\n",npeaks-2-i-1,Volts2[npeaks-2-i-1]);
  }

  Float_t Volts3[npeaks-3];//(n-3 peaks) 
  for (Int_t i=0;i<(npeaks-3);i++) {
    Volts3[npeaks-3-i-1]=Volts[npeaks-i-1];
    //printf("Volts3[%d]=%f\n",npeaks-3-i-1,Volts3[npeaks-3-i-1]);
  }

  Float_t Volts4[npeaks-4];//(n-4 peaks) 
  for (Int_t i=0;i<(npeaks-4);i++) {
    Volts4[npeaks-4-i-1]=Volts[npeaks-i-1];
    //printf("Volts4[%d]=%f\n",npeaks-4-i-1,Volts4[npeaks-4-i-1]);
  }
  
  TCanvas *c1 = new TCanvas();
  c1->SetWindowPosition(0,63);
  c1->SetWindowSize(1362,656);
  c1->Divide(1,2);

  Bool_t dowait=kFALSE; //wait betweeen fits
  Bool_t docentroid=kTRUE; //use centroid fit
  Bool_t docomp=kFALSE; //compare two files
  
  if(dowait) {
    TCanvas *c2 = new TCanvas("c2","double-click me",270,100);
    c2->cd();
    TButton *but2 = new TButton("quit ROOT",".q",.05,.1,.45,.45);
    but2->Draw();
  }

  Int_t  size = 4096*4;
  TH1I *h1 = new TH1I("h1","h1",size,0,size);//min was 300
  TTree *DataTree = (TTree*)f1->Get("DataTree");
  if(docomp) {
    TH1I *h2 = new TH1I("h2","h2",size,0,size);
    h2->SetLineColor(kRed+3);
    TTree *DataTree2 = (TTree*)f2->Get("DataTree");
  }
  
  TF1 *fit = new TF1("fit","pol1",0,size/2);
  TF1 *fit2 = new TF1("fit2","pol2",0,size/2);
  TF1 *fit3 = new TF1("fit3","pol1",0,size/2);
  fit2->SetLineColor(3);
  fit2->SetLineStyle(2);
  fit3->SetLineColor(4);
  fit3->SetLineStyle(2);
  
  Double_t zeroshift = 0;//zero voltage shift (fit offset)
  Double_t vperch = 0;//volts per channel (fit slope)
  Double_t q0 = 0;
  TSpectrum *s = 0;
  TGraph *FitGraph = 0;
  
  ofstream outfile; //offsets, ROB
  ofstream outfile2;//peaks
  ofstream outfile3;//offsets, full fit
  ofstream outfile4;//offsets, centroid ROB
  ofstream outfile5;//peaks, centroid
  
  outfile.open("Sipulser_offsets.dat");
  outfile2.open("Sipulser_peaklocations.dat");
  outfile3.open("Sipulser_offsets_full.dat");
  
  if(docentroid) {
    TF1 *fit4 = new TF1("fit4","pol1",0,size/2);
    fit4->SetLineColor(7);
    fit4->SetLineStyle(2);
    TGraph *FitGraph2 = 0;
    outfile4.open("Sipulser_offsets_centroid.dat");
    outfile5.open("Sipulser_peaklocations_centroid.dat");
    TPolyMarker *pm; //centroid peak locations  
  }
    
  TPolyMarker *pm2; //calculated intercepts
  
  for (Int_t MBID=1; MBID<3; MBID++) {  
    for (Int_t CBID=1; CBID<15; CBID++) { 
      ////for the front of the detectors, positive polarity 
      //if((CBID==1 || CBID==2 || CBID==5 || CBID==6 || CBID==9 ||CBID==10) ||(MBID==1 && (CBID==11 || CBID==12 || CBID==13 || CBID==14))){//run250 //run645 // run 1035 //run 1264
      
      ////for the back of the detectors, negative polarity 
      if((CBID==3 || CBID==4 || CBID==7 || CBID==8) || (MBID==2 && (CBID==11 || CBID==12))) {//run251 //run643 //1034 //run 1262
	
	for (Int_t ChNum=0; ChNum<16; ChNum++) {
	  //// Mask bad channels  //that can create problem in cruising calibration.
	  //if(!(MBID==2 && CBID==10 && ChNum==12))continue;
	  //if(MBID==1 && CBID<9)continue;
	  //if(MBID==1 && CBID==9 &&ChNum<9)continue;
	  //if(MBID==1 || CBID>8)continue;
	  
	  //if(MBID==1 && CBID==1 && (ChNum==4 ||ChNum ==14)) continue;//bad at front run250
	  //if((MBID==1 && CBID==8 && (ChNum==0 || ChNum==1 || ChNum==2))||(MBID==2 && CBID ==8 && (ChNum==1 || ChNum==2)))continue;//bad at back run251
	  //if(MBID == 2 && CBID == 7 && (ChNum ==10))continue;//bad at back run251
	      
	  //if(MBID==1 && CBID==1 && ChNum ==14) continue;//run 643	      
	  //if(MBID==1 && CBID==8 && (ChNum==0 || ChNum==15))continue; //run 643
	  //if(MBID ==1 && CBID ==12  && ChNum ==6) continue;//run 643
	  //if(MBID ==2 && CBID ==7  && (ChNum ==4 || ChNum ==5 || ChNum ==9 || ChNum ==10))continue;//run 643
	  //if(MBID ==2 && CBID ==12  && (ChNum ==13 || ChNum ==14))continue;//run 643

	  //runs 1034 1262
	  //if(MBID==1 && CBID==8 && ChNum==2)continue;
	  //if(MBID==2 && CBID==8 && ChNum==2)continue;

	  //runs 1264
	  if(MBID==1 && CBID==1 && ChNum==4)continue;
	  //if(MBID==1 && CBID==8 && ChNum==2)continue;
	  //if(MBID==2 && CBID==8 && ChNum==2)continue;
	      
	  c1->cd(1);	    
	  DataTree->Draw("Si.Energy>>h1",Form("Si.Energy>1 && Si.MBID==%d && Si.CBID==%d && Si.ChNum==%d",MBID,CBID,ChNum));
	  if(docomp)
	    DataTree2->Draw("Si.Energy>>h2",Form("Si.MBID==%d && Si.CBID==%d && Si.ChNum==%d",MBID,CBID,ChNum),"same");
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

	  TSpectrum *s = new TSpectrum();
	      
	  h1->SetTitle(Form("MBID %d CBID %d ChNum %d",MBID,CBID,ChNum));
	  if(docomp)
	    h2->SetTitle(Form("MBID %d CBID %d ChNum %d",MBID,CBID,ChNum));  

	  //Int_t nfound = s->Search(h1,15,"",0.25);
	  //Int_t nfound = s->Search(h1,5," nobackground",0.10);
	  h1->GetXaxis()->UnZoom();
	  //s->SetResolution(3);
	  //Int_t nfound = s->Search(h1,9," ",0.15);//run 643?
	  //Int_t nfound = s->Search(h1,10," ",0.15);//run 1034
	  Int_t nfound = s->Search(h1,10," ",0.05);//run 1264

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

	  outfile2 << MBID << "\t" << CBID << "\t" << ChNum;// << nfound;
	  for (int i=0; i<nfound; i++) {
	    outfile2 << "\t" << xpeaks[i];
	  }
	  outfile2 << endl;

	  if(docentroid) {
	    Int_t nstep=npeaks;
	    if(nfound<npeaks)
	      nstep=nfound;
	    
	    Float_t mnsp=25; //minimum space between adjacent peaks
	    mnsp=xpeaks[nstep-1]-xpeaks[0];
	    mnsp=500;//for run 1262
	    for (Int_t i=0; i<(nstep-1); i++) {//find min. peak spacing
	      if(((xpeaks[i+1]-xpeaks[i])<mnsp)&&((i+1)<npeaks))
		mnsp=xpeaks[i+1]-xpeaks[i];
	    }
	    
	    Float_t gpar[npar];
	    for (Int_t i=0; i<nstep; i++) {//fit each peak with a non-overlapping gaussian
	      h1->Fit("gaus","+q","",xpeaks[i]-mnsp/2,xpeaks[i]+mnsp/2);
	      for (Int_t j=0; j<3; j++) {
		gpar[(3*i)+j]=gaus->GetParameter(j);
	      }
	    }
	    
	    Double_t par[npar];
	    Float_t gpeaks[npeaks];//=xpeaks; 
	    Float_t gfitwide=mnsp*3;
	    Float_t gfitmin=xpeaks[0]-gfitwide;
	    Float_t gfitmax=xpeaks[nstep-1]+gfitwide;
	    
	    TString fform; //define functional form for global fit
	    fform="gaus(0)";
	    for (Int_t i=1; i<nstep; i++) {
	      fform+=Form("+gaus(%d)",3*i);
	    }
	    
	    //global fit
	    if((TF1 *)gROOT->FindObject("total"))total->Delete();
	    TF1 *total = new TF1("total",fform,gfitmin,gfitmax);
	    h1->GetXaxis()->SetRangeUser(gfitmin,gfitmax);
	    
	    total->SetLineColor(3);
	    for (Int_t i=0; i<(npar); i++) {
	      par[i]=gpar[i];
	    }
	    total->SetParameters(par);
	    h1->Fit(total,"Mq+","");
	    
	    //plot individual fits using parameters from global fit
	    TF1 **functions = new TF1*[nstep];
	    pm=new TPolyMarker(npeaks);
	    pm->SetMarkerStyle(23);
	    pm->SetMarkerSize(1.2);
	    pm->SetMarkerColor(4);
	    for (Int_t i=0;i<nstep;i++) {
	      char fnames[20];
	      sprintf(fnames,"f%d",i);
	      functions[i] = new TF1(fnames,"gaus",gfitmin,gfitmax);
	      for (Int_t j=0;j<3;j++) {
		par[j+3*i]=total->GetParameter(j+3*i);
		functions[i]->SetParameter(j,par[j+3*i]);
	      }
	      gpeaks[i]=par[1+3*i];
	      //printf("Peak %d is %f or %f (%f)\n",i,xpeaks[i],gpeaks[i],xpeaks[i]-gpeaks[i]);
	      pm->SetPoint(i,par[1+3*i],par[0+3*i]);
	      
	      functions[i]->SetLineColor(1);
	      functions[i]->SetLineStyle(2);
	      functions[i]->Draw("same");
	    }
	    pm->Draw("same");
	    
	    outfile5 << MBID << "\t" << CBID << "\t" << ChNum;
	    for (int i=0; i<nfound; i++) {
	      outfile5 << "\t" << gpeaks[i];
	    }
	    outfile5 << endl;
	  }
	  
	  if(FitGraph!=0) {
	    delete FitGraph;
	  }

	  if(nfound==npeaks) {
	    FitGraph = new TGraph(nfound,xpeaks, &(Volts[0]));
	    if(docentroid)
	      FitGraph2 = new TGraph(nfound,gpeaks, &(Volts[0]));
	  }
	  else if(nfound==(npeaks-1)) {
	    FitGraph = new TGraph(nfound,xpeaks, &(Volts1[0]));
	    if(docentroid)
	      FitGraph2 = new TGraph(nfound,gpeaks, &(Volts1[0]));
	  }
	  else if(nfound==(npeaks-2)) {
	    FitGraph = new TGraph(nfound,xpeaks, &(Volts2[0]));
	    if(docentroid)
	      FitGraph2 = new TGraph(nfound,gpeaks, &(Volts2[0]));
	  }
	  else if(nfound==(npeaks-3)) {
	    FitGraph = new TGraph(nfound,xpeaks, &(Volts3[0]));
	    if(docentroid) 
	      FitGraph2 = new TGraph(nfound,gpeaks, &(Volts3[0]));
	  }
	  else if(nfound==(npeaks-4)) {
	    FitGraph = new TGraph(nfound,xpeaks, &(Volts4[0]));
	    if(docentroid)
	      FitGraph2 = new TGraph(nfound,gpeaks, &(Volts4[0]));
	  }
	  else {
	    printf("number of peaks found = %d of %d\n", nfound, npeaks);
	    outfile  << MBID << "\t" << CBID << "\t" << ChNum << "\t"<< endl;
	    outfile3 << MBID << "\t" << CBID << "\t" << ChNum << "\t"<< endl;
	    if(docentroid) 
	      outfile4 << MBID << "\t" << CBID << "\t" << ChNum << "\t"<< endl;
	    continue;
	  }
		
	  pm2=new TPolyMarker(4);
	  pm2->SetMarkerStyle(2);
	  pm2->SetMarkerSize(1.2);
	  pm2->SetMarkerColor(6);
	      
	  //h2->Draw("same");
	  c1->cd(2);
      
	  FitGraph->SetTitle(Form("Fit MBID %d CBID %d ChNum %d",MBID,CBID,ChNum));
	  FitGraph->GetHistogram()->GetXaxis()->SetTitle("Peak positions (channel)");
	  FitGraph->GetHistogram()->GetYaxis()->SetTitle("Pulser setting (voltage)");
    	  FitGraph->Draw("AP*");
	  FitGraph->Fit("fit","qROB=0.95");
	  FitGraph->Fit("fit2","q+");
	  FitGraph->Fit("fit3","q+");

	  leg = new TLegend(0.1,0.75,0.2,0.9);
	  leg->AddEntry(fit,"pol1, ROB=95","l");
	  leg->AddEntry(fit2,"pol2","l");
	  leg->AddEntry(fit3,"pol1","l");   
	  if(docentroid) {
	    FitGraph2->Draw("* same");
	    FitGraph2->Fit("fit4","qROB=0.95");
	    leg->AddEntry(fit4,"pol1 centroid, ROB=95","l");
	  }
	  leg->Draw();
		
	  ////FitGraph->Fit("fit","E");
	  zeroshift = fit->GetParameter(0);
	  vperch = fit->GetParameter(1);
	  q0 = -zeroshift/vperch;
	  pm2->SetPoint(0,q0,0);
	  printf("          ROB fit intercept = %9.4f, offset = %f\n",q0,zeroshift);
	  //calculate root of quadratic fit
	  Float_t a,b,c;
	  a=fit2->GetParameter(2);
	  b=fit2->GetParameter(1);
	  c=fit2->GetParameter(0);
	  Float_t r1=(-b+TMath::Sqrt(TMath::Power(b,2)-4*a*c))/(2*a);
	  pm2->SetPoint(1,r1,0);
	  printf("         quad fit intercept = %9.4f\n",r1);
	  Float_t zeroshiftf = fit3->GetParameter(0);
	  Float_t vperchf = fit3->GetParameter(1);
	  Float_t q0f=-zeroshiftf/vperchf;
	  pm2->SetPoint(2,q0f,0);
	  printf("  full linear fit intercept = %9.4f, offset = %f\n",q0f,zeroshiftf);
	  if(docentroid) {
	    Float_t zeroshiftc = fit4->GetParameter(0);
	    Float_t vperchc = fit4->GetParameter(1);
	    Float_t q0c=-zeroshiftc/vperchc;
	    pm2->SetPoint(3,q0c,0);
	    printf(" centroid ROB fit intercept = %9.4f, offset = %f\n",q0c,zeroshiftc);
	  }
	  pm2->Draw("same");

	  c1->Update();
	  if(dowait) {
	    c2->Update();      
	    c2->WaitPrimitive();
	  }
	  
	  cout << zeroshift << " " << vperch << endl;;
	  cout << MBID << " " << CBID << " " << ChNum << " q0 = "<<q0<<endl;
	      
	  outfile << MBID << "\t" << CBID << "\t" << ChNum << "\t"<< zeroshift << "\t" << vperch <<endl;
	  //outfile << MBID << "\t" << CBID << "\t" << ChNum << "\t"<< zeroshift << "\t" << vperch << "\t"<< q0 <<endl;
	  outfile3 << MBID << "\t" << CBID << "\t" << ChNum << "\t"<< zeroshiftf << "\t" << vperchf <<endl;
	  if(docentroid)
	    outfile4 << MBID << "\t" << CBID << "\t" << ChNum << "\t"<< zeroshiftc << "\t" << vperchc <<endl;
	}//end channel loop
      }//end polarity check
    }//end CBID loop
  }//end MBID loop
  
  //delete FitGraph;
  //if(docentroid) delete FitGraph2;
  //delete s;
}//end SiPulser_All
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
