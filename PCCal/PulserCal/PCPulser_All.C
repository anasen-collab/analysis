/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//// ANASEN proportional counter (PC) pulser calibration.
////
//// Output file (e.g."PCpulserCal.dat") has the following columns:
//// ID, Chan, zero shift (fit offset), voltage per can (fit slope)
////
//// Usage: root -l PCPulser_All.C 
////
//// Edited by : Jon Lighthall Nov 2016
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
{
  //run250// December 4, 2015
  //const Int_t npeaks = 6;
  //Float_t Volts[npeaks] = {0.5, 1.0, 3.0, 5.0, 7.0, 9.0};
  //TFile *file1 = new TFile("/data0/nabin/ANASEN_KTM/New/evt2root/run250_NSCL11_Pulser.root");  
  //TFile *file1 = new TFile("/data0/nabin/ANASEN/ANASEN_NKJ/New/evt2root/run250_NSCL11_Pulser.root");  

  //run262// January 19, 2016
  //const Int_t npeaks = 5;
  //Float_t Volts[npeaks] = {0.05, 0.1, 0.15, 0.2, 0.25};
  //const Int_t npeaks = 4;
  //Float_t Volts[npeaks] = {0.1, 0.15, 0.2, 0.25};
  //TFile *file1 = new TFile("/data0/nabin/ANASEN/ANASEN_NKJ/ANASEN_N/262.root");  

  //run354// February 5, 2016
  //const Int_t npeaks = 5;
  //Float_t Volts[npeaks] = {0.1, 0.2, 0.3, 0.4, 0.5};
  //Float_t Volts2[4] = {0.2, 0.3, 0.4, 0.5};
  //TFile *f1 = new TFile("/home2/parker/ANASEN/LSU/evt2root_root/run354.root");
  
  //run640// May 10, 2016
  //const Int_t npeaks = 5;
  //Float_t Volts[npeaks] = {0.06, 0.12, 0.3, 0.6, 0.8};
  //Float_t Volts2[4] = {0.12, 0.3, 0.6, 0.8};
  //TFile *file1 = new TFile("/home/manasta/Desktop/parker_codes/evt2root_files/run640.root");  
  
  //run942 18Ne data May 31, 2016
  //const Int_t npeaks = 6;
  //Float_t Volts[npeaks] = {0.006, 0.012, 0.03, 0.06, 0.08,0.12};
  //Float_t Volts2[4] = {0.03, 0.06, 0.08, 0.12};
  //TFile *file1 = new TFile("/home/manasta/Desktop/parker_codes/evt2root_files/run942.root"); 
  //TFile *file1 = new TFile("/data0/manasta/evt2root_files/run942.root"); 

  //run1036 24Mg data July 11, 2016
  const Int_t npeaks = 6;
  Float_t Volts[npeaks] = {0.003, 0.006, 0.01, 0.03, 0.06, 0.1};
  TFile *file1 = new TFile("/data1/lighthall/root/run1036.root"); 
  
  const int npar=npeaks*3;

  Float_t Volts1[npeaks-1];//(n-1) peaks, replaces Volts5
  for (Int_t i=0;i<(npeaks-1);i++) {
    Volts1[npeaks-2-i]=Volts[npeaks-i-1];
    //printf("Volts1[%d]=%f\n",npeaks-2-i,Volts1[npeaks-2-i]);
  }
  
  Float_t Volts2[npeaks-2];//(n-2 peaks) replaces Volts4
  for (Int_t i=0;i<(npeaks-2);i++) {
    Volts2[npeaks-3-i]=Volts[npeaks-i-1];
    //printf("Volts2[%d]=%f\n",npeaks-3-i,Volts2[npeaks-3-i]);
  }

  TCanvas *c1 = new TCanvas();
  c1->Divide(1,2);

  Bool_t dowait=kFALSE; //wait betweeen fits
  
  if(dowait) {
    TCanvas *c2 = new TCanvas("c2","double-click me",260,100);
    c2->cd();
    //TButton *but1 = new TButton("next fit","",.05,.55,.45,.9);
    //but1->Draw();
    TButton *but2 = new TButton("quit ROOT",".q",.05,.1,.45,.45);
    but2->Draw();
  }
  
  //TH1I *h1 = new TH1I("h1","h1",256,0,4095);
  //TH1I *h1 = new TH1I("h1","h1",512,0,6000);
  Int_t  size = 4096;
  TH1I *h1 = new TH1I("h1","h1",size,0,size);
  
  TTree *DataTree = (TTree*)file1->Get("DataTree");

  TF1 *fit = new TF1("fit","pol1",0,size/2);
  TF1 *fit2 = new TF1("fit2","pol2",0,size/2);
  TF1 *fit3 = new TF1("fit3","pol1",0,size/2);
  TF1 *fit4 = new TF1("fit4","pol1",0,size/2);
  fit2->SetLineColor(3);
  fit2->SetLineStyle(2);
  fit3->SetLineColor(4);
  fit3->SetLineStyle(2);
  fit4->SetLineColor(7);
  fit4->SetLineStyle(2);
  
  TGraph *FitGraph = 0;
  TGraph *FitGraph2 = 0;
  
  ofstream outfile;
  ofstream outfile3;//offsets, full fit
  ofstream outfile4;//offsets, centroid ROB
  outfile.open("PCpulserCal.dat");
  outfile3.open("PCpulserCal_full.dat");
  outfile4.open("PCpulser_centroid.dat");
  
  TPolyMarker *pm;
  
  for (Int_t id=2; id<4; id++) {
    for (Int_t chan=0; chan<32; chan++) {
    
      if(id==3 && chan>15) continue;

      c1->cd(1);
      h1->GetXaxis()->UnZoom();
      //DataTree->Draw("ADC.Data>>h1",Form("ADC.ID==%d && ADC.ChNum==%d && ADC.Data>500 && ADC.Data<3900",id,chan));
      //DataTree->Draw("ADC.Data>>h1",Form("ADC.ID==%d && ADC.ChNum==%d && ADC.Data>180 && ADC.Data<4900",id,chan));
      //DataTree->Draw("ADC.Data>>h1",Form("ADC.ID==%d && ADC.ChNum==%d && ADC.Data>200 && ADC.Data<5000",id,chan));
      DataTree->Draw("ADC.Data>>h1",Form("ADC.ID==%d && ADC.ChNum==%d && ADC.Data>200",id,chan));
      if(h1->GetEntries()==0) continue;
      
      TSpectrum *s = new TSpectrum();
      h1->SetTitle(Form("ADC %i Chan %i",id,chan));
      
      Int_t nfound = s->Search(h1,3," ",0.05); //Search(histo,sigma,option,threshold)
      c1->Update();
      
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
	    
      for (Int_t l=0; l<nfound; l++) {
        if(nfound>npeaks)
	  {
	    cout << "Warning: Wire " << id << "(" << chan << ")" << "has too many peaks" << endl;
	    break;
	  }
        //outfile << id << "\t" << chan << "\t" << Volts[l] << "\t" << xpeaks[l] << endl;
      }

      Float_t mnsp=25; //minimum space between adjacent peaks
      mnsp=xpeaks[npeaks-1]-xpeaks[0];
      for (Int_t i=0; i<(npeaks-1); i++){//find min. peak spacing
	if(((xpeaks[i+1]-xpeaks[i])<mnsp)&&((i+1)<npeaks))
	  mnsp=xpeaks[i+1]-xpeaks[i];
      }

      Float_t gpar[npar];
      for (Int_t i=0; i<npeaks; i++) {//fit each peak with a non-overlapping gaussian
	h1->Fit("gaus","+q","",xpeaks[i]-mnsp/2,xpeaks[i]+mnsp/2);
	for (Int_t j=0; j<3; j++) {
	  gpar[(3*i)+j]=gaus->GetParameter(j);
	}
      }
            
      Double_t par[npar];
      Float_t gpeaks[npeaks];//=xpeaks; 
      Float_t gfitwide=mnsp*3;
      Float_t gfitmin=xpeaks[0]-gfitwide;
      Float_t gfitmax=xpeaks[npeaks-1]+gfitwide;

      TString fform; //define functional form for global fit
      fform="gaus(0)";
      for (Int_t i=1; i<npeaks; i++) {
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
      h1->Fit(total,"Mq+");

      //plot individual fits using parameters from global fit
      TF1 **functions = new TF1*[npeaks];
      pm=new TPolyMarker(npeaks);
      pm->SetMarkerStyle(23);
      pm->SetMarkerSize(1.2);
      pm->SetMarkerColor(4);
      for (Int_t i=0;i<npeaks;i++) {
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
      
      c1->cd(2);
      if (nfound==npeaks) {
	FitGraph = new TGraph(npeaks,xpeaks,Volts);
	FitGraph2 = new TGraph(npeaks,gpeaks,Volts);
      }
      else if(nfound==(npeaks-1)) {
	FitGraph = new TGraph(nfound,xpeaks, &(Volts1[0]));
	FitGraph2 = new TGraph(nfound,gpeaks, &(Volts1[0]));
      }
      else if(nfound==(npeaks-2)) {
	FitGraph = new TGraph(nfound,xpeaks, &(Volts2[0]));
	FitGraph2 = new TGraph(nfound,gpeaks, &(Volts2[0]));
      }
      else {
	printf("number of peaks found = %d of %d\n", nfound, npeaks);
	continue;
      }

      //TF1 *pol1 = new TF1("pol1","[0]+[1]*x",0,5000);
      //FitGraph->Fit("pol1","q");
      //FitGraph2->Fit("pol1","Mq+");
      FitGraph->SetTitle(Form("Fit ADC %i Chan %i",id,chan));
      FitGraph->Draw("A*");
      FitGraph->Fit("fit","qROB");
      FitGraph->Fit("fit2","q");
      FitGraph->Fit("fit3","q");
      
      FitGraph2->SetMarkerColor(4);
      FitGraph2->Draw("*same");
      FitGraph2->Fit("fit4","qROB");
      
      leg = new TLegend(0.1,0.75,0.2,0.9);
      leg->AddEntry(fit,"pol1, ROB","l");
      leg->AddEntry(fit2,"pol2","l");
      leg->AddEntry(fit3,"pol1","l");   
      leg->AddEntry(fit4,"pol1 centroid, ROB","l");
      leg->AddEntry(FitGraph,"peaks","p");
      leg->AddEntry(FitGraph2,"centroids","p");
      leg->Draw();
            
      c1->Update();
      if(dowait) {
	c2->Update();      
	c2->WaitPrimitive();
      }
      
      outfile  << id << "\t" << chan << "\t" << fit->GetParameter(0)  << "\t" << fit->GetParameter(1)  << endl;
      outfile3 << id << "\t" << chan << "\t" << fit3->GetParameter(0) << "\t" << fit3->GetParameter(1) << endl;
      outfile4 << id << "\t" << chan << "\t" << fit4->GetParameter(0) << "\t" << fit4->GetParameter(1) << endl;

      //TGraph *FitGraph = new TGraph(24,xpeaks,&(Volts[0]));
    }
  }
   if(dowait)
     c2->Close();
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
