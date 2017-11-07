////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ANASEN proportional counter (PC) pulser calibration.
// See readme.md for general instructions.
//
// Developed by : Jon Lighthall Nov 2016
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
{
  // 24Mg data July 11, 2016
  const Int_t npeaks = 7;
  Float_t Volts[npeaks] = {0,0.003, 0.006, 0.01, 0.03, 0.06, 0.1};
  TFile *file1 = new TFile("/data0/lighthall/root/raw/run1036.root");
  
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
  c1->SetWindowPosition(0,63);
  c1->SetWindowSize(1362,656);
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
  
  Int_t  size = 4096;
  TH1I *hist1 = new TH1I("hist1","hist1",size,0,size);
  
  TTree *DataTree = (TTree*)file1->Get("DataTree");

  TF1 *fit = new TF1("fit","pol1",0,size);
  TF1 *fit2 = new TF1("fit2","pol2",0,size);
  TF1 *fit3 = new TF1("fit3","pol1",0,size);
  TF1 *fit4 = new TF1("fit4","pol1",0,size);
  TF1 *fit5 = new TF1("fit5","pol1",0,size);
  fit2->SetLineColor(3);
  fit2->SetLineStyle(2);
  fit3->SetLineColor(4);
  fit3->SetLineStyle(2);
  fit4->SetLineColor(7);
  fit4->SetLineStyle(2);
  fit5->SetLineStyle(2);
  
  TGraph *FitGraph = 0; //for TSPectrum peaks
  TGraph *FitGraph2 = 0;//for gaussian fit peaks
  TPolyMarker *pm;
  
  ofstream outfile; //offsets, ROB fit
  ofstream outfile3;//offsets, full fit, all points equally weighted
  ofstream outfile4;//offsets, centroid ROB
  ofstream outfile2;//offsets, pedistal 
  ofstream outfile5;//offsets for zero peak
  outfile.open("saves/PCpulserCal.dat");
  outfile3.open("saves/PCpulserCal_full.dat");
  outfile4.open("saves/PCpulserCal_centroid.dat");
  outfile2.open("saves/PCpulserCal_zero.dat");
  outfile5.open("saves/PCpulserCal_zero_offset.dat");
  outfile << "ID\tChan\tVolt offset\tVolts/Chan" << endl;
  outfile3 << "ID\tChan\tVolt offset\tVolts/Chan" << endl;
  outfile4 << "ID\tChan\tVolt offset\tVolts/Chan" << endl;
  outfile2 << "ID\tChan\tVolt offset\tVolts/Chan\tPedi\tVolt" << endl;
  outfile5 << "ID\tChan\tVolt zero\tVolts/Chan" << endl;
  
  for (Int_t id=2; id<4; id++) {
    for (Int_t chan=0; chan<32; chan++) {
      //if(chan>2) continue;
      if(id==3 && chan>15) continue;
      cout << "ADC " << id << " Chan " << chan << " ";
      c1->cd(1);
      hist1->GetXaxis()->UnZoom();
      DataTree->Draw("ADC.Data>>hist1",Form("ADC.ID==%d && ADC.ChNum==%d && ADC.Data>0",id,chan));
      if(hist1->GetEntries()==0) continue;

      // Step 1: Search for peaks
      TSpectrum *s = new TSpectrum();
      hist1->SetTitle(Form("ADC %i Chan %i",id,chan));
      Int_t nfound = s->Search(hist1,3," ",0.05); //Search(histo,sigma,option,threshold)
      c1->Update();
      Float_t *xpeaks = s->GetPositionX();

      // order peaks
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
        if(nfound>npeaks) {
	    cout << "Warning: Wire " << id << "(" << chan << ")" << "has too many peaks" << endl;
	    break;
	  }
        //outfile << id << "\t" << chan << "\t" << Volts[l] << "\t" << xpeaks[l] << endl;
      }

      // Step 2: Fit each peak with a non-overlapping gaussian
      Float_t mnsp=25; //minimum space between adjacent peaks
      mnsp=xpeaks[npeaks-1]-xpeaks[0];
      for (Int_t i=0; i<(npeaks-1); i++){//find min. peak spacing
	if(((xpeaks[i+1]-xpeaks[i])<mnsp)&&((i+1)<npeaks))
	  mnsp=xpeaks[i+1]-xpeaks[i];
      }

      Float_t gpar[npar];
      for (Int_t i=0; i<npeaks; i++) {//fit each peak with a non-overlapping gaussian
	hist1->Fit("gaus","+q","",xpeaks[i]-mnsp/2,xpeaks[i]+mnsp/2);
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

      // Step 3: Calculate global fit
      if((TF1 *)gROOT->FindObject("total"))total->Delete();
      TF1 *total = new TF1("total",fform,gfitmin,gfitmax);
      hist1->GetXaxis()->SetRangeUser(gfitmin,gfitmax);

      total->SetLineColor(3);
      for (Int_t i=0; i<(npar); i++) {
	par[i]=gpar[i];
      }
      total->SetParameters(par);
      hist1->Fit(total,"Mq+");

      // plot individual fits using parameters from global fit
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
      
      // Step 4: Calculate calibration constants
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
	if(nfound==1)
	  printf(" peak zero at %d\n", xpeaks[0]);
	continue;
      }
      
      FitGraph->GetHistogram()->GetXaxis()->SetTitle("Channel positions from peak fit");
      FitGraph->GetHistogram()->GetYaxis()->SetTitle("Voltages from calibration file");
      
      FitGraph->SetTitle(Form("Fit ADC %i Chan %i",id,chan));
      FitGraph->GetYaxis()->SetRangeUser(-0.01,0.11);
      FitGraph->Draw("A*");
      FitGraph->Fit("fit","qROB");
      FitGraph->Fit("fit2","q");
      FitGraph->Fit("fit3","q");
      
      //centroid fit
      FitGraph2->SetMarkerColor(4);
      FitGraph2->Draw("*same");
      FitGraph2->Fit("fit4","qROB");//ROB may give errors
      fit->Draw("same");
      fit2->Draw("same");
      fit3->Draw("same");

      Float_t slope=fit->GetParameter(1);
      Float_t offset=fit->GetParameter(0);
      printf(" Fit parameters are: Slope = %f volts per channel, Offset = %f volts\n",slope,offset);
      printf(" Inverse fit parameters are slope %f, offset %f\n",1/slope,-offset/slope); 

      //offset=xpeaks[0]*slope;
      fit5->FixParameter(0,-xpeaks[0]*slope);
      fit5->FixParameter(1,slope);
      fit5->Draw("same");
      
      leg = new TLegend(0.1,0.6,0.2,0.9);
      leg->AddEntry(fit,"pol1, ROB","l");
      leg->AddEntry(fit2,"pol2","l");
      leg->AddEntry(fit3,"pol1","l");   
      leg->AddEntry(fit4,"pol1 centroid, ROB","l");
      leg->AddEntry(fit5,"pol1 zero","l");
      leg->AddEntry(FitGraph,"peaks","p");
      leg->AddEntry(FitGraph2,"centroids","p");
      leg->Draw();
      
      // Step 5: test fit
      printf(" Testing fit:\n");
      for (Int_t i=0; i<npeaks; i++) {
       	printf("  Peak %2d at %6.1f is %9.6f (%9.6f)\n",i,xpeaks[i],xpeaks[i]*slope+offset,(xpeaks[i]*slope+offset)-Volts[i]);
      }
            
      c1->Update();
      if(dowait) {
	c2->Update();      
	c2->WaitPrimitive();
      }
      
      outfile  << id << "\t" << chan << "\t" << fit->GetParameter(0)  << "\t" << fit->GetParameter(1)  << endl;
      outfile3 << id << "\t" << chan << "\t" << fit3->GetParameter(0) << "\t" << fit3->GetParameter(1) << endl;
      outfile4 << id << "\t" << chan << "\t" << fit4->GetParameter(0) << "\t" << fit4->GetParameter(1) << endl;
      outfile2 << id << "\t" << chan << "\t" << fit->GetParameter(0) << "\t" << fit->GetParameter(1) << "\t"<< xpeaks[0] << "\t" << xpeaks[0]*slope+fit->GetParameter(0) << endl;
      outfile5 << id << "\t" << chan << "\t" << -xpeaks[0]*slope << "\t" << fit->GetParameter(1) << endl;
    }
  }
  if(dowait)
    c2->Close();
}
