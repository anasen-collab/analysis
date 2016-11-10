/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
{
  //run250//
  //const Int_t npeaks = 6;
  //Float_t Volts[npeaks] = {0.5, 1.0, 3.0, 5.0, 7.0, 9.0};
  //TFile *file1 = new TFile("/data0/nabin/ANASEN_KTM/New/evt2root/run250_NSCL11_Pulser.root");  
  //TFile *file1 = new TFile("/data0/nabin/ANASEN/ANASEN_NKJ/New/evt2root/run250_NSCL11_Pulser.root");  

  //run262//
  //const Int_t npeaks = 5;
  //Float_t Volts[npeaks] = {0.05, 0.1, 0.15, 0.2, 0.25};
  //const Int_t npeaks = 4;
  //Float_t Volts[npeaks] = {0.1, 0.15, 0.2, 0.25};
  //TFile *file1 = new TFile("/data0/nabin/ANASEN/ANASEN_NKJ/ANASEN_N/262.root");  

  //run354//
  //const Int_t npeaks = 5;
  //Float_t Volts[npeaks] = {0.1, 0.2, 0.3, 0.4, 0.5};
  //Float_t Volts2[4] = {0.2, 0.3, 0.4, 0.5};
  //TFile *f1 = new TFile("/home2/parker/ANASEN/LSU/evt2root_root/run354.root");
  
  //run640//
  //const Int_t npeaks = 5;
  //Float_t Volts[npeaks] = {0.06, 0.12, 0.3, 0.6, 0.8};
  //Float_t Volts2[4] = {0.12, 0.3, 0.6, 0.8};
  //TFile *file1 = new TFile("/home/manasta/Desktop/parker_codes/evt2root_files/run640.root");  
  
  //run942 18Ne data
  //const Int_t npeaks = 6;
  //Float_t Volts[npeaks] = {0.006, 0.012, 0.03, 0.06, 0.08,0.12};
  //Float_t Volts2[4] = {0.03, 0.06, 0.08, 0.12};
  //TFile *file1 = new TFile("/home/manasta/Desktop/parker_codes/evt2root_files/run942.root"); 
  //TFile *file1 = new TFile("/data0/manasta/evt2root_files/run942.root"); 

  //run1036 24Mg data
  const Int_t npeaks = 6;
  Float_t Volts[npeaks] = {0.003, 0.006, 0.01, 0.03, 0.06, 0.1};
  Float_t Volts2[4] = {0.01, 0.03, 0.06, 0.1};
  TFile *file1 = new TFile("/data1/lighthall/root/run1036.root"); 
  
  TCanvas *c1 = new TCanvas();
  c1->Divide(1,2);

  TCanvas *c2 = new TCanvas("c2","click me",270,100);
  
  //TH1I *h1 = new TH1I("h1","h1",256,0,4095);
  //TH1I *h1 = new TH1I("h1","h1",512,0,6000);
  Int_t  size = 4096;
  TH1I *h1 = new TH1I("h1","h1",size,0,size);
  TTree *DataTree = (TTree*)file1->Get("DataTree");
    
  ofstream outfile;
  outfile.open("PCpulserCal.dat");

  TPolyMarker *pm;
  
  for (Int_t id=2; id<4; id++) {
    for (Int_t chan=0; chan<32; chan++) {
    
      if(id==3 && chan>15) continue;

      c1->cd(1);
      h1->GetXaxis()->UnZoom();
      //DataTree->Draw("ADC.Data>>h1",Form("ADC.ID==%d && ADC.ChNum==%d && ADC.Data>500 && ADC.Data<3900",id,chan));
      //DataTree->Draw("ADC.Data>>h1",Form("ADC.ID==%d && ADC.ChNum==%d && ADC.Data>180 && ADC.Data<4900",id,chan));
      DataTree->Draw("ADC.Data>>h1",Form("ADC.ID==%d && ADC.ChNum==%d && ADC.Data>200 && ADC.Data<5000",id,chan));
      if(h1->GetEntries()==0) continue;
      
      TSpectrum *s = new TSpectrum(npeaks+1);
      h1->SetTitle(Form("ADC %i Chan %i",id,chan));
      
      //Int_t nfound = s->Search(h1,9," ",0.09);
      //Int_t nfound = s->Search(h1,3," ",0.05); //Search(histo,sigma,option,threshold)
      Int_t nfound = s->Search(h1,3," ",0.05);
      c1->Update();
      //c1->WaitPrimitive();
      
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
      
      Float_t gpar[3*npeaks];
      for (Int_t i=0; i<npeaks; i++) {//fit each peak with a non-overlapping gaussian
	gfitc("h1",xpeaks[i],mnsp/2,"+q");
	for (Int_t j=0; j<3; j++) {
	  gpar[(3*i)+j]=gaus->GetParameter(j);
	}
      }
      
      const int npar=npeaks*3;
      Double_t par[npar];
      Float_t gpeaks[npeaks];//=xpeaks; 
      Bool_t docon=kTRUE;
      Float_t gfitwide=mnsp*3;
      Float_t gfitmin=xpeaks[0]-gfitwide;
      Float_t gfitmax=xpeaks[npeaks-1]+gfitwide;

      TString fform; //define functional form for global fit
      fform="gaus(0)";
      for (Int_t i=1; i<npeaks; i++) {
	fform+=Form("+gaus(%d)",3*i);
      }

      if((TF1 *)gROOT->FindObject("total"))total->Delete();
      TF1 *total = new TF1("total",fform,gfitmin,gfitmax);
      h1->GetXaxis()->SetRangeUser(gfitmin,gfitmax);

      total->SetLineColor(3);
      for (Int_t i=0; i<(npar); i++) {
	par[i]=gpar[i];
      }
      total->SetParameters(par);
      h1->Fit(total,"Mq+");

      //plot individual fits
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
	TGraph *graph = new TGraph(npeaks,xpeaks,Volts);
	graph->Draw("A*");
	TGraph *graph2 = new TGraph(npeaks,gpeaks,Volts);
	graph2->SetMarkerColor(4);
	graph2->Draw("*same");
      }else{
	TGraph *graph = new TGraph(4,xpeaks,Volts2);
	graph->Draw("A*");
      }
      //      TF1 *pol1 = new TF1("pol1","[0]+[1]*x",0,5000);
      graph->Fit("pol1","q");
      graph2->Fit("pol1","q+");
      //pol1->SetLineColor(4);
      c1->Update();
      //c1->WaitPrimitive();
      c2->Update();
      //c2->WaitPrimitive();
      
      outfile << id << "\t" << chan << "\t" << pol1->GetParameter(0) << "\t" << pol1->GetParameter(1) << endl;

      //TGraph *FitGraph = new TGraph(24,xpeaks,&(Volts[0]));
    }
  }
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
