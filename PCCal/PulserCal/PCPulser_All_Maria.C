/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
{
  //TFile *f1 = new TFile("/data0/nabin/ANASEN_KTM/New/evt2root/run250_NSCL11_Pulser.root");  
  //TFile *f1 = new TFile("/data0/nabin/ANASEN/ANASEN_NKJ/New/evt2root/run250_NSCL11_Pulser.root");

  //TFile *f1 = new TFile("/data0/nabin/ANASEN/ANASEN_NKJ/ANASEN_N/262.root");  
  
 //run942 18Ne data
  const Int_t npeaks = 6;
  Float_t Volts[npeaks] = {0.006, 0.012, 0.03, 0.06, 0.08,0.12};
  Float_t Volts2[4] = {0.03, 0.06, 0.08, 0.12};
  //TFile *f1 = new TFile("/home/manasta/Desktop/parker_codes/evt2root_files/run942.root"); 
  TFile *f1 = new TFile("/data0/manasta/evt2root_files/run942.root");

  //run1036 24Mg data
  //const Int_t npeaks = 6;
  //Float_t Volts[npeaks] = {0.003, 0.006, 0.01, 0.03, 0.06, 0.1};
  //Float_t Volts2[4] = {0.01, 0.03, 0.06, 0.1};
    
  TCanvas *c1 = new TCanvas();
  c1->Divide(1,2);

  TH1I *h1 = new TH1I("h1","h1",512,0,6000);
  TTree *DataTree = (TTree*)f1->Get("DataTree");
  
  ofstream outfile;
  outfile.open("junk.dat");
	
  for (Int_t id=2; id<4; id++) {
    for (Int_t chan=0; chan<32; chan++) {
    
      if(id==3 && chan>15) continue;

      c1->cd(1);
      DataTree->Draw("ADC.Data>>h1",Form("ADC.ID==%d && ADC.ChNum==%d && ADC.Data>180 && ADC.Data<4900",id,chan));
      if(h1->GetEntries()==0) continue;
      
      TSpectrum *s = new TSpectrum(npeaks+1);
      h1->SetTitle(Form("ADC %i Chan %i",id,chan));
      
      
      Int_t nfound = s->Search(h1,3," ",0.05); //Search(histo,sigma,option,threshold)
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
      c1->cd(2);
      if (nfound==6){
	TGraph *graph = new TGraph(npeaks,xpeaks,Volts);
	graph->Draw("A*");
      }else{
	TGraph *graph = new TGraph(4,xpeaks,Volts2);
	graph->Draw("A*");
      }
      TF1 *fun = new TF1("fun","[0]+[1]*x",0,5000);
      graph->Fit(fun);
      c1->Update();
      c1->WaitPrimitive();
  
      outfile << id << "\t" << chan << "\t" << fun->GetParameter(0) << "\t" << fun->GetParameter(1) << endl;

      //TGraph *FitGraph = new TGraph(24,xpeaks,&(Volts[0]));
    }
  }
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
