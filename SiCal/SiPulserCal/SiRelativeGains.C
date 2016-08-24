//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//// Relative calibration of Si gains
////
//// Output file (e.g."Sipulser_2015Dec13.dat") has the following columns:
//// MBID, CBID, ASICs_Channel, ZeroShift(offset), Voltage_per_Ch(slope)
////
//// Usage: root -l SiPulser_All.C++ (from the same directory).
////
//// Edited by : John Parker , 2016Jan22
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
void SiRelativeGains(void)
{
  //MBID CID ASICs_Chn ZeroShift VperCh

  //run250//run251 //
  const Int_t npeaks = 4;
  Float_t Volts[npeaks] = { 0.5, 1.0, 3.0, 5.0, 7.0, 9.0 };
  
  TCanvas *c1 = new TCanvas();
  c1->Divide(2,2);

  TFile *f1 = new TFile("run236.root");//front
  //TFile *f1 = new TFile("/data0/nabin/ANASEN/ANASEN_NKJ/New/evt2root/run251_NSCL11_Pulser.root");//back

  TH1I *h1 = new TH1I("h1","h1",16084,300,16384);
  //TTree *DataTree = (TTree*)f1->Get("DataTree");
  TF1 *fit = new TF1("fit","pol1",0,8192);

  Double_t zeroshift = 0;
  Double_t vperch = 0;
  Double_t q0 = 0;
  TSpectrum *s = 0;
  TGraph *FitGraph = 0;
  
  ofstream outfile;
  ofstream outfile2;

  outfile.open("X3RelativeGains012216.dat");
  //outfile2.open("Sipulser_peaklocations22Jan2016Back2.dat");
  TF1 *fun = new TF1("fun","1 - x",0,1);

  for (Int_t DetNum=4; DetNum<28; DetNum++){  
    for (Int_t FrontChNum=0; FrontChNum<4; FrontChNum++){

      for (Int_t BackChNum=0; BackChNum<4; BackChNum++){
	TH2F *hist = (TH2F*)f1->Get(Form("up_vs_down%i_%i_%i",DetNum,FrontChNum,BackChNum));
	hist->GetXaxis()->SetRange(0,120);
	hist->GetYaxis()->SetRange(0,120);
      
	//// Mask bad channels  //that can create problem in cruising calibration.
	
	c1->cd(1);
	hist->ProjectionY("cut2",18,22);
	hist->ProjectionY("cut3",
	cut2->Draw();
	
	c1->Update();

	c1->WaitPrimitive();
	if(hist->GetEntries()==0) continue;
	
	if(s!=0){
	  delete s;
	}
      
	TSpectrum *s = new TSpectrum(1);

	Int_t nfound = s->Search(h1,5," ",0.9);//9 and 0.15
	
	Float_t *xpeaks = s->GetPositionX();
	Float_t Temp=0;
	
	if(FitGraph!=0){
	  delete FitGraph;
	}
	if (nfound==6){
	  FitGraph = new TGraph(nfound,xpeaks, &(Volts[0]));
	}else if (nfound==5){
	  FitGraph = new TGraph(5,xpeaks, &(Volts5[0]));
	}else if (nfound==4){
	  FitGraph = new TGraph(4,xpeaks, &(Volts4[0]));
	}else if (nfound==7){
	  FitGraph = new TGraph(6,xpeaks, &(Volts[0]));
	}else{
	  cout << "Wrong number of peaks\n";
	  FitGraph = new TGraph(nfound,xpeaks, &(Volts[0]));
	}
	c1->cd(2);
	
	FitGraph->Draw("AP*");
	FitGraph->Fit("fit","qROB=0.95");
	////FitGraph->Fit("fit","E");
	zeroshift = fit->GetParameter(0);
	vperch = fit->GetParameter(1);
	q0 = -zeroshift/vperch;
	// cout << zeroshift << " " << vperch << endl;;
	cout << MBID << " " << CBID << " " << ChNum << " q0 = "<<q0<<endl;
	c1->Update();
	outfile << MBID << "\t" << CBID << "\t" << ChNum << "\t"<< zeroshift << "\t" << vperch <<endl;
	//outfile << MBID << "\t" << CBID << "\t" << ChNum << "\t"<< zeroshift << "\t" << vperch << "\t"<< q0 <<endl;
	outfile2 << MBID << "\t" << CBID << "\t" << ChNum;
	for (int i=0; i<nfound; i++){
	  outfile2 << "\t" << xpeaks[i];
	}
	outfile2 << endl;
      }
    }
  }

  delete FitGraph;
  delete s;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////		 

	    
	     
    
