//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Relative calibration of Si gains
//
// Output file (e.g."Sipulser_2015Dec13.dat") has the following columns:
// MBID, CBID, ASICs_Channel, ZeroShift(offset), Voltage_per_Ch(slope)
//
// Usage: root -l SiPulser_All.C++ (from the same directory).
//
// Edited by : John Parker , 2016Jan22
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//C++
#include <fstream>
#include <exception>
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
#include <TLine.h>
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
using namespace std;

void AlphaCal(void) {
  const Int_t npeaks = 3;
  const Int_t range=4096*4;
  
  //for gold at spacer zero, this is what the alpha energy should be
  //see lab book page 9
  Double_t energy[28];
  for (Int_t i=0; i<28; i++) {
    if (i==0 || i==3) {
      energy[i] = 12.03;
    }else if (i==1 || i==2) {
      energy[i] = 12.07;
    }else if (i>3 && i<16) {
      energy[i] = 12.3;
    }else if (i>15) {
      energy[i] = 12.84;
    }
  }

  Float_t Energies[npeaks] = {5,7,10};
  Float_t Energies1[npeaks-1] = {5,10};

  //TFile *f1 = new TFile("/home2/parker/ANASEN/LSU/ParkerMain_root/run_alpha0_282_284_cal022716.root");
  //TFile *f1 = new TFile("/home2/parker/ANASEN/LSU/ParkerMain_root/run417_cal.root");
  //TFile *f1 = new TFile("/home/lighthall/anasen/root/run1255-7mQ2S3_geo_init.root");
  TFile *f1 = new TFile("/home/lighthall/anasen/root/run1255-61mQ2S3_geo_init.root");
  
  if ( !f1->IsOpen() ) {
    cout << "Error: Root file does not exist\n";
    exit(EXIT_FAILURE);
  }
    
  TTree *MainTree = NULL;
  MainTree = (TTree*)f1->Get("MainTree");
  if (MainTree==NULL) {
    cout << "Tree does not exist\n";
    exit(EXIT_FAILURE);
  }
  
  ifstream infile;
  infile.open("saves/AlphaCalibration_init.dat");
  Int_t det=0,ch=0;
  Double_t dummy_slope = 0;
  Double_t slope[27];
  if (infile.is_open()) {
    infile.ignore(100,'\n');//read in dummy line
    while (!infile.eof()) {
      infile >> det >> ch >> dummy_slope;
      slope[det] = dummy_slope;
    }
  }
  else {
    cout << "Infile not opened\n";
    exit(EXIT_FAILURE);
  }
  infile.close();

  ofstream outfile;
  outfile.open("saves/AlphaCalibration.dat");
  outfile << "DetNum\tOffset\tSlope\n";
  
  TCanvas *can = new TCanvas("can","can",800,600);
  can->Divide(1,2);
  TH1F *hist = new TH1F("hist","hist",300,0,range);
  TF1 *fit = new TF1("fit","pol1",0,range);
  Double_t zeroshift = 0;
  Double_t MeVperCh = 0;
  Double_t q0 = 0;
  TSpectrum *s = 0;
  TGraph *FitGraph = 0;

  for (Int_t DetNum=0; DetNum<28; DetNum++) {
    if (DetNum==1 || DetNum==2) {
      Energies[0] = 5;
      Energies[1] = 7;
      Energies[2] = 10;
    }
    else if (DetNum>3 && DetNum<16) {
      Energies[0] = 5;
      Energies[1] = 7;
      Energies[2] = 10;    
    }
    else if (DetNum>15) {
      Energies[0] = 5;
      Energies[1] = 7;
      Energies[2] = 10;    
    }
    
    // if ( DetNum==7 || DetNum==9 ) {
    //   continue;
    // }
    Float_t *average_slope;
    can->cd(1);	 
    MainTree->Draw("EnergyBack>>hist",Form("DetID==%i && (HitType==111 || HitType==11)",DetNum),"");
    hist->SetTitle(Form("Det%i",DetNum));  

    if(hist->GetEntries()==0) {
      printf("Histogram %s has zero entries.\n",hist->GetTitle());
      continue;
    }
	      
    if(s!=0) {
      delete s;
    }

    hist->GetXaxis()->SetRangeUser(1500,6000);
    can->Update();
    TSpectrum *s = new TSpectrum();
    Int_t nfound = s->Search(hist,2," ",0.2);//9 and 0.15
   
    if(nfound <(npeaks-1)) {
      printf("DetNum %2d: peaks = %d, ",DetNum,nfound);
      printf("Less than %d peaks found. Aborting.\n",npeaks);
      outfile << DetNum << "\t"<< 0 << "\t" << 1 <<endl;
      continue;
    }

    if(nfound >npeaks) {
      printf("DetNum %2d: peaks = %d, ",DetNum,nfound);
      printf("greater than %d peaks found. Aborting.\n",npeaks);
      outfile << DetNum << "\t"<< 0 << "\t" << 1 <<endl;
      continue;
    }
    
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
    
    if (npeaks==1) {
      if (nfound != 1) {
	cout << DetNum << "  " << nfound << endl;
	cout << "No peaks found\n";
	//Double_t median = hist->GetMaximumBin()*0.00266666 + 0.8;
	Double_t median = hist->GetMaximumBin()*0.0012 + 0.8;
	cout << median << endl;
	slope[DetNum] *= (energy[DetNum]/median);
	TLine *line = new TLine(median,0,median,1000);
	line->Draw("same");
	can->Update();
	continue;
      }
      average_slope = s->GetPositionX();
      slope[DetNum] *= (energy[DetNum]/average_slope[0]);

      cout << DetNum << " " << average_slope[0] << "  " << slope[DetNum] << endl;
      outfile << DetNum << "\t" << 0 << "\t" << slope[DetNum] << endl;
      can->Update();
    }

    if(FitGraph!=0) {
      delete FitGraph;
    }

    if(nfound==npeaks) {
      FitGraph = new TGraph(nfound,xpeaks, &(Energies[0]));
    }
    else if(nfound==(npeaks-1)) {
      FitGraph = new TGraph(nfound,xpeaks, &(Energies1[0]));
    }
    can->cd(2);
    FitGraph->Draw("AP*");
    FitGraph->Fit("fit","q");
    zeroshift = fit->GetParameter(0);
    MeVperCh = fit->GetParameter(1);
    q0 = -zeroshift/MeVperCh;
    // cout << zeroshift << " " << MeVperCh << endl;;
    printf("DetNum %2d: peaks = %d, slope = %f /t offset = %f\n",DetNum,nfound,zeroshift,MeVperCh );
    can->Update();
    outfile << DetNum << "\t"<< zeroshift << "\t" << MeVperCh <<endl;
  }
  delete FitGraph;
  delete s;
}
