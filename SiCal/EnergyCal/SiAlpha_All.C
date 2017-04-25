//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//// ASICS all channels Alpha calibration.
////
//// Output file (e.g."SiAlpha_2016Jan22.dat") has the following columns:
//// MBID, CBID, ASICs_Channel, ZeroShift(offset), MeV/Ch(slope)
////
//// Usage: root -l SiAlpha_All.C++ (from the same directory).
////
//// Edited by : Nabin Rijal , 2016Jan22
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
void SiAlpha_All(void) {
  //MBID CID ASICs_Chn ZeroShift MeVperCh
  const Int_t npeaks = 3;
  const Int_t range=4096*4;

  //Float_t Energies[npeaks] = {5.42315,5.68537,6.05,6.28808,6.77803,8.74886};
  //Float_t Energies[npeaks] = {5.42315,5.68537,6.28808,6.77803,8.74886};
  Float_t Energies[npeaks] = {5,7,10};
  Float_t Energies1[npeaks] = {5,10};

  //TFile *f1 = new TFile("/data0/nabin/ANASEN/ANASEN_NKJ/ANASEN_N/253Si_Alpha_Cal.root");
  //TFile *f1 = new TFile("/data0/nabin/ANASEN/ANASEN_NKJ/ANASEN_N/Calibration/253-60_Si_Alpha_Cal.root");
  //TFile *f1 = new TFile("/home/lighthall/anasen/root/run1255-61.root");
  TFile *f1 = new TFile("/home/lighthall/anasen/root/run1255-61mQ2S3_geo_init.root");
  
  if ( !f1->IsOpen() ) {
    cout << "Error: Root file does not exist\n";
    exit(EXIT_FAILURE);
  }
    
  TTree *MainTree = NULL;
  MainTree = (TTree*)f1->Get("MainTree");
  if (MainTree==NULL){
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
  outfile.open("saves/SiAlpha.dat");
  outfile << "DetNum\tOffset\tSlope\n";

  TCanvas *can = new TCanvas();
  can->Divide(1,2);
  TH1I *hist = new TH1I("hist","hist",3000,0,range);
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
    
    can->cd(1);	    
    MainTree->Draw("EnergyBack>>hist",Form("DetID==%i && (HitType==111 || HitType==11)",DetNum),"");
    hist->SetTitle(Form("Det%i",DetNum));  
	  
    //can->WaitPrimitive();
    if(hist->GetEntries()==0) {
      printf("Histogram %s has zero entries.\n",hist->GetTitle());
      continue;
    }
	      
    if(s!=0) {
      delete s;
    }

    TSpectrum *s = new TSpectrum();
	  
    hist->GetXaxis()->SetRangeUser(1500,6000);
    can->Update();
	  
    //Int_t nfound = s->Search(hist,15,"",0.25);
    //Int_t nfound = s->Search(hist,5," nobackground",0.10);
    Int_t nfound = s->Search(hist,10," ",0.05);
    Float_t *xpeaks = s->GetPositionX();

    if(nfound <(npeaks-1)) {
      printf("DetNum %2d: peaks = %d, ",DetNum,nfound);
      printf("Less than %d peaks found. Aborting.\n",npeaks);
      outfile << DetNum << "\t"<< 0 << "\t" << 1 <<endl;
      continue;
    }

    //sort peaks in order of channle number
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
    FitGraph->Fit("fit","qROB=1");
    ////FitGraph->Fit("fit","E");
    zeroshift = fit->GetParameter(0);
    MeVperCh = fit->GetParameter(1);
    q0 = -zeroshift/MeVperCh;
    // cout << zeroshift << " " << MeVperCh << endl;;
    printf("DetNum %2d: peaks = %d, slope = %f /t offset = %f\n",DetNum,nfound,zeroshift,MeVperCh );
    //can->Update();
    outfile << DetNum << "\t"<< zeroshift << "\t" << MeVperCh <<endl;
    
  }
  delete FitGraph;
  delete s;
}
