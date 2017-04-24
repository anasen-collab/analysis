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

  TCanvas *c1 = new TCanvas();
  c1->Divide(1,2);

  //TFile *f1 = new TFile("/data0/nabin/ANASEN/ANASEN_NKJ/ANASEN_N/253Si_Alpha_Cal.root");
  //TFile *f1 = new TFile("/data0/nabin/ANASEN/ANASEN_NKJ/ANASEN_N/Calibration/253-60_Si_Alpha_Cal.root");
  TFile *f1 = new TFile("/home/lighthall/anasen/root/run1255-61.root");
  if ( !f1->IsOpen() ) {
    cout << "Error: Root file does not exist\n";
    exit(EXIT_FAILURE);
  }
  
  TH1I *h1 = new TH1I("h1","h1",3000,0,range);
  TTree *DataTree = (TTree*)f1->Get("DataTree");
  TF1 *fit = new TF1("fit","pol1",0,range);

  Double_t zeroshift = 0;
  Double_t MeVperCh = 0;
  Double_t q0 = 0;
  TSpectrum *s = 0;
  TGraph *FitGraph = 0;
  
  ofstream outfile;
  outfile.open("saves/SiAlpha.dat");
  
  for (Int_t MBID=1; MBID<3; MBID++) {  
    for (Int_t CBID=1; CBID<15; CBID++) { 
      ////for the front of the detectors //run250 || 267
      //if((CBID==1 || CBID==2 || CBID==5 || CBID==6 || CBID==9 ||CBID==10) ||(MBID==1 && (CBID==11 || CBID==12 || CBID==13 || CBID==14))) {

      ////for the back of the detectors //run251 || 268
      //if((CBID==3 || CBID==4 || CBID==7 || CBID==8) || (MBID==2 && (CBID==11 || CBID==12))) {	 
      //if(MBID==2 && (CBID==3 || CBID==4 || CBID==7 || CBID==8 || CBID==11 || CBID==12)) {	 
	
	for (Int_t ChNum=0; ChNum<16; ChNum++) {
	  //// Mask bad channels  //that can create problem in cruising calibration.

	  if (MBID==1 && CBID==1 && (ChNum==4 ||ChNum ==14)) { //bad at front 
	    continue;
	  }	      
	  if((MBID==1 && CBID==8 && (ChNum==0 || ChNum==1 || ChNum==2))||(MBID==2 && CBID ==8 && (ChNum==1 || ChNum==2))) {//bad at back 
	    continue;
	  } 
	  if(MBID == 2 && CBID == 7 && (ChNum ==4 || ChNum ==10)) {//bad at back 
	    continue;
	  }
		
	  c1->cd(1);	    
	  DataTree->Draw("Si.Energy>>h1",Form("Si.MBID==%d && Si.CBID==%d && Si.ChNum==%d",MBID,CBID,ChNum));
	 
	  //c1->WaitPrimitive();
	  if(h1->GetEntries()==0) continue;
	      
	  if(s!=0) {
	    delete s;
	  }

	  TSpectrum *s = new TSpectrum();
	       
	  h1->SetTitle(Form("MBID %d CBID %d ChNum %d",MBID,CBID,ChNum));  
	  h1->GetXaxis()->SetRangeUser(1500,6000);
	  c1->Update();
	  
	  //Int_t nfound = s->Search(h1,15,"",0.25);
	  //Int_t nfound = s->Search(h1,5," nobackground",0.10);
	  Int_t nfound = s->Search(h1,10," ",0.05);

	  Float_t *xpeaks = s->GetPositionX();
	  Float_t Temp=0;

	  if(nfound <npeaks) {
	    printf("Less than %d peaks found. Aborting.\n",npeaks);
	    continue;
	  }
	      
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

	  TGraph *FitGraph = new TGraph(nfound,xpeaks, &(Energies[0]));
	  c1->cd(2);
	      
	  FitGraph->Draw("AP*");
	  FitGraph->Fit("fit","qROB=0.95");
	  ////FitGraph->Fit("fit","E");
	  zeroshift = fit->GetParameter(0);
	  MeVperCh = fit->GetParameter(1);
	  q0 = -zeroshift/MeVperCh;
	  // cout << zeroshift << " " << MeVperCh << endl;;
	  cout << MBID << " " << CBID << " " << ChNum << " q0 = "<<q0<<endl;
	  c1->Update();
	  outfile << MBID << "\t" << CBID << "\t" << ChNum << "\t"<< zeroshift << "\t" << MeVperCh <<endl;
	  //outfile << MBID << "\t" << zeroshift << "\t" << MeVperCh <<endl;
	  //outfile << MBID << "\t" << CBID << "\t" << ChNum << "\t"<< zeroshift << "\t" << MeVperCh << "\t"<< q0 <<endl;
	}
	//}
    }
  }
  delete FitGraph;
  delete s;
}
