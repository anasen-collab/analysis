//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//// ASICS all channels Pulser calibration.
////
//// Output file (e.g."Sipulser_2015Dec13.dat") has the following columns:
//// MBID, CBID, ASICs_Channel, ZeroShift(offset), Voltage_per_Ch(slope)
////
//// Usage: root -l SiPulser_All.C++ (from the same directory).
////
//// Edited by : Nabin Rijal , 2015Dec13
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
void SiPulser_All(void)
{
  //MBID CID ASICs_Chn ZeroShift VperCh

  //run250//run251 //
  const Int_t npeaks = 6;
  Float_t Volts[npeaks] = { 0.5, 1.0, 3.0, 5.0, 7.0, 9.0 };
  Float_t Volts5[5] = { 1.0, 3.0, 5.0, 7.0, 9.0 };
  Float_t Volts4[4] = { 3.0, 5.0, 7.0, 9.0 };
  
  TCanvas *c1 = new TCanvas();
  c1->Divide(1,2);

  TFile *f1 = new TFile("run585.root");//front
  //TFile *f1 = new TFile("/data0/nabin/ANASEN/ANASEN_NKJ/New/evt2root/run251_NSCL11_Pulser.root");//back

  TH1I *h1 = new TH1I("h1","h1",16084,300,16384);
  TTree *DataTree = (TTree*)f1->Get("DataTree");
  TF1 *fit = new TF1("fit","pol1",0,8192);

  Double_t zeroshift = 0;
  Double_t vperch = 0;
  Double_t q0 = 0;
  TSpectrum *s = 0;
  TGraph *FitGraph = 0;
  
  ofstream outfile;
  ofstream outfile2;

  outfile.open("junk.dat");
  outfile2.open("junk2.dat");  

  for (Int_t MBID=1; MBID<3; MBID++)
    {  
      for (Int_t CBID=1; CBID<15; CBID++) 
	{ 
	  ////for the front of the detectors //run250
	  //if((CBID==1 || CBID==2 || CBID==5 || CBID==6 || CBID==9 ||CBID==10) ||(MBID==1 && (CBID==11 || CBID==12 || CBID==13 || CBID==14))){
	 

	    ////for the back of the detectors //run251
	    if((CBID==3 || CBID==4 || CBID==7 || CBID==8) || (MBID==2 && (CBID==11 || CBID==12))){	 
	
	
	    for (Int_t ChNum=0; ChNum<16; ChNum++) 	
	      {
	     
	      
		//// Mask bad channels  //that can create problem in cruising calibration.

		if (MBID==1 && CBID==1 && (ChNum==4 ||ChNum ==14)){ //bad at front run250
		  continue;
		}	      
		if((MBID==1 && CBID==8 && (ChNum==0 || ChNum==1 || ChNum==2))||(MBID==2 && CBID ==8 && (ChNum==1 || ChNum==2))){//bad at back run251
		  continue;
		} 
		if(MBID == 2 && CBID == 7 && (ChNum ==10)){//bad at back run251
		  continue;
		}
	      
	      
		c1->cd(1);	    
		DataTree->Draw("Si.Energy>>h1",Form("Si.MBID==%d && Si.CBID==%d && Si.ChNum==%d",MBID,CBID,ChNum));
		c1->Update();
		c1->WaitPrimitive();
		if(h1->GetEntries()==0) continue;
	      
		if(s!=0)
		  {
		    delete s;
		  }

		TSpectrum *s = new TSpectrum(npeaks+1);
		
	       
		h1->SetTitle(Form("MBID %d CBID %d ChNum %d",MBID,CBID,ChNum));  

		//Int_t nfound = s->Search(h1,15,"",0.25);
		//Int_t nfound = s->Search(h1,5," nobackground",0.10);
		Int_t nfound = s->Search(h1,20," ",0.05);//9 and 0.15

		Float_t *xpeaks = s->GetPositionX();
		Float_t Temp=0;
	      
		for(Int_t i=0;i<nfound;i++) 
		  {
		    for(Int_t j=i;j<nfound;j++) 
		      {
			if (xpeaks[j] < xpeaks[i]) 
			  {
			    Temp = xpeaks[i];
			    xpeaks[i] = xpeaks[j];
			    xpeaks[j] = Temp;  
			  }
		      }
		  }

	    
		if(FitGraph!=0)
		  {
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
    }
    
  delete FitGraph;
  delete s;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////		 

	    
	     
    
