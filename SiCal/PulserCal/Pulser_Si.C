#define Pulser_Si_cxx
/*********************************************************************************************
Code: Pulser_Si.C

Description: This code takes pulser data for the silicon detectors and gets the offsets
             and slopes for each channel in each detector by fitting the position of the
             peaks to the respective voltages.  These voltages must be provided by the user
             below in the Volts array, as well as the dimension of this array, npeaks, which
             is stands for the number of expected peaks in the histograms.

*********************************************************************************************/


#include "Pulser_Si.h"
#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TSpectrum.h>
#include <TGraph.h>
#include <TF1.h>
#include <iostream>
#include <fstream>

using namespace std;

#include "../Main/ChannelMap.h"


ChannelMap *CMAP;
TCanvas* Mon;
TH1I*** PulserHist;
TGraph*** FitGraph;

Int_t MBID[30][12];
Int_t CBID[30][12];
Int_t ChNum[30][12];

/////////////////////////////////////////////////////////////////////////////////////////////
// CONTROL PANEL
const Int_t npeaks = 6;
const Float_t Volts[npeaks] = {0.5,1,3,5,7,9};
const Float_t PeakSigma = 10.0;        
const Float_t MinPeakAmplitude = 0.05; // The treshold peak amplitud to be used in the Search() method.
const Float_t NormalChi2Value = 0.003; //
/////////////////////////////////////////////////////////////////////////////////////////////

void Pulser_Si::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

   cout << "----------------------------------------------------------------\n"
	<< "> Pulser_Si\n";


   CMAP = new ChannelMap();
   CMAP->Init("../Main/ASICS_cmap_Final2015");
   

   FitGraph = new TGraph**[24];
   PulserHist = new TH1I**[24];
   for(Int_t det=4; det<28; det++){
     FitGraph[det-4] = new TGraph*[12];
     PulserHist[det-4] = new TH1I*[12];
     for(Int_t ch=0; ch<12; ch++){
       PulserHist[det-4][ch] = new TH1I(Form("PulserHist_SX3%d_ch%d",det,ch),Form("SX3=%d ch=%d",det,ch),4096,100,16384);
       PulserHist[det-4][ch]->GetXaxis()->SetTitle("SX3 Voltage [arb. units]");     
     }
   }
   cout << "> Matrix of 30 by 12 histograms created.\n"
	<< "> Filling histograms ..." << endl;
   


   for(Int_t MB=1; MB<=2; MB++){
     for(Int_t CB=1; CB<=14; CB++){
       cout << MB << "  " << CB << endl;
       for(Int_t Ch=0; Ch<16; Ch++){
	 Int_t DetID=-1, DetCh=-1;
	 CMAP->IdentifyDetChan(MB,CB,Ch, DetID, DetCh);
	 if(DetID>=0 && DetCh>=0){
	   MBID[DetID][DetCh] = MB;
	   CBID[DetID][DetCh] = CB;
	   ChNum[DetID][DetCh] = Ch;
	 }
       }
     }
   }
   cout << "All good here\n";
}

void Pulser_Si::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

}

Bool_t Pulser_Si::Process(Long64_t entry)
{
   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // It can be passed to either Pulser_Si::GetEntry() or TBranch::GetEntry()
   // to read either all or the required parts of the data. When processing
   // keyed objects with PROOF, the object is already loaded and is available
   // via the fObject pointer.
   //
   // This function should contain the "body" of the analysis. It can contain
   // simple or elaborate selection criteria, run algorithms on the data
   // of the event and typically fill histograms.
   //
   // The processing can be stopped by calling Abort().
   //
   // Use fStatus to set the return value of TTree::Process().
   //
   // The return value is currently not used.
  
  if (entry == TMath::Nint(0.01*fChain->GetEntries())) cout << ">  1% through the data" << endl;
  if (entry == TMath::Nint(0.1*fChain->GetEntries()))  cout << "> 10% through the data" << endl;
  if (entry == TMath::Nint(0.2*fChain->GetEntries()))  cout << "> 20% through the data" << endl;
  if (entry == TMath::Nint(0.3*fChain->GetEntries()))  cout << "> 30% through the data" << endl;
  if (entry == TMath::Nint(0.4*fChain->GetEntries()))  cout << "> 40% through the data" << endl;
  if (entry == TMath::Nint(0.5*fChain->GetEntries()))  cout << "> 50% through the data" << endl;
  if (entry == TMath::Nint(0.6*fChain->GetEntries()))  cout << "> 60% through the data" << endl;
  if (entry == TMath::Nint(0.7*fChain->GetEntries()))  cout << "> 70% through the data" << endl;
  if (entry == TMath::Nint(0.8*fChain->GetEntries()))  cout << "> 80% through the data" << endl;
  if (entry == TMath::Nint(0.9*fChain->GetEntries()))  cout << "> 90% through the data" << endl;
  if (entry == TMath::Nint(1.0*fChain->GetEntries()-1))cout << "> Complete" << endl;

  // Getting data from branches   
  b_SiNhits->GetEntry(entry);
  b_Si_MBID->GetEntry(entry);
  b_Si_CBID->GetEntry(entry);
  b_Si_ChNum->GetEntry(entry);
  b_Si_Energy->GetEntry(entry);
  b_Si_Time->GetEntry(entry);
 
  Int_t DetID=-1,DetCh=-1;

  for (Int_t n=0; n<Si_Nhits; n++) 
  {
	  CMAP->IdentifyDetChan(Si_MBID[n],Si_CBID[n],Si_ChNum[n], DetID, DetCh);
	  PulserHist[DetID][DetCh]->Fill(Si_Energy[n]);
  }

   return kTRUE;
}

void Pulser_Si::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}

void Pulser_Si::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.

  Mon = new TCanvas("Mon","",1300,600);
  Mon->Divide(2);
  TF1* linear;
  Float_t chi2[28][12];
  Float_t slope[28][12];
  Float_t offset[28][12];

  for(Int_t det=4; det<28; det++)
  {    
	  for(Int_t ch=0; ch<12; ch++)
	  {
		  chi2[det-4][ch]=999;
		  slope[det-4][ch]=1;
		  offset[det-4][ch]=0;

		  Mon->cd(1);
		  TSpectrum *s = new TSpectrum(19);
		  Int_t nfound = s->Search(PulserHist[det-4][ch],PeakSigma," ",MinPeakAmplitude);
		  PulserHist[det-4][ch]->SetTitle(Form("Peaks found = %d",nfound));
		  Mon->Update();

		  cout << det << "\t" << ch << "\t" << nfound << "\t";

		  Float_t *xpeaks = s->GetPositionX();
		  Float_t Temp=0;
		  
		  //Reorder the peaks' x-coordinate from lower to higher x-value.
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

		  // This section depends on how your data looks like. 
		  // We saw three cases: nfound=17 (for the back segments), nfound=13 (for front segments)
		  // and nfound=12 (for det=6 ch=4).
		  if(nfound==0) {cout << "Zero Peaks found for Det " << det << " channel " << ch << endl;}
		  if(nfound>npeaks){nfound=npeaks;}
		  FitGraph[det][ch] = new TGraph(nfound,&xpeaks[0],&Volts[0]);

      
		  FitGraph[det][ch]->GetXaxis()->SetTitle("Measured signal [arb. units]");
		  FitGraph[det][ch]->GetXaxis()->SetLimits(0,16000);
		  FitGraph[det][ch]->GetYaxis()->SetTitle("Pulser voltage [V]");
		  FitGraph[det][ch]->GetYaxis()->SetLimits(0,10);
      
		  linear = new TF1("linear","pol1",100,16000);
		  FitGraph[det][ch]->Fit(linear,"qr");
		  chi2[det][ch] = linear->GetChisquare();
		  slope[det][ch] = linear->GetParameter(1);
		  offset[det][ch] = linear->GetParameter(0);
		  cout << offset[det][ch] << "\t" << slope[det][ch] << "\t" << chi2[det][ch] << "\t" << endl;
		  
		  FitGraph[det][ch]->SetTitle(Form("Det=%d  Ch=%d  |  Chi^2=%f",det, ch, chi2[det][ch]));
      
		  Mon->cd(2);
		  FitGraph[det][ch]->Draw("AP*");
		  Mon->Update();
		  Mon->WaitPrimitive();
	  }
  }

  ofstream outputfile;
  outputfile.open("alignchannels.100912_2");
  

  for(Int_t MB=1; MB<=2; MB++)
  {
	  for(Int_t CB=1; CB<=14; CB++)
	  {
		  for(Int_t Ch=0; Ch<16; Ch++)
		  {
			  Int_t DetID=-1, DetCh=-1;
			  CMAP->IdentifyDetChan(MB,CB,Ch, DetID, DetCh);
			  if(DetID>=0 && DetCh>=0) 
			  {
				  outputfile << MB << "\t" << CB << "\t" << Ch << "\t";
				  outputfile << offset[DetID][DetCh] << "\t" << slope[DetID][DetCh] << endl;
				  if(chi2[DetID][DetCh]<NormalChi2Value)
				  {
				  }
				  else 
				  {
					  cout << "> Warning: Det=" << DetID << " Ch=" << DetCh;
					  cout << " High chi^2 value " << chi2[DetID][DetCh] << endl;
				  }
			  }
		  }	    
      }
  }
  outputfile.close();

 
}
