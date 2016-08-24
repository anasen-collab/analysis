//Written by John Parker
//This code matches Silicon events with PC events based on the phi location of the event
//Usage: g++ -o ParkerTrack EnergyLoss.cpp ParkerTrack.cpp `root-config --cflags --glibs`
//      ./ParkerROOT DataList.txt outputfile.root cutfile_1 cutfile_2
////DataList.txt contains a list of rootfile from ParkerMain
////------------------------------------------------------
////If you are sorting lots of files and don't want to make a huge data tree, you can comment out the #define FillTree line below and only the histograms will be filled (or vice versa)

//#define TimingCut
//#define EnergyCut

#define M_PI  3.14159265358979323846264338328 // Pi 

#define FillHists
#define FillTree

#define MaxSiHits   500
#define MaxADCHits  500
#define MaxTDCHits  500

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <string>
#include <sstream>
#include <exception>
#include <stdexcept>
#include <map>
#include <algorithm>

//ROOT libraries
#include <TFile.h>
#include <TTree.h>
#include <TROOT.h>
#include <TNtuple.h>
#include <TCanvas.h>
#include <TRint.h>
#include <TObjArray.h>
#include <TGraph.h>
#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TMath.h>
#include <TList.h>
#include <TCutG.h>
#include <TLorentzVector.h>
#include <TVector3.h>

#include "/home/manasta/Desktop/parker_codes/Include/organizetree.h"
#include "/home/manasta/Desktop/parker_codes/Include/Reconstruct.h"

using namespace std;

void MyFill(string name,
	    int binsX, double lowX, double highX, double valueX);

void MyFill(string name,
	    int binsX, double lowX, double highX, double valueX,
	    int binsY, double lowY, double highY, double valueY);

TList* fhlist;
std::map<string,TH1*> fhmap;

bool Track::Si_sort_method(struct Silicon_Event a,struct Silicon_Event b){
  if(a.SiEnergy > b.SiEnergy)
    return 1;
  return 0;
};

bool Track::PC_sort_method(struct PropCounter_Event a,struct PropCounter_Event b){
  if(a.PCEnergy > b.PCEnergy)
    return 1;
  return 0;
};

int main(int argc, char* argv[]){
  TApplication *myapp=new TApplication("myapp",0,0); //Don't know what this does, but libraries won't load without it
  //call with 4 arguments :  

  if ( argc!=3){
    cout << "Error: Wrong Number of Arguments\n";
    exit(EXIT_FAILURE);
  }

  char* filename_raw  = new char [100]; // for input  .root filename
  char* filename_cal = new char [100]; // for output .root filename
  char* filename_cut = new char[100]; //for cut file
  char* filename_cut2 = new char[100]; //for cut file

  strcpy( filename_raw, argv[1] );
  strcpy( filename_cal, argv[2] );
  //strcpy( filename_cut, argv[3] );
  //strcpy( filename_cut2, argv[4] );

  
  TObjArray *RootObjects = new TObjArray();
  SiHit Si;
  PCHit PC;
  Track Old_Tr;
  Track Tr;
  Reconstruct Rec;
  Int_t RFTime,MCPTime;
  Int_t Old_RFTime,Old_MCPTime;
  Double_t IC;
  Double_t E_IC;
  Double_t Old_IC;
  Double_t Old_E_IC;
  

  Si.ReadDet = 0;
  Si.ReadHit = 0;
  PC.ReadHit = 0;
  Old_Tr.ReadTrEvent = 0;

  
  TFile *outputfile = new TFile(filename_cal,"RECREATE");
  TTree *MainTree = new TTree("MainTree","MainTree");
  MainTree->Branch("Tr.NTracks",&Tr.NTracks,"NTracks/I"); //sum of all 3 tracks
  MainTree->Branch("Tr.NTracks1",&Tr.NTracks1,"NTracks1/I"); // #of correlated Si and PC events
  MainTree->Branch("Tr.NTracks2",&Tr.NTracks2,"NTracks2/I"); // #of Si no PC
  MainTree->Branch("Tr.NTracks3",&Tr.NTracks3,"NTracks3/I"); // PC no Si
  MainTree->Branch("Tr.TrEvent",&Tr.TrEvent); 
  MainTree->Branch("Tr.SiEvent",&Tr.SiEvent);
  MainTree->Branch("Tr.PCEvent",&Tr.PCEvent);
  MainTree->Branch("Tr.BeEvent",&Tr.BeEvent);
  MainTree->Branch("Tr.AlEvent",&Tr.AlEvent);
  MainTree->Branch("RFTime",&RFTime,"RFTime/I");
  MainTree->Branch("MCPTime",&MCPTime,"MCPTime/I");
  MainTree->Branch("IC",&IC,"IC/D");
  MainTree->Branch("E_IC",&E_IC,"E_IC/D");

  TH2F *E_de = new TH2F("E_de","E_de",1000,-1,35,100,-0.1,1);
  TH2F *E_theta = new TH2F("E_theta","E_theta",600,-0.1,180,600,-1,35);

  RootObjects->Add(MainTree);
  RootObjects->Add(E_de);

  fhlist = new TList;
  RootObjects->Add(fhlist);
  Double_t avg_beam_energy = 0;

  ifstream inFileList;
  inFileList.open(filename_raw);
  if (!inFileList.is_open()){
    cout << "In File List not open\n";
    exit(EXIT_FAILURE);
  }
  string rootfile;
  char rootfile_char[100];

  

  while (!inFileList.eof()){//loop over all of the incoming root files----------------------------------------------------------------------------------------------
    getline(inFileList,rootfile);
    if (rootfile.empty()){
      //cout << "Root File does not exist: " << rootfile << endl;
      break;
    }
    strcpy(rootfile_char,rootfile.c_str());
    TFile *inputFile = new TFile(rootfile_char);//open root file and make sure it exists----------------------------------------------------------------------------
    if (!inputFile->IsOpen()){
      cout << "Root file: " << rootfile << " could not be opened.\n";
      continue;
    }
    cout << "Processing File: " << rootfile << endl;

    TTree *raw_tree = (TTree*) inputFile->Get("MainTree");

    raw_tree->SetBranchAddress("Si.NSiHits",&Si.NSiHits);
    raw_tree->SetBranchAddress("Si.Detector",&Si.ReadDet);
    raw_tree->SetBranchAddress("Si.Hit",&Si.ReadHit);
    raw_tree->SetBranchAddress("PC.NPCHits",&PC.NPCHits);
    raw_tree->SetBranchAddress("PC.Hit",&PC.ReadHit);
    raw_tree->SetBranchAddress("Tr.NTracks",&Old_Tr.NTracks);
    raw_tree->SetBranchAddress("Tr.TrackEvent",&Old_Tr.ReadTrEvent);
    raw_tree->SetBranchAddress("RFTime",&Old_RFTime);
    raw_tree->SetBranchAddress("MCPTime",&Old_MCPTime);
    raw_tree->SetBranchAddress("IC",&Old_IC);
    raw_tree->SetBranchAddress("E_IC",&Old_E_IC);
    

    
    Long64_t nentries = raw_tree->GetEntries();
    Int_t status;
    for (Long64_t i=0; i<nentries; i++){//loop over all events
      status = raw_tree->GetEvent(i);

      //cout << "Event_Number: " << i << endl;
      if (i == TMath::Nint(0.01*nentries))  cout << " 1% through the data" << endl;
      if (i == TMath::Nint(0.10*nentries))  cout << " 10% through the data" << endl;
      if (i == TMath::Nint(0.15*nentries))  cout << " 15% through the data" << endl;
      if (i == TMath::Nint(0.25*nentries))  cout << " 25% through the data" << endl;
      if (i == TMath::Nint(0.35*nentries))  cout << " 35% through the data" << endl;
      if (i == TMath::Nint(0.50*nentries))  cout << " 50% through the data" << endl;
      if (i == TMath::Nint(0.65*nentries))  cout << " 65% through the data" << endl;
      if (i == TMath::Nint(0.75*nentries))  cout << " 75% through the data" << endl;
      if (i == TMath::Nint(0.90*nentries))  cout << " 90% through the data" << endl;
      if (i == TMath::Nint(0.95*nentries))  cout << " 95% through the data" << endl;
      if (i == TMath::Nint(1.00*nentries))  cout << " 100% through the data" << endl;


//-------------------------------------
//------------Ion Chamber----------------------

      IC = Old_IC;
      E_IC = Old_E_IC;

      IC = 0; E_IC = 0;
      for(Int_t n=0; n<ADC.Nhits; n++)
	{ if(ADC.ID[n]==3 && ADC.ChNum[n]==24)
	    {
	      IC = ADC.Data[n];
	      //cout << IC << endl;
	    }
	}

      for(Int_t n=0; n<ADC.Nhits; n++)
	{ if(ADC.ID[n]==3 && ADC.ChNum[n]==28)
	    {
	      E_IC = ADC.Data[n];
	      //cout << IC << endl;
	    }
	}


      Double_t correct=0.99023;
      MCPTime = Old_MCPTime;
      RFTime = Old_RFTime;

      MyFill("Timing",400,-600,600,(MCPTime*correct-RFTime)%538);

#ifdef TimingCut   
      if ( (MCPTime*correct-RFTime)%538<68 || (MCPTime*correct-RFTime)%538>370 ){
	continue;
      }
      if ( (MCPTime*correct-RFTime)%538>100 && (MCPTime*correct-RFTime)%538<330 ){
	continue;
      }
#endif

      Tr.zeroTrack();
      //first, just copy everything from the old file into one tree
      //Track Type 1 has a silicon and pc together
      //Track Type 2 has only a silicon
      //Track Type 3 has only pc events
      for (Int_t j=0; j<Old_Tr.ReadTrEvent->size(); j++){//loop over all of the tracks	
	//if we have a good hit type set the parameters in your new tree
	Old_Tr.track_place_holder = Old_Tr.ReadTrEvent->at(j);
	if ( Old_Tr.track_place_holder.SiEnergy <= 0 ){//make sure that the track was filled
	  continue;
	}
	//copy stuff from the old file
	Tr.track_place_holder.TrackType = 1;
	Tr.track_place_holder.DetID = Old_Tr.track_place_holder.DetID;
	Tr.track_place_holder.WireID = Old_Tr.track_place_holder.WireID;
	Tr.track_place_holder.SiEnergy = Old_Tr.track_place_holder.SiEnergy;
	Tr.track_place_holder.SiZ = Old_Tr.track_place_holder.SiZ;
	Tr.track_place_holder.SiR = Old_Tr.track_place_holder.SiR;
	Tr.track_place_holder.SiPhi = Old_Tr.track_place_holder.SiPhi;
	Tr.track_place_holder.PCEnergy = Old_Tr.track_place_holder.PCEnergy;
	Tr.track_place_holder.PCZ = Old_Tr.track_place_holder.PCZ;
	Tr.track_place_holder.PCR = Old_Tr.track_place_holder.PCR;
	Tr.track_place_holder.PCPhi = Old_Tr.track_place_holder.PCPhi;
	Tr.track_place_holder.IntPoint = Old_Tr.track_place_holder.IntPoint;
	Tr.track_place_holder.PathLength = Old_Tr.track_place_holder.PathLength;
	Tr.track_place_holder.Theta = Old_Tr.track_place_holder.Theta;
	Tr.track_place_holder.Phi = Old_Tr.track_place_holder.SiPhi;
	

#ifdef FillHists
	for(Tr.track_place_holder.DetID=0; Tr.track_place_holder.DetID<4; Tr.track_place_holder.DetID++){
	  cout << "DetID: " <<Tr.track_place_holder.DetID << endl;
	E_de->Fill(Tr.track_place_holder.SiEnergy,Tr.track_place_holder.PCEnergy);
	MyFill("E_de_corrected",300,0,30,Tr.track_place_holder.SiEnergy,300,0,1,Tr.track_place_holder.PCEnergy*sin(Tr.track_place_holder.Theta));

	E_theta->Fill(Tr.track_place_holder.Theta*180/M_PI,Tr.track_place_holder.SiEnergy);
	MyFill("E_theta",300,0,180,Tr.track_place_holder.Theta*180/M_PI,300,0,35,Tr.track_place_holder.SiEnergy);
	}
#endif

	Tr.NTracks++;
	Tr.NTracks1++;
	Tr.TrEvent.push_back(Tr.track_place_holder);

      }

      for (Int_t j=0; j<Si.ReadHit->size(); j++){//loop over all silicon
	Si.hit_place_holder = Si.ReadHit->at(j);
	if ( Si.hit_place_holder.TrackType==1 ){//if they are track type 1, then the event has aleady been copied
	  continue;
	}
	if ( Si.hit_place_holder.Energy <= 0 ){//make sure that the Silicon energy was filled
	  continue;
	}
	Tr.si_place_holder.TrackType = 2;
	Tr.si_place_holder.DetID = Si.hit_place_holder.DetID;
	Tr.si_place_holder.SiEnergy = Si.hit_place_holder.Energy;
	Tr.si_place_holder.SiZ = Si.hit_place_holder.ZW;
	Tr.si_place_holder.SiR = Si.hit_place_holder.RW;
	Tr.si_place_holder.SiPhi = Si.hit_place_holder.PhiW;

	Tr.NTracks2++;
	Tr.NTracks++;
	Tr.SiEvent.push_back(Tr.si_place_holder);

      }

      sort( Tr.SiEvent.begin(), Tr.SiEvent.end(),Tr.Si_sort_method );//sort by si energy, high to low

      for (Int_t j=0; j<PC.ReadHit->size(); j++ ){//loop over pc
	PC.pc_place_holder = PC.ReadHit->at(j);
	if ( PC.pc_place_holder.TrackType == 1 ){//skip events that have already been copied
	  continue;
	}
	Tr.pc_place_holder.TrackType = 3;
	Tr.pc_place_holder.WireID = PC.pc_place_holder.WireID;
	Tr.pc_place_holder.PCEnergy = PC.pc_place_holder.Energy;
	Tr.pc_place_holder.PCZ = PC.pc_place_holder.ZW;
	Tr.pc_place_holder.PCR = PC.pc_place_holder.RW;
	Tr.pc_place_holder.PCPhi = PC.pc_place_holder.PhiW;

	Tr.NTracks3++;
	Tr.NTracks++;
	Tr.PCEvent.push_back(Tr.pc_place_holder);

      }
      sort( Tr.PCEvent.begin(), Tr.PCEvent.end(),Tr.PC_sort_method );//sort by pc energy, high to low

      //done copying-----------------------------------------------------------------------

//--------------------------------------------------------------

//-------------------------------------------------------------------


#ifdef FillTree
      MainTree->Fill();
#endif
    }
    //    delete inputFile;
  }


  outputfile->cd();
  RootObjects->Write();

  //maintreefile->Close();
  outputfile->Close();

}

void MyFill(string name,
	    int binsX, double lowX, double highX, double valueX){
  try{
    fhmap.at(name)->Fill(valueX);
  } catch(out_of_range e) {
    TH1F* newHist = new TH1F(name.c_str(),name.c_str(),
			     binsX,lowX,highX);
    newHist->Fill(valueX);
    fhlist->Add(newHist);
    fhmap[name] = newHist;
  }
}

void MyFill(string name,
	    int binsX, double lowX, double highX, double valueX,
	    int binsY, double lowY, double highY, double valueY){

  try{
    //cout << name << endl;
    fhmap.at(name)->Fill(valueX,valueY);
  } catch(std::out_of_range e) {
    TH2F* newHist = new TH2F(name.c_str(),name.c_str(),
			     binsX,lowX,highX,
			     binsY,lowY,highY);
    newHist->Fill(valueX,valueY);
    fhlist->Add(newHist);
    fhmap[name] = newHist;
  }  
}
