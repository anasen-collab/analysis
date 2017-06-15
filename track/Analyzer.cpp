//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Goal: To Analyze the Tracking and Reconstruction of the 7Be+d-> p + alpha + alpha Events  
// & to Analyze the Elastic Scattering of Deuterons to Calibrate ANASEN for 7Be+d Experiments..
// & for the other (d,p),(d,alpha)..etc..ANASEN experiments with Gas volume target
//
// //To create a dictionary:
//  rootcint -f tr_dict.cxx -c tree_structure.h LinkDef.h
//
// Usage: g++ -o Analyzer tr_dict.cxx LookUp.cpp Analyzer.cpp `root-config --cflags --glibs`
//
// ./Analyzer DataListCal.txt 282_3_4Cal5Analyzer20161102.root cut/He4.root //
//
// Uses Lookup tables instead of doing integration multiple times for Energyloss, 
// Final Energy, Initial Energy & Distance calculation.
//
// Author: Nabin Rijal, 2016 September.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#define FillTree
#define FillEdE_cor
#define CheckBasic
//#define DoCut //read in and apply cut file?
//#define DoLoss //look up energy loss?

#define DiffIP 2 //cm
#define ConvAngle 180./TMath::Pi() //when multiplied, Converts to Degree from Radian 

#define EdE
#define Be8

#define PCWireCal
#define PCPlots //for heavy hits

#define MaxSiHits   500
#define MaxADCHits  500
#define MaxTDCHits  500
#define MaxTracks   100

// ANASEN
#define pcr 3.846284509 //3.75+0.096284509; //correction for the centroid Kx applied
#define La 54.42        //Length of ANASEN gas volume.

// target positions
#define gold_pos 27.7495 //Spacer-0 all the way in, based on geometry measurements we did with Lagy at 2/22/2017
//#define gold_pos 22.9495 //Spacer-1 //4.8cm
//#define gold_pos 16.9495 //Spacer-2 //10.8cm
//#define gold_pos 12.4495 //Spacer-4 //15.3cm
//#define gold_pos  7.4495 //Spacer-5 //20.3cm
//#define gold_pos  1.5495 //Spacer-6 //26.2cm
//#define gold_pos -2.8505 //Spacer-7 //30.6cm

///////////////////Nuclear Masses ///////////////////////////////////////////////////
//nuclear masses //MeV
//NIST values
#define M_P 938.2720813
#define M_N 939.5654133
#define M_D 1875.612928
#define M_3He 2808.391586
#define M_alpha 3727.379378

// Nabin masses
#define M_5Li 4667.6163636366931 //correct
#define M_5He 4667.67970996292
#define M_6Li 5601.518452737
#define M_7Li 6533.83277448969
#define M_7Be 6534.1836677282
#define M_8Be 7454.85043438849

//---MARIA nuclear masses //MeV--------------------------------------
#define M_16O  14895.079
#define M_18Ne 16767.09917
#define M_21Na 19553.56884
#define M_24Mg 22335.79043
#define M_27Al 25126.49834

//Jon masses
#define M_17F  15832.754 
#define M_Ne20 18617.733

//reaction
#define BeamE 61.54 //MeV energy of 17F beam inside Kapton Window.
#define QValue 4.1296 //MeV 17F(a,p)20Ne

/////////////////////////////////////////////////////////////////////////////////////
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
#include <TRandom1.h>
#include <TList.h>
#include <TCutG.h>
#include <TLorentzVector.h>
#include <TVector3.h>

#include "tree_structure.h"
#include "LookUp.h"

using namespace std;
////////////////////////////////////////////////////////////////////////////////////
Int_t FindMaxPC(Double_t phi, PCHit& PC);

void MyFill(string name,int binsX, double lowX, double highX, double valueX);

void MyFill(string name,int binsX, double lowX, double highX, double valueX,
	    int binsY, double lowY, double highY, double valueY);

TList* fhlist;
std::map<string,TH1*> fhmap;
////////////////////////////////////////////////////////////////////////////////////
bool Track::Tr_Sisort_method(struct TrackEvent a,struct TrackEvent b){
  if(a.SiEnergy > b.SiEnergy)
    return 1;
  return 0;
};
////////////////////////////////////////////////////////////////////////////////////
bool Track::Tr_PCsort_method(struct TrackEvent c,struct TrackEvent d){
  if(c.PCEnergy > d.PCEnergy)
    return 1;
  return 0;
};
////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[]) { 
  //Don't know what this does, but libraries won't load without it
  TApplication *myapp=new TApplication("myapp",0,0); 
  
  Int_t numarg=3;
#ifdef DoCut
  numarg=4;
  printf("expecting cut file...\n");
#endif
  if (argc!=numarg) {
    cout << "Error: Wrong Number of Arguments\n";
    exit(EXIT_FAILURE);
  }

  char* file_raw  = new char [300]; // for input .root file
  char* file_cal = new char [300]; // for output .root file

  strcpy( file_raw, argv[1] );
  strcpy( file_cal, argv[2] );

  cout << "Command is " << argv[0] << endl;
  cout << "Input file is " << argv[1] << endl;
  cout << "Output file is "<< argv[2] << endl;
  
#ifdef DoCut
  //////////////// CUTS ////////////////////

  char* file_cut1 = new char[300]; //for allcut
  strcpy( file_cut1, argv[3] );
 
  //He4 cut
  TFile *cut_file1 = new TFile(file_cut1);
  if (!cut_file1->IsOpen()){
    cout << "Cut file1: " << file_cut1 << " could not be opened.\n";
    exit(EXIT_FAILURE);
  }
  TCutG *cut1 = NULL;
  cut1 = (TCutG*)cut_file1->Get("He4");
  
  if (cut1 == NULL){
    cout << "Cut1 does not exist\n";
    exit(EXIT_FAILURE);
  }
  cout << argv[3] << endl;
#endif
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  TObjArray *RootObjects = new TObjArray();
  SiHit Si;
  PCHit PC;
  ///CsIHit CsI;
  Track Tr; 
  Int_t RFTime, MCPTime;
  Int_t Old_RFTime,Old_MCPTime;

  Si.ReadDet = 0;
  Si.ReadHit = 0;
  PC.ReadHit = 0;
  //CsI.ReadHit = 0;

#ifdef DoLoss
  LookUp *E_Loss_7Be = new LookUp("/data0/nabin/Vec/Param/Be7_D2_400Torr_20160614.eloss",M_7Be);
  LookUp *E_Loss_alpha = new LookUp("/data0/nabin/Vec/Param/He4_D2_400Torr_20160614.eloss",M_alpha);
  LookUp *E_Loss_proton = new LookUp("/data0/nabin/Vec/Param/P_D2_400Torr_20160614.eloss",M_P);
  LookUp *E_Loss_deuteron = new LookUp("/data0/nabin/Vec/Param/D2_D2_400Torr_20160614.eloss",M_D); 
  LookUp *E_Loss_3He = new LookUp("/data0/nabin/Vec/Param/He3_D2_400Torr.eloss",M_3He); 
  //LookUp *E_Loss_Li6 = new LookUp("/data0/nabin/Vec/Param/D2_D2_400Torr_20160614.eloss",M_D); 

  E_Loss_7Be->InitializeLookupTables(30.0,200.0,0.01,0.04);
  //E_Loss_7Be->InitializeLookupTables(20.0,60.0,0.01,0.04);//not stopped--65.0

  E_Loss_alpha->InitializeLookupTables(40.0,1300.0,0.01,0.04);
  //E_Loss_alpha->InitializeLookupTables(40.0,2000.0,0.01,0.04);
  //E_Loss_alpha->InitializeLookupTables(30.0,1000.0,0.01,0.1);

  E_Loss_proton->InitializeLookupTables(40.0,16000.0,0.05,0.01);
  //E_Loss_proton->InitializeLookupTables(40.0,20000.0,0.05,0.01);
  //E_Loss_proton->InitializeLookupTables(30.0,10000.0,0.05,0.1);

  E_Loss_deuteron->InitializeLookupTables(30.0,6000.0,0.02,0.04);
  //E_Loss_deuteron->InitializeLookupTables(25.0,4000.0,0.02,0.04);
  //E_Loss_deuteron->InitializeLookupTables(20.0,2000.0,0.02,0.1); 

  E_Loss_3He->InitializeLookupTables(20.0,4000.0,0.02,0.04);
#endif
  ///////////////////////////////////////////////////////////////////////////////////////////////////

  TFile *outputfile = new TFile(file_cal,"RECREATE");
  TTree *MainTree = new TTree("MainTree","MainTree");
  //++
  MainTree->Branch("Tr.NTracks",&Tr.NTracks,"NTracks/I");
  MainTree->Branch("Tr.NTracks1",&Tr.NTracks1,"NTracks1/I");
  MainTree->Branch("Tr.NTracks2",&Tr.NTracks2,"NTracks2/I");
  MainTree->Branch("Tr.NTracks3",&Tr.NTracks3,"NTracks3/I");
  MainTree->Branch("Tr.TrEvent",&Tr.TrEvent);
  MainTree->Branch("RFTime",&RFTime,"RFTime/I");
  MainTree->Branch("MCPTime",&MCPTime,"MCPTime/I");
  
  RootObjects->Add(MainTree);

  fhlist = new TList;
  RootObjects->Add(fhlist);  
  ////////////////////////////////////////////////////  
  Int_t count_2A_1P=0, count_2A=0, count_1A_1P=0;
  Int_t count_1track=0, count_2tracks=0, count_3tracks =0,count_morethan3tracks=0;
  Int_t count_3He=0;
  Int_t count_3He_NT2;
  ////////////////////////////////////////////////////
  //============  Read the root files from the Data list ==============  
  ifstream inFileList;
  inFileList.open(file_raw);

  if (!inFileList.is_open()){
    cout << "In File List not open\n";
    exit(EXIT_FAILURE);
  }
  
  string rootfile;
  char rootfile_char[100];

  while (!inFileList.eof()){//===================loop over all of the incoming root files================
    getline(inFileList,rootfile);
    if (rootfile.empty()){
      //cout << "Root File does not exist: " << rootfile << endl;
      break;
    }
    strcpy(rootfile_char,rootfile.c_str());
    TFile *inputFile = new TFile(rootfile_char);//open root file and make sure it exists-----
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
    //raw_tree->SetBranchAddress("RFTime",&RFTime);
    //raw_tree->SetBranchAddress("MCPTime",&MCPTime);
    raw_tree->SetBranchAddress("RFTime",&Old_RFTime);
    raw_tree->SetBranchAddress("MCPTime",&Old_MCPTime);
 
    Long64_t nentries = raw_tree->GetEntries();
    cout<<"nentries = "<<nentries<<endl;

    Int_t status;
    for (Long64_t i=0; i<nentries; i++){//====================loop over all events=================
      //cout<<" i =  "<<i<<endl;

      status = raw_tree->GetEvent(i);
      //std::cout << "\rDone: " << i*100./nentries << "%          " << std::flush;
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
      ///////////////////////////////////////////////////////////////////////////////////////////////////

      Double_t correct=1.004009623;
      Double_t wrap=272.0975;//546//538
      MCPTime = Old_MCPTime;
      RFTime = Old_RFTime;
      Double_t tof=MCPTime-RFTime;
      Double_t tiemc=MCPTime*correct-RFTime;
      Int_t tbins=600;
      if(MCPTime > 0 && RFTime>0) {
	MyFill("tof",tbins,0,tbins,tof);
	MyFill("tofw",tbins,0,TMath::Nint(wrap),fmod(tof,wrap));
	MyFill("tofwc",tbins,0,TMath::Nint(wrap),fmod((MCPTime*correct-RFTime),wrap));
      }
      
#ifdef MCP_RF_Cut    
      if(MCPTime > 0 && RFTime>0) {
	//cout<<"   RFTime  == "<<RFTime<<"   MCPTime  =="<<MCPTime<<endl;
	
	MyFill("MCP_RF_Wrapped",2*tbins,-tbins,tbins,tof%wrap);

	if( ((tof%wrap)<47) || ((tof%wrap)>118  && (tof%wrap)<320) || (tof%wrap)>384 ) {
	  //if( ((tof% wrap)<60) || ((tof% wrap)>110  && (tof% wrap)<325) || (tof% wrap)>380 ){
	  continue;
	}else{	  
	}
      }else{
	continue;
      }
#endif
      ///////////////////////////////////////////////////////////////////////////////////////////////////
      Tr.zeroTrack();
      Int_t GoodPC = -1;         
      /////////////////////////////////////////////////////////////////////////////////////////////////////    
      //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
      //cout<<"Si.ReadHit->size() = "<<Si.ReadHit->size()<<endl;  
      MyFill("Si_ReadHit_size",500,0,50,Si.ReadHit->size());  

      for (Int_t j=0; j<Si.ReadHit->size(); j++) {//loop over all silicon
	
	Si.hit_obj = Si.ReadHit->at(j);	//if we have a good hit type set the parameters in your new tree

	if ( Si.hit_obj.Energy <= 0 ) {
	  continue;
	}
	else {

	  GoodPC = FindMaxPC(Si.hit_obj.PhiW, PC);

	  if (GoodPC > -1){//if a PC is found do Tracking

	    PC.pc_obj = PC.ReadHit->at(GoodPC);

	    Tr.ZeroTr_obj();

	    Tr.track_obj.TrackType = 1;
	  
	    Tr.track_obj.SiEnergy = Si.hit_obj.Energy;
	    Tr.track_obj.SiPhi = Si.hit_obj.PhiW;	  
	    Tr.track_obj.SiZ = Si.hit_obj.ZW;
	    Tr.track_obj.SiR = Si.hit_obj.RW;
	    Tr.track_obj.DetID = Si.hit_obj.DetID;
	    Tr.track_obj.SiBCh = Si.hit_obj.BackChannel;
	    Tr.track_obj.HitType = Si.hit_obj.HitType;

	    Tr.track_obj.PCEnergy = PC.pc_obj.Energy;
	    Tr.track_obj.PCPhi = PC.pc_obj.PhiW;	  
	    Tr.track_obj.PCZ = PC.pc_obj.ZW;
	    Tr.track_obj.PCR = PC.pc_obj.RW;
	    Tr.track_obj.WireID = PC.pc_obj.WireID;	
	 
	    // eliminate Si and Wire from further tracking
	    Si.ReadHit->at(j).Energy = -1000;
	    PC.ReadHit->at(GoodPC).Energy = -10;

	    Tr.NTracks++;//total no. of tracks...Good or bad
	    Tr.NTracks1++;//good tracks...PC & Si both
	  }else{
	    continue;
	  }
	  Tr.TrEvent.push_back(Tr.track_obj);
	}   
      }     
      sort( Tr.TrEvent.begin(), Tr.TrEvent.begin()+Tr.NTracks1,Tr.Tr_Sisort_method );    
      //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    
      //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      for (Int_t k=0; k<Si.ReadHit->size(); k++){//loop over all silicon
	
	Si.hit_obj = Si.ReadHit->at(k);

	if ( (Si.hit_obj.Energy == -1000) || (Si.hit_obj.Energy <= 0)){ //make sure that the Silicon energy was filled
	  continue;
	}else{
	  Tr.ZeroTr_obj();
	  Tr.track_obj.TrackType = 2;	
	  Tr.track_obj.SiEnergy = Si.hit_obj.Energy;
	  Tr.track_obj.SiPhi = Si.hit_obj.PhiW;
	  Tr.track_obj.SiZ = Si.hit_obj.ZW;	  
	  Tr.track_obj.SiR = Si.hit_obj.RW;	  
	  Tr.track_obj.DetID = Si.hit_obj.DetID;
	  Tr.track_obj.SiBCh = Si.hit_obj.BackChannel;
	  Tr.track_obj.HitType = Si.hit_obj.HitType;

	  Tr.NTracks2++;//Only Si && no PC
	  Tr.NTracks++;//total no. of tracks...Good or bad
	  Tr.TrEvent.push_back(Tr.track_obj);
	}
      }      
      sort( Tr.TrEvent.begin()+Tr.NTracks1, Tr.TrEvent.begin()+Tr.NTracks1+Tr.NTracks2,Tr.Tr_Sisort_method );  
      //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    
      //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      //cout<<"PC.ReadHit->size() = "<<PC.ReadHit->size()<<endl;
      MyFill("PC_ReadHit_size",500,0,50,PC.ReadHit->size());  

      for (Int_t l=0; l<PC.ReadHit->size(); l++ ){//loop over pc
	
	PC.pc_obj = PC.ReadHit->at(l);

	if((PC.ReadHit->at(l).Energy == -10) || (PC.ReadHit->at(l).Energy <= 0)){
	  continue;
	}else{	
	  Tr.ZeroTr_obj();  
	  Tr.track_obj.TrackType = 3;

	  Tr.track_obj.PCEnergy = PC.pc_obj.Energy;
	  Tr.track_obj.WireID = PC.pc_obj.WireID;	
	  Tr.track_obj.PCZ = PC.pc_obj.ZW;
	  Tr.track_obj.PCR = PC.pc_obj.RW;
	  Tr.track_obj.PCPhi = PC.pc_obj.PhiW;

	  Tr.NTracks3++;//Only PC && no Si
	  Tr.NTracks++;//total no. of tracks...Good && bad both
	  Tr.TrEvent.push_back(Tr.track_obj);
	}
      }        
      sort(Tr.TrEvent.begin()+Tr.NTracks1+Tr.NTracks2, Tr.TrEvent.end(),Tr.Tr_PCsort_method );      
      //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      //////////////////////////////////////////////////////////////////////////////////////
      //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      //cout<<"Tr.NTracks = "<<Tr.NTracks<<" Tr.NTracks1 = "<<Tr.NTracks1<<" Tr.NTracks2 = "<<Tr.NTracks2<<" Tr.NTracks3 = "<<Tr.NTracks3<<endl;      
    
      /////////////////////////////////////////////////////////////////////////////////////////////
      ////////////// checking for the heavy hit &/or cross talk in the wire ///////////////////////
      /////////////////////////////////////////////////////////////////////////////////////////////    
#ifdef PCPlots
      Int_t pct;
      for(Int_t pc=0; pc<Tr.NTracks; pc++) {
	//cout<<"   Check 10 " <<Tr.NTracks<<endl;
	
	//All PCWire vs Energy
	MyFill("WireID_vs_PCEnegy",25,0,24,Tr.TrEvent[pc].WireID,500,0,2,Tr.TrEvent[pc].PCEnergy);
	
	for(Int_t pca=0; pca<Tr.NTracks1; pca++) {
	  
	  //PCWire with track in channel 12 & other tracks & non-tracks in other channels
	  MyFill("WireID_mod1_vs_PCEnegy",
		 25,0,24,(((Int_t)Tr.TrEvent[pc].WireID-(Int_t)Tr.TrEvent[pca].WireID +12)%24),
		 500,0,2,Tr.TrEvent[pc].PCEnergy);
	}
      }
#endif
      /////////////////////////////////////////////////////////////////////////////////////////////



      //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ///////////////////////////////////For the Tracking///////////////////////////////////    
      //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++        
      //reconstruction variables
      Double_t m = 0, b = 0; 

      for(Int_t p=0; p<Tr.NTracks1;p++) {
	
	//CCCCCCCCCCCCCCCCCCCCCCCCCC///Let's put some checks //2016July28 CCCCCCCCCCC	
#ifdef CheckBasic
	if(Tr.TrEvent[p].SiZ < -4.0 || Tr.TrEvent[p].PCZ < -4.0){
	  continue;
	}
	if(Tr.TrEvent[p].PCR>4.0 || Tr.TrEvent[p].PCR<3.5){
	  continue;
	}
	if(Tr.TrEvent[p].SiR>11.0 || Tr.TrEvent[p].SiR<4.0){
	  continue;
	}
#endif
	//CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC	  
	//if(Tr.TrEvent[p].SiZ > 0 && Tr.TrEvent[p].PCZ <= 0){
	////cout<<" SiZ ="<<Tr.TrEvent[p].SiZ<<" PCZ ="<<Tr.TrEvent[p].PCZ<<endl;
	//}
	//cout<<"SiR = "<<Tr.TrEvent[p].SiR<<" PCR ="<<Tr.TrEvent[p].PCR<<" SiZ ="<<Tr.TrEvent[p].SiZ<<" PCZ ="<<Tr.TrEvent[p].PCZ<<endl;

	m = (Tr.TrEvent[p].PCR-Tr.TrEvent[p].SiR)/(Tr.TrEvent[p].PCZ-Tr.TrEvent[p].SiZ);
	b = (Tr.TrEvent[p].PCR - m*Tr.TrEvent[p].PCZ);
	Tr.TrEvent[p].IntPoint = -b/m;
	////cout<<" IntPoint = "<<Tr.TrEvent[p].IntPoint<<endl;
	//CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#ifdef CheckBasic
	if(Tr.TrEvent[p].IntPoint<-5.0 || Tr.TrEvent[p].IntPoint>55.0){
	  continue;
	}
#endif
	//CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	/////////////////////Calculating Theta and PathLength/////////////////////////////////////////////
	// Theta is the angle between particle and beam, but our beam points in the negative z-direction 
	if ((Tr.TrEvent[p].IntPoint - Tr.TrEvent[p].SiZ) > 0){
	  // This is forward-scattering, so we have theta b/w 0 and 90 
	  Tr.TrEvent[p].Theta = atan(Tr.TrEvent[p].SiR/(Tr.TrEvent[p].IntPoint - Tr.TrEvent[p].SiZ));
	  //Tr.TrEvent[p].Theta = TMath::Pi() - atan(Tr.TrEvent[p].SiR/(Tr.TrEvent[p].IntPoint - Tr.TrEvent[p].SiZ));
	  Tr.TrEvent[p].PathLength = Tr.TrEvent[p].SiR/sin(Tr.TrEvent[p].Theta);
	 
	  //cout<<" Tr.TrEvent[p].Theta1 =  "<<Tr.TrEvent[p].Theta*ConvAngle<<" Tr.TrEvent[p].PathLength1 = "<<Tr.TrEvent[p].PathLength<<endl;
	}
	else if ((Tr.TrEvent[p].IntPoint - Tr.TrEvent[p].SiZ) < 0){
	  Tr.TrEvent[p].Theta = TMath::Pi() + atan(Tr.TrEvent[p].SiR/(Tr.TrEvent[p].IntPoint - Tr.TrEvent[p].SiZ));
	  //Tr.TrEvent[p].Theta = atan(Tr.TrEvent[p].SiR/(Tr.TrEvent[p].IntPoint - Tr.TrEvent[p].SiZ));
	  Tr.TrEvent[p].PathLength = Tr.TrEvent[p].SiR/sin(Tr.TrEvent[p].Theta);
	
	  //if(Tr.TrEvent[p].Theta>0){
	  //cout<<" Tr.TrEvent[p].Theta2 =  "<<Tr.TrEvent[p].Theta*ConvAngle<<" Tr.TrEvent[p].PathLength2 = "<<Tr.TrEvent[p].PathLength<<endl;
	  //}
	}
	else{
	  Tr.TrEvent[p].Theta = TMath::Pi()/2;
	  Tr.TrEvent[p].PathLength = Tr.TrEvent[p].SiR;

	  //cout<<" Tr.TrEvent[p].Theta3 =  "<<Tr.TrEvent[p].Theta*ConvAngle<<" Tr.TrEvent[p].PathLength3 = "<<Tr.TrEvent[p].PathLength<<endl;
	}
	//////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef DoLoss
	if(Tr.TrEvent[p].IntPoint >0.0 && Tr.TrEvent[p].IntPoint<54.0) {
	  Tr.TrEvent[p].EnergyLoss = E_Loss_7Be->GetEnergyLoss(BeamE,(La-Tr.TrEvent[p].IntPoint));
	  //Tr.TrEvent[p].BeamEnergy = BeamE - Tr.TrEvent[p].EnergyLoss;
	  
	  if((La-Tr.TrEvent[p].IntPoint)>0.0 && (La-Tr.TrEvent[p].IntPoint)<54.0) {
	    Tr.TrEvent[p].BeamEnergy = E_Loss_7Be->GetLookupEnergy(BeamE,(La-Tr.TrEvent[p].IntPoint));
	  }	
	}
#endif
	////////////////////////////////////////////////////////////////////////////////////////
      }//end of for loop Tracking.
      ////////////////////////////////////////////////////////////////////////////////////////
      
#ifdef PCWireCal      
      Double_t tantheta;
           
      for(Int_t s=0; s<Tr.NTracks1;s++) {

	//determine PC position from Silicon position and gold position
	tantheta = Tr.TrEvent[s].SiR/(Tr.TrEvent[s].SiZ - gold_pos);
	Tr.TrEvent[s].PCZ_Ref = pcr/tantheta+gold_pos;
	//cout<<"tantheta = "<< tantheta <<" PCZ_Ref = "<<Tr.TrEvent[s].PCZ_Ref<<endl;

	Int_t pcbins=400;
	
	MyFill(Form("PCZ_Ref%i",Tr.TrEvent[s].WireID),pcbins,1,30,Tr.TrEvent[s].PCZ_Ref);
	// before PCWIRECAL applied
	MyFill(Form("PCZ_vs_Z%i",Tr.TrEvent[s].WireID),
	       pcbins,-1.5,1.5,Tr.TrEvent[s].PCZ,
	       pcbins,1.0,30.0,Tr.TrEvent[s].PCZ_Ref); 

	// after PCWIRECAL applied
	MyFill(Form("PCZ_vs_Zc%i",Tr.TrEvent[s].WireID),
	       pcbins,-1.0,30.0,Tr.TrEvent[s].PCZ,
	       pcbins,-1.0,30.0,Tr.TrEvent[s].PCZ_Ref); 

	MyFill("PCZ_vs_Zc",
	       pcbins,-1.0,30.0,Tr.TrEvent[s].PCZ,
	       pcbins,-1.0,30.0,Tr.TrEvent[s].PCZ_Ref);
	
	if(Tr.TrEvent[s].DetID<16 && Tr.TrEvent[s].DetID>-1) {
	  MyFill(Form("PCZ_vs_Zc_q3_r1%i",Tr.TrEvent[s].WireID),
		 pcbins,-1.0,30.0,Tr.TrEvent[s].PCZ,
		 pcbins,-1.0,30.0,Tr.TrEvent[s].PCZ_Ref);
	  MyFill("PCZ_vs_Zc_q3_r1",
		 pcbins,-1.0,30.0,Tr.TrEvent[s].PCZ,
		 pcbins,-1.0,30.0,Tr.TrEvent[s].PCZ_Ref);
	}

	if(Tr.TrEvent[s].DetID<28 && Tr.TrEvent[s].DetID>3) {
	  MyFill(Form("PCZ_vs_Zc_r1_r2%i",Tr.TrEvent[s].WireID),pcbins,-1.0,30.0,Tr.TrEvent[s].PCZ,pcbins,-1.0,30.0,Tr.TrEvent[s].PCZ_Ref);
	}
	
	if (Tr.TrEvent[s].SiEnergy>9.4 && Tr.TrEvent[s].SiEnergy< 9.7) {//this cut is to clean up calibration data... 	  
	  
	  MyFill(Form("PCZ_Refg%i",Tr.TrEvent[s].WireID),
		 pcbins,1.0,30.0,Tr.TrEvent[s].PCZ_Ref);
	  MyFill(Form("PCZ_vs_Zg%i",Tr.TrEvent[s].WireID),
		 pcbins,-1.5,1.5,Tr.TrEvent[s].PCZ,
		 pcbins,1.0,30.0,Tr.TrEvent[s].PCZ_Ref);
	  MyFill(Form("PCZ_vs_Zgc%i",Tr.TrEvent[s].WireID),
		 pcbins,-1.0,30.0,Tr.TrEvent[s].PCZ,
		 pcbins,-1.0,30.0,Tr.TrEvent[s].PCZ_Ref); // after PCWIRECAL applied
	  /////////////////////////////////////////////////////////////////	 
	}
	//MyFill(Form("PCZ_vs_Z%i",Tr.TrEvent[s].WireID),100,-2.0,2.0,Tr.TrEvent[s].PCZ,750,0,30,Tr.TrEvent[s].PCZ_Ref);
      }
#endif

      /////////////////////////////////////////////////////////////////////////////////////
#ifdef FillEdE_cor
      Float_t demin=-0.01;
      Float_t demax=0.25;
      Int_t debins=600;
      for(Int_t q=0; q<Tr.NTracks1;q++) {
	MyFill("E_de",debins,-1,29,Tr.TrEvent[q].SiEnergy,debins,demin,demax,Tr.TrEvent[q].PCEnergy);
	//MyFill("E_de_corrected",debins,-1,35,Tr.TrEvent[q].SiEnergy,debins,demin,demax,Tr.TrEvent[q].PCEnergy*Tr.TrEvent[q].PathLength);
	MyFill("E_de_corrected",debins,-1,29,Tr.TrEvent[q].SiEnergy,debins,demin,demax,Tr.TrEvent[q].PCEnergy *sin(Tr.TrEvent[q].Theta));	  
	///MyFill(Form("PCZ%i",PC.WireID[GoodPC]),300,-10,50,Tr.TrEvent[q].PCZ);	    
	MyFill("InteractionPoint",300,-10,50,Tr.TrEvent[q].IntPoint);	
      }
#endif

      ///////////////////////// DO your analysis here.... /////////////////////////////////////// 
   


      //////////////////////////////////////////////////////////////////////////////////////////// 
#ifdef FillTree
      MainTree->Fill();
#endif
      //////////////////////////////////////////////////////////////////////////////
    }
    ////////////////////////////////////////////////////////////////////////////////   
  } 
  //////////////////////////////////////////////////////////////////////////////////
  outputfile->cd();
  RootObjects->Write(); 
  outputfile->Close();
}//end of Main
////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////
Float_t phidiff ( Float_t phi1,Float_t phi2)
{
  return (fmodf((fabs(phi1-phi2) + 2*TMath::Pi()), 2*TMath::Pi()));
}
/////////////////////////////////////////////////////////////////////////////////////
/*
// Finds a maximum PC within a given phi range
// Nabin Rijal, June 2016

  Int_t FindMaxPC(Double_t phi, PCHit& PC){
  Int_t GoodPC = -1;
  Double_t MaxPC = -10;
  //Double_t MinPhi = 0.2619;
  Double_t MinPhi = 0.5238;

  for (int k=0; k<PC.NPCHits; k++){//loop over the pc hits
    //if the PC falls in a range of phi then it is possible correlated
    //we find the maximum energy on the pc
    PC.pc_obj = PC.ReadHit->at(k);

    //if (PC.pc_obj.TrackType == 1){
    //continue;
    //}

    if ( (fabs(PC.pc_obj.PhiW-phi) <= MinPhi) || ((2*TMath::Pi() - fabs(PC.pc_obj.PhiW-phi)) <= MinPhi) ) {
      if ( PC.pc_obj.Energy >= MaxPC ){
	MaxPC = PC.pc_obj.Energy;
	GoodPC = k;
      }
    }
  }
  return GoodPC;
}
*/

///////////////////////////////////////////////////////////////////////////////////
// If there are more than one silicon firing within the range of given phi, 
// picks the one with closer phi and assigns the another pc hit to the next silicon.

// Nabin Rijal, September 18, 2016

Int_t FindMaxPC(Double_t phi, PCHit& PC){
 
  Int_t MaxPCindex = -1,NexttoMaxPCindex =-1;
  Double_t MaxPC = -10, NexttoMaxPC = -10;

  //Double_t MinPhi = 0.2619;
  Double_t MinPhi = 0.5238;

  for (int k=0; k<PC.NPCHits; k++){//loop over the pc hits 
    //if the PC falls in a range of phi then it is possible correlated
    //we find the maximum energy on the pc
  
    PC.pc_obj = PC.ReadHit->at(k);   

    if ( phidiff(PC.pc_obj.PhiW,phi) < MinPhi) {      

      if ( PC.pc_obj.Energy >= MaxPC ){
	NexttoMaxPC = MaxPC;
	NexttoMaxPCindex = MaxPCindex;
	MaxPC = PC.pc_obj.Energy;
	MaxPCindex = k;

      }
      else if (PC.pc_obj.Energy >= NexttoMaxPC ){
	NexttoMaxPC = PC.pc_obj.Energy;
	NexttoMaxPCindex = k;
      }else{
	continue;
      }
    }    
  }  
  if (NexttoMaxPCindex>0){
    // there is an ambiguity: pick the one with smaller DeltaPhi
    if (phidiff( PC.ReadHit->at(NexttoMaxPCindex).PhiW,phi) > phidiff(PC.ReadHit->at(NexttoMaxPCindex).PhiW,phi)){
      return(NexttoMaxPCindex);
    }
    else
      return(MaxPCindex);      
  }  
  return MaxPCindex; 
}
/////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////
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
////=================================================================================
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
/////////////////////////////////////////////////////////////////////////////////////
