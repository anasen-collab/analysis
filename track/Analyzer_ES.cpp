//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Goal: To Analyze the Elastic Scattering of Deuterons to Calibrate ANASEN for 7Be+d Experiments..
// & for the other (d,p),(d,alpha)..etc..ANASEN experiments with Gas volume target
//
// //To create a dictionary:
//  rootcint -f tr_dict.cxx -c ../Include/tree_structure.h LinkDef.h
//
// Usage: g++ -o Analyzer_ES tr_dict.cxx LookUp.cpp Analyzer_ES.cpp `root-config --cflags --glibs`
//
// ./Analyzer_ES DataListCal.txt 2430Cal5Analyzer20170303.root cut/D2.root //
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

#define DiffIP 2 //cm
#define ConvAngle 57.27272727 //when multiplied, Converts to Degree from Radian 

#define EdE
#define Be8

#define PCWireCal

#define MaxSiHits   500
#define MaxADCHits  500
#define MaxTDCHits  500
#define MaxTracks   100

#define BeamE 19.6 //Energy of 7Be beam inside Kapton Window.
#define pcr 3.846284509 //3.75+0.096284509; //correction for the centroid Kx applied
#define La 55.0545   //Length of ANASEN gas volume as measured 2/22/2017 with Lagy

///////////////////Nuclear Masses ///////////////////////////////////////////////////
//nuclear masses //MeV
//NIST values
#define M_P 938.2720813
#define M_N 939.5654133
#define M_D2 1875.612928
#define M_3He 2808.391586
#define M_alpha 3727.379378

#define M_Be8 7454.85043438849
#define M_Li5 4667.6163636366931 //correct
//#define M_Li5 4665.7163636366931 //correction of -1.90 MeV applied

#define M_Li6 5601.518452737

#define M_Be7 6534.1836677282 

#define M_Li7 6533.83277448969
#define M_He5 4667.67970996292

#define QValue 
#define gold_pos 28.9
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

#include "../include/tree_structure.h"
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

  if (argc!=4) {
    cout << "Error: Wrong Number of Arguments\n";
    exit(EXIT_FAILURE);
  }

  char* file_raw  = new char [100]; // for input .root file
  char* file_cal = new char [100]; // for output .root file

  strcpy( file_raw, argv[1] );
  strcpy( file_cal, argv[2] );
  //
  //////////////// CUTS ////////////////////
  //
  char* file_cut1 = new char[100]; 
  strcpy( file_cut1, argv[3] );
 

  //D2 cut
  TFile *cut_file1 = new TFile(file_cut1);
  if (!cut_file1->IsOpen()){
    cout << "Cut file1: " << file_cut1 << " could not be opened.\n";
    exit(EXIT_FAILURE);
  }
  TCutG *cut1 = NULL;
  cut1 = (TCutG*)cut_file1->Get("D2");
  
  if (cut1 == NULL){
    cout << "Cut1 does not exist\n";
    exit(EXIT_FAILURE);
  } 
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  TObjArray *RootObjects = new TObjArray();
  SiHit Si;
  PCHit PC;

  Track Tr; 
  Int_t RFTime, MCPTime; 

  Si.ReadDet = 0;
  Si.ReadHit = 0;
  PC.ReadHit = 0;

  LookUp *E_Loss_7Be = new LookUp("/data0/nabin/Vec/Param/Be7_D2_400Torr_20160614.eloss",M_Be7);
  LookUp *E_Loss_deuteron = new LookUp("/data0/nabin/Vec/Param/D2_D2_400Torr_20160614.eloss",M_D2); 
  E_Loss_7Be->InitializeLookupTables(30.0,200.0,0.01,0.04);
  E_Loss_deuteron->InitializeLookupTables(30.0,6000.0,0.02,0.04); 
  ///////////////////////////////////////////////////////////////////////////////////////////////////

  TFile *outputfile = new TFile(file_cal,"RECREATE");
  TTree *MainTree = new TTree("MainTree","MainTree");
  //++
  MainTree->Branch("Tr.NTracks",&Tr.NTracks,"NTracks/I");
  MainTree->Branch("Tr.NTracks1",&Tr.NTracks1,"NTracks1/I");
  MainTree->Branch("Tr.NTracks2",&Tr.NTracks2,"NTracks2/I");
  MainTree->Branch("Tr.NTracks3",&Tr.NTracks3,"NTracks3/I");
  MainTree->Branch("Tr.TrEvent",&Tr.TrEvent);   
  
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
    raw_tree->SetBranchAddress("RFTime",&RFTime);
    raw_tree->SetBranchAddress("MCPTime",&MCPTime);
 
    Long64_t nentries = raw_tree->GetEntries();
    cout<<"nentries = "<<nentries<<endl;

    Int_t status;
    for (Long64_t i=0; i<nentries; i++){//====================loop over all events=================
      //cout<<" i =  "<<i<<endl;

      status = raw_tree->GetEvent(i);

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
#ifdef MCP_RF_Cut    
      if(MCPTime > 0 && RFTime>0){
	//cout<<"   RFTime  == "<<RFTime<<"   MCPTime  =="<<MCPTime<<endl;
	MyFill("MCP_RF_Wrapped",400,-600,600,(MCPTime-RFTime)%538);

	if( (((MCPTime - RFTime)% 538)<47) || (((MCPTime - RFTime)% 538)>118  && ((MCPTime - RFTime)% 538)<320) || ((MCPTime - RFTime)% 538)>384 ){
	  //if( (((MCPTime - RFTime)% 538)<60) || (((MCPTime - RFTime)% 538)>110  && ((MCPTime - RFTime)% 538)<325) || ((MCPTime - RFTime)% 538)>380 ){
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
    
      //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ///////////////////////////////////For the Tracking///////////////////////////////////    
      //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++        
      //reconstruction variables
      Double_t m = 0, b = 0; 

      for(Int_t p=0; p<Tr.NTracks1;p++){
	
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
	if(Tr.TrEvent[p].IntPoint >0.0 && Tr.TrEvent[p].IntPoint<54.0){
	  Tr.TrEvent[p].EnergyLoss = E_Loss_7Be->GetEnergyLoss(BeamE,(La-Tr.TrEvent[p].IntPoint));
	  //Tr.TrEvent[p].BeamEnergy = BeamE - Tr.TrEvent[p].EnergyLoss;
	    
	  if((La-Tr.TrEvent[p].IntPoint)>0.0 && (La-Tr.TrEvent[p].IntPoint)<54.0){
	    Tr.TrEvent[p].BeamEnergy = E_Loss_7Be->GetLookupEnergy(BeamE,(La-Tr.TrEvent[p].IntPoint));
	  }	
	}
	////////////////////////////////////////////////////////////////////////////////////////
      }//end of for loop Tracking.
      ////////////////////////////////////////////////////////////////////////////////////////
      ///////////////////////////////////////////////////////////////////////////////////////////
      ///////////////////////////// Elastic scattering of Deuterons /////////////////////////////
      ///////////////////////////////////////////////////////////////////////////////////////////    
      Float_t E_deut_rxn =0.0;
      Float_t Energy_7Be_D2=0.0;
      Float_t Theta_7Be_D2=0.0;

      for(Int_t c=0; c<Tr.NTracks1;c++){
	 
	if (cut1->IsInside(Tr.TrEvent[c].SiEnergy,Tr.TrEvent[c].PCEnergy*sin(Tr.TrEvent[c].Theta))){

	  MyFill("D2_E_si",500,0,20,Tr.TrEvent[c].SiEnergy);
	  MyFill("D2_E_si_vs_Theta",500,0,200,Tr.TrEvent[c].Theta*ConvAngle,500,0,15,Tr.TrEvent[c].SiEnergy);
	  MyFill("D2_E_si_vs_IntPoint",500,0,100,Tr.TrEvent[c].IntPoint,500,0,15,Tr.TrEvent[c].SiEnergy);

	  if(Tr.TrEvent[c].DetID<4 && Tr.TrEvent[c].DetID>-1){
	    MyFill("D2_E_si_Q3",500,0,20,Tr.TrEvent[c].SiEnergy);
	    MyFill("D2_E_si_vs_Theta_Q3",500,0,200,Tr.TrEvent[c].Theta*ConvAngle,500,0,15,Tr.TrEvent[c].SiEnergy);
	    MyFill("D2_E_si_vs_IntPoint_Q3",500,0,100,Tr.TrEvent[c].IntPoint,500,0,15,Tr.TrEvent[c].SiEnergy);
	  }else if(Tr.TrEvent[c].DetID<16 && Tr.TrEvent[c].DetID>3){
	    MyFill("D2_E_si_SX3_1",500,0,20,Tr.TrEvent[c].SiEnergy);
	    MyFill("D2_E_si_vs_Theta_SX3_1",500,0,200,Tr.TrEvent[c].Theta*ConvAngle,500,0,15,Tr.TrEvent[c].SiEnergy);
	    MyFill("D2_E_si_vs_IntPoint_SX3_1",500,0,100,Tr.TrEvent[c].IntPoint,500,0,15,Tr.TrEvent[c].SiEnergy);
	  }else if(Tr.TrEvent[c].DetID<28 && Tr.TrEvent[c].DetID>15){
	    MyFill("D2_E_si_SX3_2",500,0,20,Tr.TrEvent[c].SiEnergy);
	    MyFill("D2_E_si_vs_Theta_SX3_2",500,0,200,Tr.TrEvent[c].Theta*ConvAngle,500,0,15,Tr.TrEvent[c].SiEnergy);
	    MyFill("D2_E_si_vs_IntPoint_SX3_2",500,0,100,Tr.TrEvent[c].IntPoint,500,0,15,Tr.TrEvent[c].SiEnergy);
	  }

	  if(Tr.TrEvent[c].SiEnergy>0.0 && Tr.TrEvent[c].SiEnergy<20.0 && Tr.TrEvent[c].PathLength>0.0 && Tr.TrEvent[c].PathLength< 100.0){
	    E_deut_rxn = E_Loss_deuteron->GetLookupEnergy(Tr.TrEvent[c].SiEnergy,(-Tr.TrEvent[c].PathLength));
	    MyFill("D2_E_rxn",500,0,20,E_deut_rxn);


	    //Be-7 Energy Calculation from the Elastic scattering:
	    Energy_7Be_D2 = ( (M_Be7+M_D2)*(M_Be7+M_D2)*E_deut_rxn/(4*M_Be7*M_D2*cos(Tr.TrEvent[c].Theta)*cos(Tr.TrEvent[c].Theta)));
	    //Be-7 Theta Calculation from the Elastic scattering:
	    Theta_7Be_D2 = asin((sin(Tr.TrEvent[c].Theta)*sqrt(M_D2/M_Be7))/sqrt((Energy_7Be_D2/E_deut_rxn)-1));

	    MyFill("D2_7Be_Energy",1000,0,25,Energy_7Be_D2);
	    MyFill("D2_7Be_Energy_VS_BeamEnergy",1000,0,25,Energy_7Be_D2,1000,0,25,Tr.TrEvent[c].BeamEnergy);

  
	    if(Tr.TrEvent[c].DetID<4 && Tr.TrEvent[c].DetID>-1){
	      MyFill("D2_7Be_Energy_Q3",1000,0,25,Energy_7Be_D2);
	      MyFill("D2_7Be_Energy_VS_BeamEnergy_Q3",1000,0,25,Energy_7Be_D2,1000,0,25,Tr.TrEvent[c].BeamEnergy);	     
	    }else if(Tr.TrEvent[c].DetID<16 && Tr.TrEvent[c].DetID>3){
	      MyFill("D2_7Be_Energy_SX3_1",1000,0,25,Energy_7Be_D2);
	      MyFill("D2_7Be_Energy_VS_BeamEnergy_SX3_1",1000,0,25,Energy_7Be_D2,1000,0,25,Tr.TrEvent[c].BeamEnergy);	     
	    }else if(Tr.TrEvent[c].DetID<28 && Tr.TrEvent[c].DetID>15){
	      MyFill("D2_7Be_Energy_SX3_2",1000,0,25,Energy_7Be_D2);
	      MyFill("D2_7Be_Energy_VS_BeamEnergy_SX3_2",1000,0,25,Energy_7Be_D2,1000,0,25,Tr.TrEvent[c].BeamEnergy);	     
	    }
	  }
	}
      }
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
