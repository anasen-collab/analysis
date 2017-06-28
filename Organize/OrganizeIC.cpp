//Written by John Parker--28 Apr, 2016
//in progress
//Usage-slightly more complicated
//This file creates a tree with more than one dynamic dimension
//To do so, we use vectors and have to create a dictionary
//The tree structure is defined in ../Include/organizetree.h
//To create a dictionary:
//  rootcint -f organize_dictionary.cxx -c ../Include/organizetree.h LinkDef.h
//To compile the code:
//  g++ -o OrganizeIC organize_dictionary.cxx OrganizeIC.cpp `root-config --cflags --glibs`
//To run the code:
//  ./OrganizeIC ../../../../../../data0/manasta/evt2root_files_new/runXXX.root ../../../../../../data0/manasta/OrganizeRaw_files/runXXX.root
//
//Obviously need to change headers to point to where you have them and all of the arguments of the ChannelMap functions
//
//
//-------------------------------------------------------------------------------
//Tree structure
//The raw root file is organized by channel number
//This code first organizes the data by detector number
//Si.NSiHits is the number of silicon detectors that fired
//Si.Detector is a branch that has the detector information broken down by channel
//     i.e. Si.Detector.BackMult is the back multiplicity
//     Note: The QQQ front channel information is placed in Si.Detector.Up*
//           The Down information should be empty
//     DetID, UpMult, DownMult, BackMult and HitType are 1D vectors of size NSiHits
//     The others are 2D vectors of size [NSiHits][UpMult] (or DownMult or BackMult)
//     Note: The HitType is just a shorthand convention for express back and front multiplicity
//           For QQQ, hit type is BackMult:FrontMult (i.e. back mult=2, front mult = 1 is hit type 21)
//           For SX3, hit type is Back:Up:Down (i.e. back mult=2, up mult=1, down mult=0 is 210)
//     The energy is saved so that each step of the energy calibration can be seen and is named in a self explanatory manner, I think.
//-----------------------------------------
//Si.Hit is a branch that condenses the information in Si.Detector into 1 hit
//This is very much under development...
//Each detector hit is treated as one hit, regardless of the channel multiplicities
//    All of the leafs are size NSiHits
//    NHitsInDet is always 1 (could be more if we want to)
//    The QQQ and SX3 are handled differently.
//    Right now, the SX3 is handled case by case based off of hit type
//    The QQQ is treated in a general manner, which I will describe here (briefly) (SortSilicon.h)
//      First, All adjacent channels are compressed into 'clusters' (the energies are added back)
//         An artificial channel number is assigned to each cluster by weighting the channels by the energy
//          i.e. cluster channel = (Ch1*E1 + Ch2*E2)/(E1+E2); (front and back are done separately)
//         A new hit type is assigned based on the compressed multiplicities
//          i.e. If the old hit type was 21 but the back channels were adjacent, we now assign a hit type=11
//      Second, the total front energy is compared with the back energy. If they agree within 10% we accept this event.
//            The energy of this event is given as the total back energy
//            The positions are determined by channel number (front channel->R, back channel->Phi)
//            An artificial channel number is assigned in the same manner as above, weighting for energies
//             i.e. ChNum = (ClusterCh1*E1 + ClusterCh2*E2 + ...)/TotalE
//        If they do not agree within 10%, we check to see if any of the individual clusters match
//          i.e. For a new hit type 21, noise in one back channel could cause event the back energies to not match the front energies. But one of the back channels might match very well. If so, we match 1 to 1 and ignore the noise. (flag = 1)
//          Once we find a match, we stop looking
//        If we still do not find a match, we look to see if removing one front channel allows the front energy to match the back energy. If so, we ignore the problem front channel and treat the rest as we do above.
//            And same for the back (flag = 2)
//        If we still do not find a match, we give up and treat it as if they match
//            We could cut these events out altogether, but the front and back energies are saved, so we can do the cut whenever
//      Third, we fill the tree...
//         Flag exists so that I can flag events that meet certain conditions (that I've flagged above...)
//         TrackType is initially set as 2 (I will explain below)
//         Everything else should be self explanatory-ish
//------------------------------------------------
//Now, we handle the PC
//We require the up and down to fire within a range, or we throw it out
//PC.Hit.TrackType = 3
//------------------------------------------------
//Now we Track...
//This can be turned on or off by '#define Tracking'
//For each Si hit, it looks for a PC that fired close by (as defined by phi)
//   If more than one pc fired, it takes the largest
//   If a pc is found, the Si.Hit.TrackType is changed to 1 for this silicon
//      and the PC.Hit.TrackType is changed to 1 for the corresponding PC
//      The energy and position for the Si and PC are copied into the Tr.Event
//      Based off of the Z and R position, an interaction point, path length, and theta are found.
//   If a pc is NOT found, then the Si.Hit.TrackType remains 2
//   The TrackType exists so that I know that a given Silicon has already been associated with a PC
//      I don't want to use the same Silicon or PC twice...
//   If you are doing a proton or alpha scattering off of gold, then #define PC_Pos_Cal and set the gold position in the relevant part of the code. It will create a branch PCZ_Ref that corresponds to the expected PC position based off of the Silicon and gold position.
//---------------------------------------------------------------------------------------------------------
//Useful ROOT stuff
//I haven't defined many histos. Feel free to do so
//I have two MyFill functions defined (for 1d and 2d histos)--very convenient
//Or you can define histos in a standard way
//I draw things from a root terminal
//E-DE:
//  MainTree->Draw("Tr.TrEvent.PCEnergy*TMath::Sin(Tr.TrEvent.Theta):Tr.TrEvent.SiEnergy>>hist(300,0,30,300,0,0.8)","","colz")
//    This draws the entire E-DE for all wires and Si
// If you want to draw it for one wire (or for any condition), you can put a condition as the second argument
//   MainTree->Draw("Tr.TrEvent.PCEnergy*TMath::Sin(Tr.TrEvent.Theta):Tr.TrEvent.SiEnergy>>hist(300,0,30,300,0,0.8)","Tr.TrEvent.WireID==3","colz")
//Proton or alpha cal
//  MainTree->Draw("Tr.Event.PCZ:Tr.Event.PCZ_Ref>>hist(bins...)","Tr.Event.WireID==1","colz")


////////////////////////////////////////////////////////////////////////////////////////////
//// Removed the randomization from the SiEnergy_Pulser calculation
////
//// Add a new histogram on the SX3s "down_vs_up_divideBack%i_front%i" in case it is need for the relative calibration of the detectors.
//// If necessary create also for the QQQs. 
////
//// Created an output file to dump the excluded Entry# and the SiHit to check the consistency of the code due to the randomizations
////
//// For the moment filling the histograms using Energy(Raw) > 0 not EnergyCal > 0
////
//// Edited by Maria Anastasiou, 2016Sept21
////




#define M_PI  3.14159265358979323846264338328 // Pi 
#define Tracking
#define PC_Pos_Cal
//#define FillHists
//#define TimingCut

#define MaxADCHits  500
#define MaxTDCHits  500
#define MaxTracks   100

#include <TH2.h>
#include <TStyle.h>
#include <TMath.h>
#include <TRandom3.h>
#include <stdexcept>
#include <map>
#include <fstream>
#include <string>
//#include <vector>

#include <TObjArray.h>
#include <TList.h>
#include <TFile.h>
//#include <TROOT.h>
#include <TTree.h>
#include <TApplication.h>

#define MaxPCHits 24
#define NPCWires  24

using namespace std;

void MyFill(string name,
	    int binsX, double lowX, double highX, double valueX,
	    int binsY, double lowY, double highY, double valueY);
void MyFill(string name,
	    int binsX, double lowX, double highX, double valueX);

#include "/home/manasta/Desktop/parker_codes/Include/ChannelMap.h"
#include "/home/manasta/Desktop/parker_codes/Include/2016_detclass.h"
#include "/home/manasta/Desktop/parker_codes/Include/SortSilicon.h"

Int_t FindMaxPC(Double_t phi, PCHit& PC);

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

  SiHit Si;
  PCHit PC;
  Track Tr;
  Int_t RFTime,MCPTime;
  Double_t IC;
  Double_t E_IC;
  SortSilicon SiSort;

  //read in command line arguments
  char* filename_histout = new char [100];//output root file
  char* filename_callist = new char [100];//input file list
  
  strcpy( filename_callist, argv[1] );
  strcpy( filename_histout, argv[2] );

  TFile *outputFile = new TFile(filename_histout,"RECREATE");
  //define tree and branches
  TTree *MainTree = new TTree("MainTree","MainTree");
  
  MainTree->Branch("Si.NSiHits",&Si.NSiHits,"NSiHits/I");
  MainTree->Branch("Si.Detector",&Si.Detector);
  MainTree->Branch("Si.Hit",&Si.Hit);
  MainTree->Branch("PC.NPCHits",&PC.NPCHits,"NPCHits/I");
  MainTree->Branch("PC.Hit",&PC.Hit);
  MainTree->Branch("Tr.NTracks",&Tr.NTracks,"NTracks/I");
  MainTree->Branch("Tr.TrackEvent", &Tr.TrEvent);
  MainTree->Branch("RFTime",&RFTime,"RFTime/I");
  MainTree->Branch("MCPTime",&MCPTime,"MCPTime/I");
  MainTree->Branch("IC",&IC,"IC/D");
  MainTree->Branch("E_IC",&E_IC,"E_IC/D");

  TObjArray *RootObjects = new TObjArray();
  RootObjects->Add(MainTree);

  //Initialize variables and channel map
  //-------------------------------------------------------------------------------------------------------------------------------------------
  TRandom3 *ChRandom = new TRandom3();
  //TRandom3 *PhiRandom = new TRandom3();
  //TRandom3 *RRandom = new TRandom3();

  fhlist = new TList;
  RootObjects->Add(fhlist);

  ChannelMap *CMAP;
  CMAP = new ChannelMap();
  //Initialization of the main channel map
  CMAP->Init("/home/manasta/Desktop/parker_codes/CalParamFiles/maria/ASICS_cmap_06292016",
	     "/home/manasta/Desktop/parker_codes/CalParamFiles/maria/alignchannels_09122016.dat",
	     "/home/manasta/Desktop/parker_codes/CalParamFiles/maria/AlphaCalibration_09132016.dat",
	     "/home/manasta/Desktop/parker_codes/CalParamFiles/maria/X3RelativeGains09202016_Step3_maskzero.dat",
	     "/home/manasta/Desktop/parker_codes/CalParamFiles/maria/QQQRelativeGains09182016_Step2.dat");
  CMAP->FinalInit("/home/manasta/Desktop/parker_codes/CalParamFiles/FinalFix012516.dat","/home/manasta/Desktop/parker_codes/CalParamFiles/maria/X3geometry_09132016.dat");
  CMAP->LoadQQQ3FinalFix("/home/manasta/Desktop/parker_codes/CalParamFiles/QQQ3FinalFix.012216");
  CMAP->InitWorldCoordinates("/home/manasta/Desktop/parker_codes/CalParamFiles/maria/NewWorld_030316.dat");
  
  CMAP->InitPCADC("/home/manasta/Desktop/parker_codes/CalParamFiles/maria/NewPCMap");
  CMAP->InitPCCalibration("/home/manasta/Desktop/parker_codes/CalParamFiles/maria/PCpulserCal2016July27.dat");
  CMAP->InitPCWireCal("/home/manasta/Desktop/parker_codes/CalParamFiles/maria/PCWireCal_030416.dat");

  //save the calibration information in the root file, so that the info doesn't get lost
  //This part is still under construction
  SETTINGS *set = new SETTINGS;
  TTree *settings_tree = new TTree("settings_tree","Settings Tree");
  settings_tree->Branch("settings_values",&set->set_val,"X3Pulser_Offset[24][12]/D:X3Pulser_Slope[24][12]/D:X3RelativeGains[24][12]/D:X3FinalFix[24][12]/D:X3Geometry_Up[24][4][4]/D:X3Geometry_Down[24][4][4]/D:QQQPulser_Offset[4][32]/D:QQQPulser_Slope[4][32]/D:QQQRelativeGains[4][32]/D:QQQFinalFix[4][32]/D:SiAlphaCal[28]/D");

  for (Int_t i=0; i<24; i++){
    for (Int_t j=0; j<12; j++){
      CMAP->GetZeroShift(i+4,j,set->set_val.X3Pulser_Offset[i][j],set->set_val.X3Pulser_Slope[i][j]);
      CMAP->GetX3MeVPerChannel1(i+4,j,set->set_val.X3RelativeGains[i][j]);   
      CMAP->GetX3FinalEnergyOffsetInMeV(i+4,j,set->set_val.X3FinalFix[i][j]);
    }
  }
  for (Int_t i=0; i<4; i++){
    for (Int_t j=0; j<32; j++){
      CMAP->GetZeroShift(i,j,set->set_val.QQQPulser_Offset[i][j],set->set_val.QQQPulser_Slope[i][j]);
      CMAP->GetQQQ3MeVPerChannel1(i,j,set->set_val.QQQRelativeGains[i][j]);
      CMAP->GetQQQ3FinalEnergyOffsetInMeV(i,j,set->set_val.QQQFinalFix[i][j]);
    }
  }
  for (Int_t i=0; i<28; i++){
    CMAP->GetX3MeVPerChannel2(i,0,set->set_val.SiAlphaCal[i]); 
  }
  for (Int_t i=0; i<24; i++){
    for (Int_t j=0;j<4; j++){
      for (Int_t k=0; k<4; k++){
	set->set_val.X3Geometry_Up[i][j][k] = CMAP->EdgeUp[i][j][k];
	set->set_val.X3Geometry_Down[i][j][k] = CMAP->EdgeDown[i][j][k];
      }
    }
  }

  settings_tree->Fill();

  TFile *inputFile = new TFile(filename_callist);//open root file and make sure it exists----------------------------------------------------------------------------
  if (!inputFile->IsOpen()){
    cout << "Root file: " << filename_callist << " could not be opened.\n";
    exit(EXIT_FAILURE);
  }
  //set branch addresses so that we can access the data tree--------------------------------------------------------------------------------------------------------
  ASICHit Si_Old;
  CAENHit ADC;
  CAENHit TDC;
  
  TTree *input_tree = (TTree*) inputFile->Get("DataTree");
  input_tree->SetBranchAddress("Si.Nhits",&Si_Old.Nhits);
  input_tree->SetBranchAddress("Si.MBID",Si_Old.MBID);
  input_tree->SetBranchAddress("Si.CBID",Si_Old.CBID);
  input_tree->SetBranchAddress("Si.ChNum",Si_Old.ChNum);
  input_tree->SetBranchAddress("Si.Energy",Si_Old.Energy);
  input_tree->SetBranchAddress("Si.Time",Si_Old.Time);
  
  input_tree->SetBranchAddress("ADC.Nhits",&ADC.Nhits);
  input_tree->SetBranchAddress("ADC.ID",ADC.ID);
  input_tree->SetBranchAddress("ADC.ChNum",ADC.ChNum);
  input_tree->SetBranchAddress("ADC.Data",ADC.Data);
  input_tree->SetBranchAddress("TDC.Nhits",&TDC.Nhits);
  input_tree->SetBranchAddress("TDC.ID",TDC.ID);
  input_tree->SetBranchAddress("TDC.ChNum",TDC.ChNum);
  input_tree->SetBranchAddress("TDC.Data",TDC.Data);
  //-----------------------------------------------------------------------------------------
  //define variables to hold data
  //I have separate arrays for each of the steps of the Silicon calibration process
  //These arrays are arranged so that the detector number is in the first dimension and the channel number is in the second dimension
  //   For the SX3s, the real detector number is the index + 4
  //   For QQQ it is one to one
  //If a detector/channel fired in an event, the corresponding row/column will be filled with the data
  //   If not, the value will be zero
  //After the arrays are filled, we loop over each row (detector number) and count the multiplicity of positive elements in that row
  //   This gives us the number of channels that fired in a detector (which we break into up/down/back multiplicity)
  Int_t DN = -1;
  Int_t DetCh=-1;
  
  Double_t ZeroShift=0,VperCh=0,slope_rel=0,slope_alpha=0,FinalShift=0;
  
  Double_t SiEnergy[NumX3][MaxX3Ch],SiTime[NumX3][MaxX3Ch];
  Double_t SiEnergy_Pulser[NumX3][MaxX3Ch],SiEnergy_Rel[NumX3][MaxX3Ch],SiEnergy_Cal[NumX3][MaxX3Ch];
  Double_t QQQ3Energy[NumQQQ3][MaxQQQ3Ch],QQQ3Time[NumQQQ3][MaxQQQ3Ch];
  Double_t QQQ3Energy_Pulser[NumQQQ3][MaxQQQ3Ch],QQQ3Energy_Rel[NumQQQ3][MaxQQQ3Ch],QQQ3Energy_Cal[NumQQQ3][MaxQQQ3Ch];
  
  ChRandom->SetSeed();  //always set the seed

  //PC place holder variables
  Int_t Side,WireID;
  Double_t PCDownstream[NPCWires],PCUpstream[NPCWires];
  Double_t PCDownVoltage[NPCWires],PCUpVoltage[NPCWires];
  Int_t ConvTest=0;
  Double_t Vcal;

  Int_t status = 0;
  Int_t counter=0;

  //create an output file to dump the excluded Entry# and the SiHit to check the consistency of the code due to the randomizations
  ofstream outfile;
  outfile.open("dumb_entries.txt");

  Long64_t nentries = input_tree->GetEntries();
  cout << nentries << endl;
  for (Long64_t global_evt=0; global_evt<nentries; global_evt++){//loop over all entries in tree---------------------------------------------------------------------------
    //cout << global_evt << endl;
    status = input_tree->GetEvent(global_evt);

    if (global_evt == TMath::Nint(0.01*nentries))  cout << " 1% through the data" << endl;
    if (global_evt == TMath::Nint(0.10*nentries))  cout << " 10% through the data" << endl;
    if (global_evt == TMath::Nint(0.15*nentries))  cout << " 15% through the data" << endl;
    if (global_evt == TMath::Nint(0.25*nentries))  cout << " 25% through the data" << endl;
    if (global_evt == TMath::Nint(0.35*nentries))  cout << " 35% through the data" << endl;
    if (global_evt == TMath::Nint(0.50*nentries))  cout << " 50% through the data" << endl;
    if (global_evt == TMath::Nint(0.65*nentries))  cout << " 65% through the data" << endl;
    if (global_evt == TMath::Nint(0.75*nentries))  cout << " 75% through the data" << endl;
    if (global_evt == TMath::Nint(0.90*nentries))  cout << " 90% through the data" << endl;
    if (global_evt == TMath::Nint(0.95*nentries))  cout << " 95% through the data" << endl;
    if (global_evt == TMath::Nint(1.00*nentries))  cout << " 100% through the data" << endl;

    //Initialize all variables at the beginning of the loop
    Si.zeroSiHit();
    PC.zeroPCHit();
    Tr.zeroTrack();

    for (int i=0; i<NumQQQ3; i++){
      for (int j=0; j<MaxQQQ3Ch; j++){
	QQQ3Energy[i][j] = 0;
	QQQ3Energy_Pulser[i][j] = 0;
	QQQ3Energy_Rel[i][j] = 0;
	QQQ3Energy_Cal[i][j] = 0;
	QQQ3Time[i][j] = 0;
      }
    }
    for (int i=0; i<NumX3; i++){
      for (int j=0; j<MaxX3Ch; j++){
	SiEnergy[i][j] = 0;
	SiEnergy_Pulser[i][j] = 0;
	SiEnergy_Rel[i][j] = 0;
	SiEnergy_Cal[i][j] = 0;
	SiTime[i][j] = 0;
      }
    }
    //Done initializing------------------------------------------------------------------------------------------------------------------------------------------------------------------
    //PC stuff
    Int_t Identifier=-1;
  
    //************************************************************************************************
    // Identifier provides information on the detector type recorded by ADC.
    // 0 - no identification was possible, this channel is not described in the ADCChannels.xxxx file
    // 1 - Gas Proportional counter data
    // 2 - CsI(Tl) scintillator
    // 3 - Ionization chamber at the downstream end
    // 4 - 16x16 detector at the end of the Ionization chamber for beam diagnostics
    //*************************************************************************************************

    // Initialize variables named above.
    for (Int_t i=0; i<NPCWires; i++){
      PCDownstream[i]  = 0;
      PCUpstream[i]    = 0;
      PCDownVoltage[i] = 0.0;
      PCUpVoltage[i]   = 0.0;
    }
    //  cout << "Event in" << endl;
      
    // Make sure your ADC.Nhits is within bounds.
    if (ADC.Nhits>MaxADCHits) ADC.Nhits=MaxADCHits;

    for (Int_t n=0; n<ADC.Nhits; n++){
      // Identify which detector type we have a hit in.
      CMAP->IdentifyADC(ADC.ID[n],ADC.ChNum[n],Identifier);
      switch (Identifier){
      case 0:
	break;
      case 1:
	CMAP->IdentifyWire(ADC.ID[n],ADC.ChNum[n],WireID,Side);
	ConvTest =  CMAP->ConvertToVoltage(ADC.ID[n],ADC.ChNum[n],ADC.Data[n],Vcal);
	if (Side==1) {
	  PCDownstream[WireID] = (Double_t)ADC.Data[n];
	  if (ConvTest) PCDownVoltage[WireID] = Vcal;
	}
	if (Side==2) {
	  PCUpstream[WireID]   = (Double_t)ADC.Data[n];
	  if (ConvTest) PCUpVoltage[WireID]  = Vcal;
	}
	//   cout << "ConvTest: " << ConvTest << endl;
	ConvTest = 0;
	break;
      case 2:
	break;
      case 3:
	break;
      case 4:
	break;
      default:
	break;
      }// End switch

    }// End loop over ADC.Nhits
    //cout << Si16En << " " << Ion0 << " " << Ion1 << " " << Ion2 <<endl;

    Double_t XWPC,YWPC,ZWPC,RWPC,PhiWPC;
    Double_t Energy8, Energy9;
    Int_t Good8 = 0, Good9 = 0;
    for ( Int_t i=0; i<NPCWires; i++ ){
      PC.zeroPlaceHolder();
      if ( PCDownstream[i]>100 && PCUpstream[i]>100 && PCDownstream[i]<4000 && PCUpstream[i]<4000 ){
	PC.pc_place_holder.WireID = i;
	PC.pc_place_holder.TrackType = 3;
	PC.pc_place_holder.Down = PCDownstream[i];
	PC.pc_place_holder.Up = PCUpstream[i];
	PC.pc_place_holder.DownVoltage = PCDownVoltage[i];
	PC.pc_place_holder.UpVoltage = PCUpVoltage[i];
	if (PCDownVoltage[i]>0 && PCUpVoltage[i]>0){
	  PC.pc_place_holder.Energy = PC.pc_place_holder.DownVoltage + PC.pc_place_holder.UpVoltage;
	  PC.pc_place_holder.Z = (PC.pc_place_holder.UpVoltage - PC.pc_place_holder.DownVoltage)/PC.pc_place_holder.Energy;
	  if ( PC.pc_place_holder.WireID == 8 ){
	    Energy8 = PC.pc_place_holder.Energy;
	    Good8 = 1;
	  }else if ( PC.pc_place_holder.WireID == 9 ){
	    Energy9 = PC.pc_place_holder.Energy;
	    Good9 = 1;
	  }
	}

	CMAP->GetPCWorldCoordinates(PC.pc_place_holder.WireID,PC.pc_place_holder.Z,XWPC,YWPC,ZWPC,RWPC,PhiWPC);
	PC.pc_place_holder.XW = XWPC;
	PC.pc_place_holder.YW = YWPC;
	PC.pc_place_holder.ZW = ZWPC;
	PC.pc_place_holder.RW = RWPC;
	PC.pc_place_holder.PhiW = PhiWPC;

	PC.Hit.push_back(PC.pc_place_holder);
	PC.NPCHits++;
      }
    }
    if (Good8 == 1 && Good9 == 1){
      MyFill("PCEnergy_Wire8_vs_9",300,0,1,Energy8,300,0,1,Energy9);
    }
 

    //--------------MCP/RF TIMING------------------------
   
    if (TDC.Nhits>MaxTDCHits) TDC.Nhits=MaxTDCHits;
  

    Double_t correct=1.004009623; //by M.A.
    RFTime = 0;
    MCPTime = 0;
  
  
    for (Int_t n=0; n<TDC.Nhits; n++){
      //TDC ID == 2?
      if( TDC.ID[n]==2 ){
	if (TDC.ChNum[n]==0)  RFTime   = TDC.Data[n];
	if (TDC.ChNum[n]==7)  MCPTime  = TDC.Data[n];
      }
    }


    MyFill("Timing",600,1,600,fmod((MCPTime*correct-RFTime),546)); //by M.A.

    
    //MyFill("Timing",400,-600,600,(MCPTime-RFTime)%538);
    //Can define the TimingCut at the top of the program
#ifdef TimingCut   
      if ( (MCPTime-RFTime)%538<71 || (MCPTime-RFTime)%538>370 ){
	continue;
      }
      if ( (MCPTime-RFTime)%538>100 && (MCPTime-RFTime)%538<330 ){
	continue;
      }
#endif  

//-------------------------------------
//------------Ion Chamber----------------------
      IC = 0; E_IC = 0;
      for(Int_t n=0; n<ADC.Nhits; n++)
	{ if(ADC.ID[n]==3 && ADC.ChNum[n]==24)
	    {
	      IC = (Double_t)ADC.Data[n];
	      //cout << IC << endl;
	    }
	}

      for(Int_t n=0; n<ADC.Nhits; n++)
	{ if(ADC.ID[n]==3 && ADC.ChNum[n]==28)
	    {
	      E_IC = (Double_t)ADC.Data[n];
	      //cout << IC << endl;
	    }
	}
      
   


    //end PC stuff et al-----------------------------------------------------------------------------------------------------------------
    // Make sure you dont have too many hits.
    if (Si_Old.Nhits>MaxSiHits) Si_Old.Nhits=MaxSiHits;
    //cout << "New Event\n";
    for (Int_t n=0; n<Si_Old.Nhits; n++){//loop over all Si hits--very similar to Main.C from Jeff
      ZeroShift  = 0;
      VperCh     = 1;
      slope_rel  = 1;
      slope_alpha = 1;
      FinalShift = 0;

      //cout << Si_Old.MBID[n] << "  " << Si_Old.CBID[n] << "  " << Si_Old.ChNum[n] << endl;

      CMAP->IdentifyDetChan(Si_Old.MBID[n],Si_Old.CBID[n],Si_Old.ChNum[n], DN, DetCh);
      CMAP->AlignASICsChannels(Si_Old.MBID[n],Si_Old.CBID[n],Si_Old.ChNum[n], ZeroShift, VperCh);
      //this is where we fill the placeholder arrays that we defined above
      if(DN>27) {continue;}
      if(DN>3){
	CMAP->GetX3MeVPerChannel1(DN,DetCh,slope_rel);  
	CMAP->GetX3MeVPerChannel2(DN,DetCh,slope_alpha);  
	CMAP->GetX3FinalEnergyOffsetInMeV(DN,DetCh,FinalShift);
	SiEnergy[DN-4][DetCh] = (Double_t)Si_Old.Energy[n];     // includes the RAW DATA 
	//SiEnergy_Pulser[DN-4][DetCh] = (Double_t)Si_Old.Energy[n]+ChRandom->Rndm()-0.5+ZeroShift/VperCh;   //// the randomization was not correct
	SiEnergy_Pulser[DN-4][DetCh] = (Double_t)Si_Old.Energy[n] + ZeroShift/VperCh; 
       	SiEnergy_Rel[DN-4][DetCh] = SiEnergy_Pulser[DN-4][DetCh]*slope_rel;
	SiEnergy_Cal[DN-4][DetCh] = SiEnergy_Rel[DN-4][DetCh]*slope_alpha;
	SiEnergy_Cal[DN-4][DetCh] += FinalShift;
	SiTime[DN-4][DetCh]   = (Double_t)Si_Old.Time[n];
      }else{
	CMAP->GetQQQ3MeVPerChannel1(DN,DetCh,slope_rel);
	CMAP->GetQQQ3MeVPerChannel2(DN,DetCh,slope_alpha);
	CMAP->GetQQQ3FinalEnergyOffsetInMeV(DN,DetCh,FinalShift);
	QQQ3Energy[DN][DetCh] = (Double_t)Si_Old.Energy[n];
	//QQQ3Energy_Pulser[DN][DetCh] = (Double_t)Si_Old.Energy[n]+ChRandom->Rndm()-0.5+ZeroShift/VperCh;   //// the randomization was not correct
	QQQ3Energy_Pulser[DN][DetCh] = (Double_t)Si_Old.Energy[n] + ZeroShift/VperCh;
	QQQ3Energy_Rel[DN][DetCh] = QQQ3Energy_Pulser[DN][DetCh]*slope_rel;
	QQQ3Energy_Cal[DN][DetCh] = QQQ3Energy_Rel[DN][DetCh]*slope_alpha;
	QQQ3Energy_Cal[DN][DetCh] += FinalShift;
	QQQ3Time[DN][DetCh] = (Double_t)Si_Old.Time[n];
      }
    }// End loop over X3_Nhits------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /////////////////-----------------------------SX3s------------------------------------------------------------------------------
    for (int i=0; i<NumX3; i++){//Determine multiplicities for SX3 detector: forward and back--------------------------------------------------------------------------------------------------
      //loop over all SX3s and count up/down/back multiplicities and fill the detector place holder if we have an energy>0
      Si.zeroPlaceHolder();
      for (int j=0; j<MaxX3Ch; j++){
	//if ( (SiEnergy_Cal[i][j] > 0.0)){
	  if ( (SiEnergy[i][j] > 0.0)){                           //Fill when the Raw Energy is >0  
	  if (j<4){
	    Si.det_place_holder.BackChNum.push_back(j);
	    Si.det_place_holder.EnergyBack_Raw.push_back(SiEnergy[i][j]);
	    Si.det_place_holder.EnergyBack_Pulser.push_back(SiEnergy_Pulser[i][j]);
	    Si.det_place_holder.EnergyBack_Rel.push_back(SiEnergy_Rel[i][j]);
	    Si.det_place_holder.EnergyBack_Cal.push_back(SiEnergy_Cal[i][j]);
	    Si.det_place_holder.TimeBack.push_back(SiTime[i][j]);
	  }else if(j>3 && j<8){
	    Si.det_place_holder.UpChNum.push_back(j-4);
	    Si.det_place_holder.EnergyUp_Raw.push_back(SiEnergy[i][j]);
	    Si.det_place_holder.EnergyUp_Pulser.push_back(SiEnergy_Pulser[i][j]);
	    Si.det_place_holder.EnergyUp_Rel.push_back(SiEnergy_Rel[i][j]);
	    Si.det_place_holder.EnergyUp_Cal.push_back(SiEnergy_Cal[i][j]);
	    Si.det_place_holder.TimeUp.push_back(SiTime[i][j]);
	  }else if(j>7 && j<12){
	    Si.det_place_holder.DownChNum.push_back(j-8);
	    Si.det_place_holder.EnergyDown_Raw.push_back(SiEnergy[i][j]);
	    Si.det_place_holder.EnergyDown_Pulser.push_back(SiEnergy_Pulser[i][j]);
	    Si.det_place_holder.EnergyDown_Rel.push_back(SiEnergy_Rel[i][j]);
	    Si.det_place_holder.EnergyDown_Cal.push_back(SiEnergy_Cal[i][j]);
	    Si.det_place_holder.TimeDown.push_back(SiTime[i][j]);
	  }
	}
      }

      if ( Si.det_place_holder.EnergyDown_Cal.size()!=0 || Si.det_place_holder.EnergyUp_Cal.size()!=0 || Si.det_place_holder.EnergyBack_Cal.size()!=0 ){
	//cout << "Sorted Data\n";
	Si.det_place_holder.UpMult = Si.det_place_holder.EnergyUp_Cal.size();
	Si.det_place_holder.DownMult = Si.det_place_holder.EnergyDown_Cal.size();
	Si.det_place_holder.BackMult = Si.det_place_holder.EnergyBack_Cal.size();
	Si.det_place_holder.DetID = i+4;
	//cout << Si.det_place_holder.DetID << endl;
	Si.det_place_holder.HitType = Si.det_place_holder.BackMult*100 + Si.det_place_holder.UpMult*10 + Si.det_place_holder.DownMult;

	SiSort.ProcessSX3_2(&Si,CMAP,RFTime);//ProcessSX3_2 treats them case by case (hit type 111, 110, 101, 211)---------SX3_1 treats them like the QQQs in a general manner (but not very well)

	Si.Detector.push_back(Si.det_place_holder);
	Si.Hit.push_back(Si.hit_place_holder);

	Si.NSiHits++;
	

	//---Histos after the energy calibration

	/*if (Si.hit_place_holder.HitType==111){
	  MyFill(Form("back_vs_front%i",Si.hit_place_holder.DetID),300,0,30,Si.hit_place_holder.EnergyFront,300,0,30,Si.hit_place_holder.EnergyBack);
	  MyFill(Form("back_vs_front%i_front%i",Si.det_place_holder.DetID,Si.det_place_holder.UpChNum[0]),300,0,30,Si.det_place_holder.EnergyUp_Cal[0]+Si.det_place_holder.EnergyDown_Cal[0],300,0,30,Si.det_place_holder.EnergyBack_Cal[0]);
	  MyFill(Form("back_vs_front%i_back%i",Si.det_place_holder.DetID,Si.det_place_holder.BackChNum[0]),300,0,30,Si.det_place_holder.EnergyUp_Cal[0]+Si.det_place_holder.EnergyDown_Cal[0],300,0,30,Si.det_place_holder.EnergyBack_Cal[0]);
	  MyFill(Form("back_vs_front%i_%i_%i",Si.det_place_holder.DetID,Si.det_place_holder.UpChNum[0],Si.det_place_holder.BackChNum[0]),300,0,30,Si.det_place_holder.EnergyUp_Cal[0]+Si.det_place_holder.EnergyDown_Cal[0],300,0,30,Si.det_place_holder.EnergyBack_Cal[0]);
	  MyFill(Form("down_vs_up%i_front%i",Si.det_place_holder.DetID,Si.det_place_holder.UpChNum[0]),300,0,30,Si.det_place_holder.EnergyUp_Cal[0],300,0,30,Si.det_place_holder.EnergyDown_Cal[0]);
	  }*/


	//---Histos before energy calibration
	
	if (Si.hit_place_holder.HitType==111){
	  MyFill(Form("back_vs_front%i",Si.hit_place_holder.DetID),512,0,16384,Si.hit_place_holder.EnergyFront,512,0,16384,Si.hit_place_holder.EnergyBack);
	  MyFill(Form("back_vs_front%i_front%i",Si.det_place_holder.DetID,Si.det_place_holder.UpChNum[0]),512,0,16384,Si.det_place_holder.EnergyUp_Cal[0]+Si.det_place_holder.EnergyDown_Cal[0],512,0,16384,Si.det_place_holder.EnergyBack_Cal[0]);
	  MyFill(Form("back_vs_front%i_back%i",Si.det_place_holder.DetID,Si.det_place_holder.BackChNum[0]),512,0,16384,Si.det_place_holder.EnergyUp_Cal[0]+Si.det_place_holder.EnergyDown_Cal[0],512,0,16384,Si.det_place_holder.EnergyBack_Cal[0]);
	  MyFill(Form("back_vs_front%i_%i_%i",Si.det_place_holder.DetID,Si.det_place_holder.UpChNum[0],Si.det_place_holder.BackChNum[0]),512,0,16384,Si.det_place_holder.EnergyUp_Cal[0]+Si.det_place_holder.EnergyDown_Cal[0],512,0,16384,Si.det_place_holder.EnergyBack_Cal[0]);
	  MyFill(Form("down_vs_up%i_front%i",Si.det_place_holder.DetID,Si.det_place_holder.UpChNum[0]),512,0,16384,Si.det_place_holder.EnergyUp_Cal[0],512,0,16384,Si.det_place_holder.EnergyDown_Cal[0]);
	  MyFill(Form("down_vs_up_divideBack%i_front%i",Si.det_place_holder.DetID,Si.det_place_holder.UpChNum[0]),100,0,1,(Si.det_place_holder.EnergyUp_Cal[0]/Si.det_place_holder.EnergyBack_Cal[0]),100,0,1,(Si.det_place_holder.EnergyDown_Cal[0]/Si.det_place_holder.EnergyBack_Cal[0]));
	  }
	
	for ( Int_t hits=0; hits<PC.NPCHits; hits++ ){
	  MyFill("PCPhi_vs_SiPhi_SX3",300,0,8,Si.hit_place_holder.PhiW,300,0,8,PC.Hit.at(hits).PhiW);
	}

      }
    }

    //////////////////////////--------------------------QQQs-----------------------------------------------------------------------------------------------
    for (int i=0; i<NumQQQ3; i++){//Determine multiplicities for QQQ detector: forward and back-------------------------------------------------------------------------------------------------------
     //loop over all QQQs and count up (front)/back multiplicities and fill the detector place holder if we have an energy>0
      Si.zeroPlaceHolder();
      for (int j=0; j<MaxQQQ3Ch; j++){
	//if ( (QQQ3Energy_Cal[i][j] > 0)){
	  if ( (QQQ3Energy[i][j] > 0)){
	  if (j<16){
	    Si.det_place_holder.BackChNum.push_back(j);
	    Si.det_place_holder.EnergyBack_Raw.push_back(QQQ3Energy[i][j]);
	    Si.det_place_holder.EnergyBack_Pulser.push_back(QQQ3Energy_Pulser[i][j]);
	    Si.det_place_holder.EnergyBack_Rel.push_back(QQQ3Energy_Rel[i][j]);
	    Si.det_place_holder.EnergyBack_Cal.push_back(QQQ3Energy_Cal[i][j]);
	    Si.det_place_holder.TimeBack.push_back(QQQ3Time[i][j]);
	    if (i==1 && j==14){
	      MyFill("PCEnergy_Wire8_vs_9_cut",300,0,1,Energy8,300,0,1,Energy9);
	    }
	  }else if(j>15){
	    Si.det_place_holder.UpChNum.push_back(j-16);
	    Si.det_place_holder.EnergyUp_Raw.push_back(QQQ3Energy[i][j]);
	    Si.det_place_holder.EnergyUp_Pulser.push_back(QQQ3Energy_Pulser[i][j]);
	    Si.det_place_holder.EnergyUp_Rel.push_back(QQQ3Energy_Rel[i][j]);
	    Si.det_place_holder.EnergyUp_Cal.push_back(QQQ3Energy_Cal[i][j]);
	    Si.det_place_holder.TimeUp.push_back(QQQ3Time[i][j]);
	  }
	}
      }
      Si.det_place_holder.UpMult = Si.det_place_holder.EnergyUp_Cal.size();
      Si.det_place_holder.BackMult = Si.det_place_holder.EnergyBack_Cal.size();

      Si.det_place_holder.DetID = i;
      Si.det_place_holder.HitType = Si.det_place_holder.BackMult*10 + Si.det_place_holder.UpMult;
      if ( Si.det_place_holder.EnergyUp_Cal.size()!=0 && Si.det_place_holder.EnergyBack_Cal.size()!=0 ){
	//we require both front and back of qqq to fire to sort the data
	//cout << "Sorted Data\n";

	//ProcessQQQ_1 treats QQQ in a general way (recommended)
	//ProcessQQQ_2 only does hit type 11, 12, and 21 (assuming adjacent channels fired)
	SiSort.ProcessQQQ_2(&Si,CMAP,RFTime);
	Si.Detector.push_back(Si.det_place_holder);
	Si.Hit.push_back(Si.hit_place_holder);
	Si.NSiHits++;

	//---Histos after the energy calibration
	
	/*if (Si.det_place_holder.HitType==11){
	  MyFill(Form("back_vs_front%i",Si.hit_place_holder.DetID),300,0,30,Si.hit_place_holder.EnergyFront,300,0,30,Si.hit_place_holder.EnergyBack);
	  MyFill(Form("back_vs_front%i_front%i",Si.det_place_holder.DetID,Si.det_place_holder.UpChNum[0]),300,0,30,Si.det_place_holder.EnergyUp_Cal[0],300,0,30,Si.det_place_holder.EnergyBack_Cal[0]);
	  MyFill(Form("back_vs_front%i_back%i",Si.det_place_holder.DetID,Si.det_place_holder.BackChNum[0]),300,0,30,Si.det_place_holder.EnergyUp_Cal[0],300,0,30,Si.det_place_holder.EnergyBack_Cal[0]);
	  MyFill(Form("back_vs_front%i_%i_%i",Si.det_place_holder.DetID,Si.det_place_holder.UpChNum[0],Si.det_place_holder.BackChNum[0]),300,0,30,Si.det_place_holder.EnergyUp_Cal[0],300,0,30,Si.det_place_holder.EnergyBack_Cal[0]);
	  }*/

	//---Histos before the energy calibration

	if (Si.det_place_holder.HitType==11){
	  MyFill(Form("back_vs_front%i",Si.hit_place_holder.DetID),512,0,16384,Si.hit_place_holder.EnergyFront,512,0,16384,Si.hit_place_holder.EnergyBack);
	  MyFill(Form("back_vs_front%i_front%i",Si.det_place_holder.DetID,Si.det_place_holder.UpChNum[0]),512,0,16384,Si.det_place_holder.EnergyUp_Cal[0],512,0,16384,Si.det_place_holder.EnergyBack_Cal[0]);
	  MyFill(Form("back_vs_front%i_back%i",Si.det_place_holder.DetID,Si.det_place_holder.BackChNum[0]),512,0,16384,Si.det_place_holder.EnergyUp_Cal[0],512,0,16384,Si.det_place_holder.EnergyBack_Cal[0]);
	  MyFill(Form("back_vs_front%i_%i_%i",Si.det_place_holder.DetID,Si.det_place_holder.UpChNum[0],Si.det_place_holder.BackChNum[0]),512,0,16384,Si.det_place_holder.EnergyUp_Cal[0],512,0,16384,Si.det_place_holder.EnergyBack_Cal[0]);
	  }

	for ( Int_t hits=0; hits<PC.NPCHits; hits++ ){
	  MyFill("PCPhi_vs_SiPhi_QQQ",300,0,8,Si.hit_place_holder.PhiW,300,0,8,PC.Hit.at(hits).PhiW);
	}

	
      }else if ( Si.det_place_holder.EnergyUp_Cal.size()!=0 || Si.det_place_holder.EnergyBack_Cal.size()!=0 ){
	//if up or down fired, we write the unsorted data to the tree, but we don't process the hit
	Si.Detector.push_back(Si.det_place_holder);
	Si.Hit.push_back(Si.hit_place_holder);
	Si.NSiHits++;
      }
    }
    
    if (Si.NSiHits != Si.Detector.size()){//these should always equal. If not, something is wrong
      cout << Si.NSiHits << "  " << Si.Detector.size() << endl;
    }
    if (Si.NSiHits != Si.Hit.size()){
      cout << Si.NSiHits << "  " << Si.Detector.size() << endl;
    }

    ////////////////////----------Tracking--------------------------------------------------------------------------------------------
#ifdef Tracking
    //associate Silicon hits with PC hits
    Int_t GoodPC = -1;
    Double_t m,b;
    for ( Int_t hits=0; hits<Si.NSiHits; hits++ ){
      //Find PC close in phi to Silicon that fired
      //If more than one found, take maximum signal
      //Mask PC hit by setting PC.Hit.TrackType = 1
      GoodPC = FindMaxPC(Si.Hit.at(hits).PhiW,PC);
      Tr.zeroPlaceHolder();
      if (GoodPC > -1){//if a PC is found do stuff
	//fill Si parameters
	Si.Hit.at(hits).TrackType = 1;
	PC.Hit.at(GoodPC).TrackType = 1;
	Tr.track_place_holder.TrackType = 1;
	Tr.track_place_holder.HitType = Si.Hit.at(hits).HitType;
	Tr.track_place_holder.DetID = Si.Hit.at(hits).DetID;
	Tr.track_place_holder.SiZ = Si.Hit.at(hits).ZW;
	Tr.track_place_holder.SiR = Si.Hit.at(hits).RW;
	Tr.track_place_holder.SiPhi = Si.Hit.at(hits).PhiW;
	Tr.track_place_holder.SiEnergy = Si.Hit.at(hits).Energy;
	//fill PC parameters
	Tr.track_place_holder.PCZ = PC.Hit.at(GoodPC).ZW;
	Tr.track_place_holder.PCR = PC.Hit.at(GoodPC).RW;
	Tr.track_place_holder.PCPhi = PC.Hit.at(GoodPC).PhiW;
	Tr.track_place_holder.PCEnergy = PC.Hit.at(GoodPC).Energy;
	Tr.track_place_holder.WireID = PC.Hit.at(GoodPC).WireID;
	//reconstruct interaction point
	m = (Tr.track_place_holder.PCR-Tr.track_place_holder.SiR)/(Tr.track_place_holder.PCZ-Tr.track_place_holder.SiZ);
	b = Tr.track_place_holder.PCR - m*Tr.track_place_holder.PCZ;
	Tr.track_place_holder.IntPoint = -b/m;
	if ( (Tr.track_place_holder.IntPoint - Tr.track_place_holder.SiZ) > 0){
	  Tr.track_place_holder.Theta = atan(Tr.track_place_holder.SiR/(Tr.track_place_holder.IntPoint - Tr.track_place_holder.SiZ));
	  Tr.track_place_holder.PathLength = Tr.track_place_holder.SiR/sin(Tr.track_place_holder.Theta);
	}else if ( (Tr.track_place_holder.IntPoint - Tr.track_place_holder.SiZ) < 0){
	  Tr.track_place_holder.Theta = atan(Tr.track_place_holder.SiR/(Tr.track_place_holder.IntPoint - Tr.track_place_holder.SiZ)) + TMath::Pi();
	  Tr.track_place_holder.PathLength = Tr.track_place_holder.SiR/sin(Tr.track_place_holder.Theta);
	}else{
	  Tr.track_place_holder.Theta = TMath::Pi()/2;
	  Tr.track_place_holder.PathLength = Tr.track_place_holder.SiR;
	}
	

	//--------------------E_dE histograms--------------------------------------------------------
	
       	for(int i=0;i<NumQQQ3; i++){
	MyFill("E_de_QQQ",300,0,30,Tr.track_place_holder.SiEnergy,300,0,1,Tr.track_place_holder.PCEnergy*sin(Tr.track_place_holder.Theta));
	MyFill("E_theta_QQQ",300,0,180,Tr.track_place_holder.Theta*180/M_PI,300,0,35,Tr.track_place_holder.SiEnergy);
	}
	
	
	for(int i=0;i<NumX3; i++){
	MyFill("E_de_SX3",300,0,30,Tr.track_place_holder.SiEnergy,300,0,1,Tr.track_place_holder.PCEnergy*sin(Tr.track_place_holder.Theta));
	MyFill("E_theta_SX3",300,0,180,Tr.track_place_holder.Theta*180/M_PI,300,0,35,Tr.track_place_holder.SiEnergy);
	}


	//if you are doing a proton or alpha cal run, then you need to input the proper gold position here
	//this calculates where the PC should have fired, based off of the gold and silicon positions
	//compare where it should have fired to the measured Z to calibrate

#ifdef PC_Pos_Cal
	Double_t gold_pos = 28.5; //all the way in
	//Double_t gold_pos = 24.156; // spacer 1= 4.8514 cm   //these distances are wrong, have to redo it M.A.
	//Double_t gold_pos = 18.085; // spacer 2= 10.8712 cm
	//Double_t gold_pos = 16.815; // spacer 3= 12.1412 cm
	//Double_t gold_pos = 13.665; // spacer 4= 15.2908 cm
	//Double_t gold_pos = 8.611; // spacer 5= 20.3454 cm
	Double_t mpc = (Tr.track_place_holder.SiZ - gold_pos)/Tr.track_place_holder.SiR;
	Double_t bpc = Tr.track_place_holder.SiZ - mpc*Tr.track_place_holder.SiR;

	Tr.track_place_holder.PCZ_Ref = mpc*3.75 + bpc;
	MyFill(Form("PCZ_vs_Z_nocal%i",Tr.track_place_holder.WireID),300,0,50,Tr.track_place_holder.PCZ_Ref,300,-1,1,PC.Hit.at(GoodPC).Z);
#endif

	Tr.NTracks++;
      }
      Tr.TrEvent.push_back(Tr.track_place_holder);

    }
#endif

    //////////////-----------------end of tracking---------------------------------------/////////////////////////////

    
    
    if (Si.NSiHits>0 ){    
      MainTree->Fill();
    }
    else{
      
      //---Dump the excluded events in a file:
          outfile << "event number = " << global_evt << " " << "Si_Old.Nhits= " << Si_Old.Nhits << endl;
          counter++; 
        }
   
  }
  
    outfile.close(); 
    cout << "event not count = " << counter << endl;
    cout << "entries = " << nentries << endl;
    cout << "hist = " << nentries - counter << endl;

outputFile->cd();
cout << "change to output file directory" << endl;
//fhlist->Write();
RootObjects->Write();
cout << "Write Root Objects" << endl;
outputFile->Close();
  cout << "Close Outputfile\n";
  return 0;
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
    fhmap.at(name)->Fill(valueX,valueY);
  } catch(out_of_range e) {
    TH2F* newHist = new TH2F(name.c_str(),name.c_str(),
			     binsX,lowX,highX,
			     binsY,lowY,highY);
    newHist->Fill(valueX,valueY);
    fhlist->Add(newHist);
    fhmap[name] = newHist;
  }
}

Int_t FindMaxPC(Double_t phi, PCHit& PC){
  Int_t GoodPC = -1;
  Double_t MaxPC = -10;
  Double_t MinPhi = 0.2619;
  for (int k=0; k<PC.NPCHits; k++){//loop over the pc hits
    //if the PC falls in a range of phi then it is possible correlated
    //we find the maximum energy on the pc
    PC.pc_place_holder = PC.Hit.at(k);

    if (PC.pc_place_holder.TrackType == 1){
      continue;
    }

    if ( (fabs(PC.pc_place_holder.PhiW-phi) <= MinPhi) || ((2*TMath::Pi() - fabs(PC.pc_place_holder.PhiW-phi)) <= MinPhi) ) {
      if ( PC.pc_place_holder.Energy >= MaxPC ){
	MaxPC = PC.pc_place_holder.Energy;
	GoodPC = k;
      }
    }
  }
  return GoodPC;
}
