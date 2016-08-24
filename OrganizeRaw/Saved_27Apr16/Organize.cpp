//Adapted Jeff's Main.C to compile with g++
//Usage: g++ -o ParkerMain ParkerMain.cpp `root-config --cflags --glibs`
//      ./ParkerMain DataList.txt outputfile.root
//DataList.txt contains a list of root files with a DataTree

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
#include <algorithm>
//#include <vector>

#include "/home2/parker/ANASEN/LSU/Include/ChannelMap.h"
#include "/home2/parker/ANASEN/LSU/Include/2015_detclass.h"
#include "/home2/parker/ANASEN/LSU/Include/organizetree.h"
//#include "organize_dictionary.h"
#include <TObjArray.h>
#include <TList.h>
#include <TFile.h>
//#include <TROOT.h>
#include <TTree.h>
#include <TApplication.h>

#define MaxPCHits 24
#define NPCWires  24

using namespace std;

struct data {
  Double_t Energy;
  Double_t Time;
  Double_t Channel;
};

bool Esort_method(struct data a,struct data b){
  if(a.Energy>b.Energy)
    return 1;
  return 0;
};

bool Csort_method(struct data a,struct data b){
  if(a.Channel<b.Channel)
    return 1;
  return 0;
};

void MyFill(string name,
	    int binsX, double lowX, double highX, double valueX,
	    int binsY, double lowY, double highY, double valueY);
void MyFill(string name,
	    int binsX, double lowX, double highX, double valueX);

void ProcessQQQ(SiHit *Si,ChannelMap *CMAP, Int_t RFTime);

TList* fhlist;
std::map<string,TH1*> fhmap;

int main(int argc, char* argv[]){
  
  TApplication *myapp=new TApplication("myapp",0,0); //Don't know what this does, but libraries won't load without it
  SiHit Si;
  PCHit PC;
  Track Tr;
  Int_t RFTime,MCPTime;

  //read in command line arguments
  char* filename_histout = new char [100];//output root file
  char* filename_callist = new char [100];//input file list
  
  strcpy( filename_callist, argv[1] );
  strcpy( filename_histout, argv[2] );

  TFile *outputFile = new TFile(filename_histout,"RECREATE");
  TTree *MainTree = new TTree("MainTree","MainTree");
  
  MainTree->Branch("Si.NSiHits",&Si.NSiHits,"NSiHits/I");
  MainTree->Branch("Si.Detector",&Si.Detector);
  MainTree->Branch("Si.Hit",&Si.Hit);

  MainTree->Branch("PC.NPCHits",&PC.NPCHits,"NPCHits/I");
  MainTree->Branch("PC.WireID",PC.WireID,"WireID[NPCHits]/I");
  MainTree->Branch("PC.Down",PC.Down,"Down[NPCHits]/D");
  MainTree->Branch("PC.Up",PC.Up,"Up[NPCHits]/D");
  MainTree->Branch("PC.DownVoltage",PC.DownVoltage,"DownVoltage[NPCHits]/D");
  MainTree->Branch("PC.UpVoltage",PC.UpVoltage,"UpVoltage[NPCHits]/D");
  MainTree->Branch("PC.Energy",PC.Energy,"Energy[NPCHits]/D");
  MainTree->Branch("PC.Zp",PC.Zp,"Zp[NPCHits]/D");
  MainTree->Branch("PC.XWp",PC.XWp,"XWp[NPCHits]/D");
  MainTree->Branch("PC.YWp",PC.YWp,"YWp[NPCHits]/D");
  MainTree->Branch("PC.ZWp",PC.ZWp,"ZWp[NPCHits]/D");
  MainTree->Branch("PC.RWp",PC.RWp,"RWp[NPCHits]/D");
  MainTree->Branch("PC.PhiWp",PC.PhiWp,"PhiWp[NPCHits]/D");

  MainTree->Branch("RFTime",&RFTime,"RFTime/I");
  MainTree->Branch("MCPTime",&MCPTime,"MCPTime/I");


  TObjArray *RootObjects = new TObjArray();
  RootObjects->Add(MainTree);

  //Initialize variables and channel map
  //-------------------------------------------------------------------------------------------------------------------------------------------
  TRandom3 *XRandom = new TRandom3();
  TRandom3 *ChRandom = new TRandom3();
  //TRandom3 *PhiRandom = new TRandom3();
  //TRandom3 *RRandom = new TRandom3();
  TRandom3 *ZRandom = new TRandom3();

  fhlist = new TList;
  RootObjects->Add(fhlist);

  SETTINGS *set = new SETTINGS;
  TTree *settings_tree = new TTree("settings_tree","Settings Tree");
  settings_tree->Branch("settings_values",&set->set_val,"X3Pulser_Offset[24][12]/D:X3Pulser_Slope[24][12]/D:X3RelativeGains[24][12]/D:X3FinalFix[24][12]/D:X3Geometry_Up[24][4][4]/D:X3Geometry_Down[24][4][4]/D:QQQPulser_Offset[4][32]/D:QQQPulser_Slope[4][32]/D:QQQRelativeGains[4][32]/D:QQQFinalFix[4][32]/D:SiAlphaCal[28]/D");
  RootObjects->Add(settings_tree);
  ChannelMap *CMAP;
  CMAP = new ChannelMap();
  //Initialization of the main channel map
  //CMAP->Init("/home2/parker/ANASEN/LSU/CalParamFiles/ASICS_cmap_022716");
  //CMAP->InitPCADC("/home2/parker/ANASEN/LSU/CalParamFiles/NewPCMap");
  CMAP->Init("/home2/parker/ANASEN/LSU/CalParamFiles/ASICS_cmap_022716","/home2/parker/ANASEN/LSU/CalParamFiles/alignchannels_012216.txt","/home2/parker/ANASEN/LSU/CalParamFiles/AlphaCalibration_022716.dat","/home2/parker/ANASEN/LSU/CalParamFiles/X3RelativeGains031516.dat","/home2/parker/ANASEN/LSU/CalParamFiles/QQQRelativeGains020216.dat");//most updated cal files 03/07/2016
  //CMAP->Init("/home2/parker/ANASEN/LSU/CalParamFiles/ASICS_cmap_022716","/home2/parker/ANASEN/LSU/CalParamFiles/alignchannels_012216.txt","/home2/parker/ANASEN/LSU/CalParamFiles/AlphaCalibration_022716.dat","/home2/parker/ANASEN/LSU/CalParamFiles/X3RelativeGains_Slope1.dat","/home2/parker/ANASEN/LSU/CalParamFiles/QQQRelativeGains020216.dat");
   //CMAP->FinalInit("/home2/parker/ANASEN/LSU/CalParamFiles/FinalFix012516.dat","/home2/parker/ANASEN/LSU/CalParamFiles/X3geometry_020416.dat");
  CMAP->FinalInit("/home2/parker/ANASEN/LSU/CalParamFiles/FinalFix012516.dat","/home2/parker/ANASEN/LSU/CalParamFiles/X3geometry_032816.dat");
  CMAP->LoadQQQ3FinalFix("/home2/parker/ANASEN/LSU/CalParamFiles/QQQ3FinalFix.012216");
  CMAP->InitWorldCoordinates("/home2/parker/ANASEN/LSU/CalParamFiles/NewWorld_030316.dat");
  
  CMAP->InitPCADC("/home2/parker/ANASEN/LSU/CalParamFiles/NewPCMap");
  
  //Mesytec Shaper
  CMAP->InitPCCalibration("/home2/parker/ANASEN/LSU/CalParamFiles/PCPulser042016.dat");
  CMAP->InitPCWireCal("/home2/parker/ANASEN/LSU/CalParamFiles/PCWireCal_030416.dat");

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

  //open up list of root files
  ifstream inFileList;
  inFileList.open(filename_callist);
  if (!inFileList.is_open()){
    cout << "In File List not open\n";
    exit(EXIT_FAILURE);
  }
  string rootfile;
  char rootfile_char[100];

  Int_t mult111_bad = 0;

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

    Int_t DN = -1;
    Int_t DetCh=-1;

    Double_t ZeroShift=0,VperCh=0,slope_rel=0,slope_alpha=0,FinalShift=0;

    Double_t SiEnergy[NumX3][MaxX3Ch],SiTime[NumX3][MaxX3Ch];
    Double_t SiEnergy_Pulser[NumX3][MaxX3Ch],SiEnergy_Rel[NumX3][MaxX3Ch],SiEnergy_Cal[NumX3][MaxX3Ch];
    Double_t QQQ3Energy[NumQQQ3][MaxQQQ3Ch],QQQ3Time[NumQQQ3][MaxQQQ3Ch];
    Double_t QQQ3Energy_Pulser[NumQQQ3][MaxQQQ3Ch],QQQ3Energy_Rel[NumQQQ3][MaxQQQ3Ch],QQQ3Energy_Cal[NumQQQ3][MaxQQQ3Ch];

    XRandom->SetSeed();
    ChRandom->SetSeed();
    //PhiRandom->SetSeed();
    //RRandom->SetSeed();
    ZRandom->SetSeed();

    //Double_t QQQR;
    //Double_t QQQPhi;

    //const Double_t InnerRadius = 5.01; //cm4.81
    //const Double_t OuterRadius = 10.1; //cm
    //const Double_t ArcLength = 0.96; //cm (distance between strips at outer radius) 
    //const Double_t StripAngle = 2*TMath::ASin(ArcLength/OuterRadius/2); // angle spanned by a single strip
    //const Double_t RingPitch = (OuterRadius-InnerRadius)/16;


    Int_t Side,WireID;
    Double_t PCDownstream[NPCWires],PCUpstream[NPCWires];
    Double_t PCDownVoltage[NPCWires],PCUpVoltage[NPCWires];
    Int_t ConvTest=0;
    Double_t Vcal;

    Int_t status = 0;
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
      MyFill("ADCNhits",101,0,100,ADC.Nhits);

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
      for ( Int_t i=0; i<NPCWires; i++ ){
	if ( PCDownstream[i]>100 && PCUpstream[i]>100 && PCDownstream[i]<4000 && PCUpstream[i]<4000 ){
	  PC.WireID[PC.NPCHits] = i;
	  PC.Down[PC.NPCHits] = PCDownstream[i];
	  PC.Up[PC.NPCHits] = PCUpstream[i];
	  PC.DownVoltage[PC.NPCHits] = PCDownVoltage[i];
	  PC.UpVoltage[PC.NPCHits] = PCUpVoltage[i];
	  if (PCDownVoltage[i]>0 && PCUpVoltage[i]>0){
	    PC.Energy[PC.NPCHits] = PC.DownVoltage[PC.NPCHits] + PC.UpVoltage[PC.NPCHits];
	    PC.Zp[PC.NPCHits] = (PC.UpVoltage[PC.NPCHits] - PC.DownVoltage[PC.NPCHits])/PC.Energy[PC.NPCHits];
	  }

	  CMAP->GetPCWorldCoordinates(PC.WireID[PC.NPCHits],PC.Zp[PC.NPCHits],XWPC,YWPC,ZWPC,RWPC,PhiWPC);
	  PC.XWp[PC.NPCHits] = XWPC;
	  PC.YWp[PC.NPCHits] = YWPC;
	  PC.ZWp[PC.NPCHits] = ZWPC;
	  PC.RWp[PC.NPCHits] = RWPC;
	  PC.PhiWp[PC.NPCHits] = PhiWPC;

	  PC.NPCHits++;
	}
      }

      if (TDC.Nhits>MaxTDCHits) TDC.Nhits=MaxTDCHits;
  
      RFTime = 0;
      MCPTime = 0;
  
  
      for (Int_t n=0; n<TDC.Nhits; n++){
	//TDC ID == 2?
	if( TDC.ID[n]==2 ){
	  if (TDC.ChNum[n]==0)  RFTime   = TDC.Data[n];
	  if (TDC.ChNum[n]==7)  MCPTime  = TDC.Data[n];
	}
      }

#ifdef TimingCut   
      if ( (MCPTime-RFTime)%538<71 || (MCPTime-RFTime)%538>370 ){
	continue;
      }
      if ( (MCPTime-RFTime)%538>100 && (MCPTime-RFTime)%538<330 ){
	continue;
      }
#endif   
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

	if(DN>27) {continue;}
	if(DN>3){
	  CMAP->GetX3MeVPerChannel1(DN,DetCh,slope_rel);  
	  CMAP->GetX3MeVPerChannel2(DN,DetCh,slope_alpha);  
	  CMAP->GetX3FinalEnergyOffsetInMeV(DN,DetCh,FinalShift);
	  SiEnergy[DN-4][DetCh] = (Double_t)Si_Old.Energy[n];
	  SiEnergy_Pulser[DN-4][DetCh] = (Double_t)Si_Old.Energy[n]+ChRandom->Rndm()-0.5+ZeroShift/VperCh;
	  SiEnergy_Rel[DN-4][DetCh] = SiEnergy_Pulser[DN-4][DetCh]*slope_rel;
	  SiEnergy_Cal[DN-4][DetCh] = SiEnergy_Rel[DN-4][DetCh]*slope_alpha;
	  SiEnergy_Cal[DN-4][DetCh] += FinalShift;
	  SiTime[DN-4][DetCh]   = (Double_t)Si_Old.Time[n];
	}else{
	  CMAP->GetQQQ3MeVPerChannel1(DN,DetCh,slope_rel);
	  CMAP->GetQQQ3MeVPerChannel2(DN,DetCh,slope_alpha);
	  CMAP->GetQQQ3FinalEnergyOffsetInMeV(DN,DetCh,FinalShift);
	  QQQ3Energy[DN][DetCh] = (Double_t)Si_Old.Energy[n];
	  QQQ3Energy_Pulser[DN][DetCh] = (Double_t)Si_Old.Energy[n]+ChRandom->Rndm()-0.5+ZeroShift/VperCh;
	  QQQ3Energy_Rel[DN][DetCh] = QQQ3Energy_Pulser[DN][DetCh]*slope_rel;
	  QQQ3Energy_Cal[DN][DetCh] = QQQ3Energy_Rel[DN][DetCh]*slope_alpha;
	  QQQ3Energy_Cal[DN][DetCh] += FinalShift;
	  QQQ3Time[DN][DetCh] = (Double_t)Si_Old.Time[n];
	}
      }// End loop over X3_Nhits------------------------------------------------------------------------------------------------------------------------------------------------------------------

      for (int i=0; i<NumX3; i++){//Determine multiplicities for SX3 detector: forward and back--------------------------------------------------------------------------------------------------
	Si.zeroPlaceHolder();
	for (int j=0; j<MaxX3Ch; j++){
	  if ( (SiEnergy[i][j] > 0.0)){
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
	//if ( Si.DownMult[Si.NSiHits]!=0 || Si.UpMult[Si.NSiHits]!=0 || Si.BackMult[Si.NSiHits]!=0 ){
	if ( Si.det_place_holder.EnergyDown_Cal.size()!=0 || Si.det_place_holder.EnergyUp_Cal.size()!=0 || Si.det_place_holder.EnergyBack_Cal.size()!=0 ){
	  //cout << "Sorted Data\n";
	  Si.det_place_holder.UpMult = Si.det_place_holder.EnergyUp_Cal.size();
	  Si.det_place_holder.DownMult = Si.det_place_holder.EnergyDown_Cal.size();
	  Si.det_place_holder.BackMult = Si.det_place_holder.EnergyBack_Cal.size();
	  Si.det_place_holder.DetID = i+4;
	  //cout << Si.det_place_holder.DetID << endl;
	  Si.det_place_holder.HitType = Si.det_place_holder.BackMult*100 + Si.det_place_holder.UpMult*10 + Si.det_place_holder.DownMult;

	  Si.hit_place_holder.NHitsInDet = 1;
	  Si.hit_place_holder.DetID = i+4;
	  Si.hit_place_holder.HitType = Si.det_place_holder.BackMult*100 + Si.det_place_holder.UpMult*10 + Si.det_place_holder.DownMult;
	  if ( Si.det_place_holder.HitType == 111 ){
	    if ( Si.det_place_holder.DownChNum[0]==Si.det_place_holder.UpChNum[0] ){
	      Si.hit_place_holder.EnergyBack = Si.det_place_holder.EnergyBack_Cal[0];
	      Si.hit_place_holder.EnergyFront = Si.det_place_holder.EnergyUp_Cal[0] + Si.det_place_holder.EnergyDown_Cal[0];
	      Si.hit_place_holder.Energy = Si.det_place_holder.EnergyBack_Cal[0];
 	      Si.hit_place_holder.Time = Si.det_place_holder.TimeBack[0];
 	      Si.hit_place_holder.RFSubtract = Si.hit_place_holder.Time - RFTime;
	      Double_t FinalZUp=0,FinalZUpCal=-10, FinalZDown=0, FinalZDownCal=-10;
	      FinalZDown = 2*Si.det_place_holder.EnergyDown_Cal[0]/Si.det_place_holder.EnergyBack_Cal[0]-1;
	      CMAP->PosCal(Si.hit_place_holder.DetID,Si.det_place_holder.DownChNum[0],Si.det_place_holder.BackChNum[0],FinalZDown,FinalZDownCal);
	      FinalZUp = 1-2*Si.det_place_holder.EnergyUp_Cal[0]/Si.det_place_holder.EnergyBack_Cal[0];
	      CMAP->PosCal(Si.hit_place_holder.DetID,Si.det_place_holder.UpChNum[0],Si.det_place_holder.BackChNum[0],FinalZUp,FinalZUpCal);
	      if (FinalZDownCal > -1 && FinalZUpCal > -1){
		Si.hit_place_holder.Z = FinalZUpCal;
		Si.hit_place_holder.X = XRandom->Rndm()+(3-Si.det_place_holder.UpChNum[0]);
	      }else{
		Si.hit_place_holder.HitType += 1000;
	      }
	      Double_t xw=0, yw=0, zw=0, rw=0, phiw=0;
	      CMAP->GetX3WorldCoordinates(Si.hit_place_holder.DetID,Si.hit_place_holder.X,Si.hit_place_holder.Z,xw,yw,zw,rw,phiw);
	      Si.hit_place_holder.XW = xw;
	      Si.hit_place_holder.YW = yw;
	      Si.hit_place_holder.ZW = zw;
	      Si.hit_place_holder.RW= rw;
	      Si.hit_place_holder.PhiW = phiw;
	    }else{
	      Si.hit_place_holder.HitType += 1000;
	    }
	  }else if ( (Si.det_place_holder.HitType == 110) || (Si.det_place_holder.HitType == 101) ){
	    Si.hit_place_holder.EnergyBack = (Si.det_place_holder.EnergyBack_Cal[0]);
	    Si.hit_place_holder.Energy = (Si.det_place_holder.EnergyBack_Cal[0]);
	    Si.hit_place_holder.Time = (Si.det_place_holder.TimeBack[0]);
	    Si.hit_place_holder.RFSubtract = (Si.hit_place_holder.Time - RFTime);
	    Double_t FinalZUp=0,FinalZUpCal=-10, FinalZDown=0, FinalZDownCal=-10;
	    if ( Si.det_place_holder.HitType == 110 ){
	      Si.hit_place_holder.EnergyFront = (Si.det_place_holder.EnergyUp_Cal[0]);
	      FinalZUp = 1-2*Si.det_place_holder.EnergyUp_Cal[0]/Si.det_place_holder.EnergyBack_Cal[0];
	      CMAP->PosCal(Si.hit_place_holder.DetID,Si.det_place_holder.UpChNum[0],Si.det_place_holder.BackChNum[0],FinalZUp,FinalZUpCal);
	      if (FinalZUpCal > -1){
		Si.hit_place_holder.Z = (FinalZUpCal);
		Si.hit_place_holder.X = (XRandom->Rndm()+(3-Si.det_place_holder.UpChNum[0]));
	      }else{
		Si.hit_place_holder.HitType += 1000;
	      }
	    }else{
	      Si.hit_place_holder.EnergyFront = (Si.det_place_holder.EnergyDown_Cal[0]);
	      FinalZDown = 2*Si.det_place_holder.EnergyDown_Cal[0]/Si.det_place_holder.EnergyBack_Cal[0]-1;
	      CMAP->PosCal(Si.hit_place_holder.DetID,Si.det_place_holder.DownChNum[0],Si.det_place_holder.BackChNum[0],FinalZDown,FinalZDownCal);
	      if (FinalZDownCal > -1){
		Si.hit_place_holder.Z = (FinalZDownCal);
		Si.hit_place_holder.X = (XRandom->Rndm()+(3-Si.det_place_holder.DownChNum[0]));
	      }else{
		Si.hit_place_holder.HitType += 1000;
	      }
	    }
	    Double_t xw=0, yw=0, zw=0, rw=0, phiw=0;
	    CMAP->GetX3WorldCoordinates(Si.hit_place_holder.DetID,Si.hit_place_holder.X,Si.hit_place_holder.Z,xw,yw,zw,rw,phiw);
	    Si.hit_place_holder.XW = (xw);
	    Si.hit_place_holder.YW = (yw);
	    Si.hit_place_holder.ZW = (zw);
	    Si.hit_place_holder.RW = (rw);
	    Si.hit_place_holder.PhiW = (phiw);
	  }else if ( Si.det_place_holder.HitType == 211 ){
	    if ( Si.det_place_holder.DownChNum[0]==Si.det_place_holder.UpChNum[0] ){
	      if ( fabs(Si.det_place_holder.BackChNum[1]-Si.det_place_holder.BackChNum[0]) == 1 ){
		Si.hit_place_holder.EnergyBack = (Si.det_place_holder.EnergyBack_Cal[0] + Si.det_place_holder.EnergyBack_Cal[1] );
		Si.hit_place_holder.EnergyFront = (Si.det_place_holder.EnergyUp_Cal[0] + Si.det_place_holder.EnergyDown_Cal[0]);
		Si.hit_place_holder.Energy = (Si.det_place_holder.EnergyBack_Cal[0] + Si.det_place_holder.EnergyBack_Cal[1]);
		Si.hit_place_holder.Time = (Si.det_place_holder.TimeBack[0]);
		Si.hit_place_holder.RFSubtract = (Si.hit_place_holder.Time - RFTime);

		Si.hit_place_holder.Z = (1.875*(3-Si.det_place_holder.BackChNum[0]) - 0.1 + 0.2*ZRandom->Rndm());//assuming a mm inbetween strips
		Si.hit_place_holder.X = (XRandom->Rndm()+(3-Si.det_place_holder.DownChNum[0]));

		Double_t xw=0, yw=0, zw=0, rw=0, phiw=0;
		CMAP->GetX3WorldCoordinates(Si.hit_place_holder.DetID,Si.hit_place_holder.X,Si.hit_place_holder.Z,xw,yw,zw,rw,phiw);
		Si.hit_place_holder.XW = (xw);
		Si.hit_place_holder.YW = (yw);
		Si.hit_place_holder.ZW = (zw);
		Si.hit_place_holder.RW = (rw);
		Si.hit_place_holder.PhiW = (phiw);
	      }else{
		Si.hit_place_holder.HitType += 1000;
	      }
	    }else{
	      Si.hit_place_holder.HitType += 1000;
	    }
	  }else{
	    Si.hit_place_holder.HitType += 1000;
	  }

	  Si.Detector.push_back(Si.det_place_holder);
	  Si.Hit.push_back(Si.hit_place_holder);

	  Si.NSiHits++;
	}
      }

      for (int i=0; i<NumQQQ3; i++){//Determine multiplicities for QQQ detector: forward and back-------------------------------------------------------------------------------------------------------
	Si.zeroPlaceHolder();
	for (int j=0; j<MaxQQQ3Ch; j++){
	  if ( (QQQ3Energy[i][j] > 0)){
	    if (j<16){
	      Si.det_place_holder.BackChNum.push_back(j);
	      Si.det_place_holder.EnergyBack_Raw.push_back(QQQ3Energy[i][j]);
	      Si.det_place_holder.EnergyBack_Pulser.push_back(QQQ3Energy_Pulser[i][j]);
	      Si.det_place_holder.EnergyBack_Rel.push_back(QQQ3Energy_Rel[i][j]);
	      Si.det_place_holder.EnergyBack_Cal.push_back(QQQ3Energy_Cal[i][j]);
	      Si.det_place_holder.TimeBack.push_back(QQQ3Time[i][j]);
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
	  //cout << "Sorted Data\n";

	  ProcessQQQ(&Si,CMAP,RFTime);
	  Si.Detector.push_back(Si.det_place_holder);
	  Si.Hit.push_back(Si.hit_place_holder);
	  Si.NSiHits++;

	}else if ( Si.det_place_holder.EnergyUp_Cal.size()!=0 || Si.det_place_holder.EnergyBack_Cal.size()!=0 ){
	  Si.Detector.push_back(Si.det_place_holder);
	  Si.Hit.push_back(Si.hit_place_holder);
	  Si.NSiHits++;
	}
      }

      if (Si.NSiHits != Si.Detector.size()){
	cout << Si.NSiHits << "  " << Si.Detector.size() << endl;
      }
      if (Si.NSiHits != Si.Hit.size()){
	cout << Si.NSiHits << "  " << Si.Detector.size() << endl;
      }
      //------------------------------------------------------------------------------------------------------------------------------
      //if (global_evt>10){
	//exit(EXIT_FAILURE);
      //}



      if (Si.NSiHits>0 ){    
	MainTree->Fill();
      }
    }
    cout << "Done with File: " << rootfile << endl;
    //delete inputFile;

  }

  outputFile->cd();
  cout << "change to output file directory" << endl;
  //fhlist->Write();
  RootObjects->Write();
  cout << "Write Root Objects" << endl;
  outputFile->Close();
  cout << "Close Outputfile\n";
  return 0;
}


void ProcessQQQ(SiHit *Si,ChannelMap *CMAP, Int_t RFTime){
  Double_t QQQR;
  Double_t QQQPhi;

  const Double_t InnerRadius = 5.01; //cm4.81
  const Double_t OuterRadius = 10.1; //cm
  const Double_t ArcLength = 0.96; //cm (distance between strips at outer radius) 
  const Double_t StripAngle = 2*TMath::ASin(ArcLength/OuterRadius/2); // angle spanned by a single strip
  const Double_t RingPitch = (OuterRadius-InnerRadius)/16;

  TRandom3 *PhiRandom = new TRandom3();
  TRandom3 *RRandom = new TRandom3();

  PhiRandom->SetSeed();
  RRandom->SetSeed();

  //cout << Si->det_place_holder.DetID << endl;

  data data_place_holder;

  vector<data> front;
  vector<data> back;
  vector<data> front_compress;
  vector<data> back_compress;

  front.clear();
  back.clear();
  front_compress.clear();
  back_compress.clear();

  //copy contents of det_place_holder into new arrays, get rid of bad energies
  for ( Int_t i=0; i<Si->det_place_holder.UpMult; i++ ){
    data_place_holder.Channel = (Double_t)Si->det_place_holder.UpChNum[i];
    data_place_holder.Energy = Si->det_place_holder.EnergyUp_Cal[i];
    data_place_holder.Time = Si->det_place_holder.TimeUp[i];
    if (data_place_holder.Energy>0){//if the energy is less than or equal to zero, we want to discard it
      front.push_back(data_place_holder);
    }
  }
  for ( Int_t i=0; i<Si->det_place_holder.BackMult; i++ ){
    data_place_holder.Channel = (Double_t)Si->det_place_holder.BackChNum[i];
    data_place_holder.Energy = Si->det_place_holder.EnergyBack_Cal[i];
    data_place_holder.Time = Si->det_place_holder.TimeBack[i];
    if (data_place_holder.Energy>0){//if the energy is less than or equal to zero, we want to discard it
      back.push_back(data_place_holder);
    }
  }
  if (front.size()==0 || back.size()==0){
    return;
  }
  //-----------------------------------------------------------------
  //sort channels from lowest to hightest
  sort( front.begin(), front.end(), Csort_method );
  sort( back.begin(), back.end(), Csort_method );
  //-------------------------------------------------------------------
  //compress adjacent channels (add back energies)
  //front
  data_place_holder.Energy = front.at(0).Energy;
  data_place_holder.Time = front.at(0).Time;
  data_place_holder.Channel = front.at(0).Channel;
  for ( Int_t i=1; i<front.size(); i++ ){
    if ( (front.at(i).Channel-front.at(i-1).Channel) == 1 ){//if front channels are adjacent
      data_place_holder.Channel = (front.at(i).Channel*front.at(i).Energy + data_place_holder.Channel*data_place_holder.Energy);
      data_place_holder.Energy = front.at(i).Energy + data_place_holder.Energy;
      data_place_holder.Channel = data_place_holder.Channel/data_place_holder.Energy;
      data_place_holder.Time = (front.at(i).Time + data_place_holder.Time)/2;
    }else{
      front_compress.push_back(data_place_holder);
      data_place_holder.Energy = front.at(i).Energy;
      data_place_holder.Time = front.at(i).Time;
      data_place_holder.Channel = front.at(i).Channel;
    }
  }
  front_compress.push_back(data_place_holder);
  //back
  data_place_holder.Energy = back.at(0).Energy;
  data_place_holder.Time = back.at(0).Time;
  data_place_holder.Channel = back.at(0).Channel;
  for ( Int_t i=1; i<back.size(); i++ ){
    if ( (back.at(i).Channel-back.at(i-1).Channel) == 1 ){//if back channels are adjacent
      data_place_holder.Channel = (back.at(i).Channel*back.at(i).Energy + data_place_holder.Channel*data_place_holder.Energy);
      data_place_holder.Energy = back.at(i).Energy + data_place_holder.Energy;
      data_place_holder.Channel = data_place_holder.Channel/data_place_holder.Energy;
      data_place_holder.Time = (back.at(i).Time + data_place_holder.Time)/2;
    }else{
      back_compress.push_back(data_place_holder);
      data_place_holder.Energy = back.at(i).Energy;
      data_place_holder.Time = back.at(i).Time;
      data_place_holder.Channel = back.at(i).Channel;
    }
  }
  back_compress.push_back(data_place_holder);
  //---------------------------------------------------------------------------------
  //sort energies, largest to smallest
  sort(front_compress.begin(), front_compress.end(), Esort_method );
  sort(back_compress.begin(), back_compress.end(), Esort_method );
  //------------------------------------------------------------------------------------
  //determine total energy in front and back
  //determine front and back channels
  //for each channel fired, weight for the energy deposited in the channel
  Double_t front_sum = 0;
  Double_t back_sum = 0;
  Double_t avg_time = 0;//just want a time for now
  Double_t avg_front_ch = 0;
  Double_t avg_back_ch = 0;
  for (Int_t i=0; i<front_compress.size(); i++){
    avg_front_ch += front_compress.at(i).Channel*front_compress.at(i).Energy;
    front_sum += front_compress.at(i).Energy;
  }
  for (Int_t i=0; i<back_compress.size(); i++){
    avg_back_ch += back_compress.at(i).Channel*back_compress.at(i).Energy;
    back_sum += back_compress.at(i).Energy;
    avg_time = back_compress.at(i).Time;
  }
  if ( front_sum<=0 || back_sum<=0 ){
    return;//if either the front or back energies don't make since, just return
  }
  avg_front_ch = avg_front_ch/front_sum;
  avg_back_ch = avg_back_ch/back_sum;

  if ( fabs(1-front_sum/back_sum)>0.1 ){//if the front energy doesn't match the back energy, do stuff
    Int_t found_good = 0;
    //first, compare each front hit to each back hit.
    //if an individual hit works, take it and move on
    for (Int_t i=0; i<front_compress.size(); i++){
      for (Int_t j=0; j<back_compress.size(); j++){
	if ( fabs(1-front_compress.at(i).Energy/back_compress.at(j).Energy)<0.1 ){
	  found_good++;
	  front_sum = front_compress.at(i).Energy;
	  back_sum = back_compress.at(j).Energy;
	  avg_front_ch = front_compress.at(i).Channel;
	  avg_back_ch = back_compress.at(j).Channel;
	  avg_time = back_compress.at(j).Time;
	  Si->hit_place_holder.Flag = 1;
	  break;
	}
      }
      if (found_good==1){
	break;
      }
    }
    if (found_good==0){//if no matches were found above, try something else
      //subtract each front energy from the total and see if that matches the total back
      Double_t new_front_sum = 0;
      Double_t new_front_ch = 0;
      for (Int_t i=0; i<front_compress.size(); i++){
	new_front_sum = front_sum - front_compress.at(i).Energy;
	new_front_ch = (avg_front_ch*front_sum - front_compress.at(i).Channel*front_compress.at(i).Energy)/new_front_sum;
	if ( fabs(1-new_front_sum/back_sum)<0.1 ){
	  found_good++;
	  front_sum = new_front_sum;
	  avg_front_ch = new_front_ch;
	  Si->hit_place_holder.Flag = 2;
	  break;
	}
      }

      if (found_good==0){
	//if subtracting each front channel didn't work, now try the back
	Double_t new_back_sum = 0;
	Double_t new_back_ch = 0;
	for (Int_t i=0; i<back_compress.size(); i++){
	  new_back_sum = back_sum - back_compress.at(i).Energy;
	  new_back_ch = (avg_back_ch*back_sum - back_compress.at(i).Channel*back_compress.at(i).Energy)/new_back_sum;
	  if ( fabs(1-new_back_sum/front_sum)<0.1 ){
	    found_good++;
	    back_sum = new_back_sum;
	    avg_back_ch = new_back_ch;
	    Si->hit_place_holder.Flag = 2;
	    break;
	  }
	}
      }

    }

  }

  //-----------------------------------------------------------------------------------------
  //determine radius and angle
  QQQR = OuterRadius - (avg_front_ch+RRandom->Rndm() )*RingPitch;
  QQQPhi = (avg_back_ch+PhiRandom->Rndm() )*StripAngle;
  //-------------------------------------------------------------
  //fill tree
  Si->hit_place_holder.NHitsInDet = 1;
  Si->hit_place_holder.DetID = Si->det_place_holder.DetID;
  Si->hit_place_holder.HitType = (back_compress.size()*10 + front_compress.size());
  Si->hit_place_holder.FrontChannel = avg_front_ch;
  Si->hit_place_holder.BackChannel = avg_back_ch;
  Si->hit_place_holder.EnergyFront = (front_sum);
  Si->hit_place_holder.EnergyBack = (back_sum);
  Si->hit_place_holder.Energy = (back_sum);
  //cout << back_compress.at(0).Time << endl;
  Si->hit_place_holder.Time = (avg_time);
  Si->hit_place_holder.RFSubtract = (Si->hit_place_holder.Time - RFTime);
  Si->hit_place_holder.X = (QQQR*TMath::Cos(QQQPhi));
  Si->hit_place_holder.Y = (QQQR*TMath::Sin(QQQPhi));

  Double_t xw=0, yw=0, rw=0, phiw=0;
  CMAP->GetQQQ3WorldCoordinates(Si->hit_place_holder.DetID,Si->hit_place_holder.X,Si->hit_place_holder.Y,xw,yw,rw,phiw);
  Si->hit_place_holder.XW = (xw);
  Si->hit_place_holder.YW = (yw);
  Si->hit_place_holder.RW = (rw);
  Si->hit_place_holder.PhiW = (phiw);
  
  if ( Si->hit_place_holder.DetID==1 || Si->hit_place_holder.DetID==2 ){
    Si->hit_place_holder.Z = (0.762);
    Si->hit_place_holder.ZW = (0.762);
  }else if (Si->hit_place_holder.DetID==0 || Si->hit_place_holder.DetID==3){
    Si->hit_place_holder.Z = (0);
    Si->hit_place_holder.ZW = (0);
  }else{
    //cout << QQQ3DetNum << endl;
    //cout << "This shouldn't ever happen\n";
  }

  delete PhiRandom;
  delete RRandom;

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

