//////////////////////////////////////////////////////////////////////////////////////
#include <TROOT.h>
#include <vector>
using namespace std;
//////////////////////////////////////////////////////////////////////////////////////
#define MaxPCHits 24
#define NPCWires  24
#define MaxSiHits 500
#define MaxTracks 100
#define MaxCsIHits 500
//////////////////////////////////////////////////////////////////////////////////////
class SiHit {  // This class is for Silicon hits

 public:
  SiHit(){};

  Int_t NSiHits;
//-------------------------------------------------------------
  struct SortByDetector {
    Int_t DetID;

    Int_t UpMult;    
    Int_t DownMult;
    Int_t FrontMult;
    Int_t BackMult;

    Int_t HitType;

    vector<Int_t> UpChNum;
    vector<Int_t> DownChNum;
    vector<Int_t> FrontChNum;
    vector<Int_t> BackChNum;

    //Raw Energy --in channels
    vector<Double_t> EUp_Raw;
    vector<Double_t> EDown_Raw;   
    vector<Double_t> EFront_Raw;
    vector<Double_t> EBack_Raw;

    //Energy after pulser calibration --in channels
    vector<Double_t> EUp_Pulser;
    vector<Double_t> EDown_Pulser;  
    vector<Double_t> EFront_Pulser;
    vector<Double_t> EBack_Pulser;

    //Energy after relative calibration --in channels
    vector<Double_t> EUp_Rel;
    vector<Double_t> EDown_Rel;    
    vector<Double_t> EFront_Rel;
    vector<Double_t> EBack_Rel;  

    //Energy after Energy calibration --in MeV
    vector<Double_t> EUp_Cal;
    vector<Double_t> EDown_Cal;   
    vector<Double_t> EFront_Cal;
    vector<Double_t> EBack_Cal;

    vector<Double_t> TUp;
    vector<Double_t> TDown;    
    vector<Double_t> TFront;
    vector<Double_t> TBack;

  }det_obj;
  //-------------------------------------------------------------
  struct SortByHit {
    Int_t NHitsInDet;
    Int_t DetID;
    
    Int_t HitType;
    Int_t FrontChannel, BackChannel;
    Double_t EnergyBack, EnergyFront;
    Double_t Energy, Time;
    Double_t X, Y, Z;
    Double_t ZUp_Dummy,ZDown_Dummy;
    Double_t XW, YW, ZW;
    Double_t RW, PhiW;

  }hit_obj;
//-------------------------------------------------------------
  vector<SortByDetector> Detector;
  vector<SortByDetector> *ReadDet;
  vector<SortByHit> Hit;  
  vector<SortByHit> *ReadHit;
  //vector<SortByHit>::iterator it;
  //-------------------------------------------------------------
  void ZeroSi_obj(){
    det_obj.DetID = -1;
    det_obj.UpMult = 0;
    det_obj.DownMult = 0;
    det_obj.BackMult = 0;
    det_obj.FrontMult = 0;
    det_obj.HitType = 0;

    det_obj.UpChNum.clear();
    det_obj.DownChNum.clear();
    det_obj.FrontChNum.clear();
    det_obj.BackChNum.clear();

    det_obj.EUp_Raw.clear();
    det_obj.EDown_Raw.clear();   
    det_obj.EFront_Raw.clear();
    det_obj.EBack_Raw.clear();

    det_obj.EUp_Pulser.clear();
    det_obj.EDown_Pulser.clear();    
    det_obj.EFront_Pulser.clear();
    det_obj.EBack_Pulser.clear();

    det_obj.EUp_Rel.clear();
    det_obj.EDown_Rel.clear();  
    det_obj.EFront_Rel.clear();
    det_obj.EBack_Rel.clear();   

    det_obj.EUp_Cal.clear();
    det_obj.EDown_Cal.clear();    
    det_obj.EFront_Cal.clear();
    det_obj.EBack_Cal.clear();

    det_obj.TUp.clear();
    det_obj.TDown.clear();    
    det_obj.TFront.clear();
    det_obj.TBack.clear();

    hit_obj.NHitsInDet = 0;
    hit_obj.DetID = -1;
    
    hit_obj.HitType = 0;
    hit_obj.FrontChannel = -1;
    hit_obj.BackChannel = -1;
    hit_obj.EnergyBack = 0;
    hit_obj.EnergyFront = 0;
    hit_obj.Energy = 0;
    hit_obj.Time = 0;
    hit_obj.X = 0;
    hit_obj.Y = 0;
    hit_obj.Z = sqrt(-1);
    hit_obj.ZUp_Dummy = sqrt(-1);
    hit_obj.ZDown_Dummy = sqrt(-1);
    hit_obj.XW = 0;
    hit_obj.YW = 0;
    hit_obj.ZW = sqrt(-1);
    hit_obj.RW = 0;
    hit_obj.PhiW = sqrt(-1);
    //hit_obj.RFSubtract = 0;
  };
//-------------------------------------------------------------
  void zeroSiHit(){
    NSiHits = 0;
    Detector.clear();
    Hit.clear();  

    ZeroSi_obj();
  };
  //-------------------------------------------------------------
};
//////////////////////////////////////////////////////////////////////////////////////
class PCHit { // This class is for Proportional Counter hits
 
 public:
  PCHit(){};

  Int_t NPCHits;
  //-------------------------------------------------------------
  struct SortByPC {
    Int_t WireID;
    Double_t Down, Up;
    Double_t DownVoltage, UpVoltage;
    Double_t Energy;
    Double_t Z;
    Double_t XW, YW, ZW;
    Double_t RW, PhiW;
    Int_t TrackType;

  }pc_obj;
  //-------------------------------------------------------------
  vector<SortByPC> Hit;
  vector<SortByPC> *ReadHit;
  //-------------------------------------------------------------
  void ZeroPC_obj(){
    pc_obj.WireID = -1;
    pc_obj.Down = 0;
    pc_obj.Up = 0;
    pc_obj.DownVoltage = 0;
    pc_obj.UpVoltage = 0;
    pc_obj.Energy = sqrt(-1);
    pc_obj.Z = sqrt(-1);
    pc_obj.XW = 0;
    pc_obj.YW = 0;
    pc_obj.ZW = 0;
    pc_obj.RW = 0;
    pc_obj.PhiW = sqrt(-1);
    pc_obj.TrackType = 0;
  };
  //-------------------------------------------------------------
  void zeroPCHit(){
    NPCHits = 0;
    Hit.clear();
    ZeroPC_obj();
  };
  //-------------------------------------------------------------
};
/////////////////////////////////////////////////////////////////////////////////////
class Track { //This Class is for Tracking

 public:
  Track(){};

  Int_t NTracks;
  Int_t NTracks1;
  Int_t NTracks2;
  Int_t NTracks3;
  //------------------------------------------------------
  struct TrackEvent {

    Int_t TrackType;
    Int_t HitType;
    
    Double_t SiEnergy;
    Double_t SiZ;
    Double_t SiX;
    Double_t SiY;
    Double_t SiR;
    Double_t SiPhi;
    Int_t DetID;    

    Double_t PCEnergy;    
    Double_t PCZ;
    Double_t PCX; 
    Double_t PCY;
    Double_t PCR;
    Double_t PCPhi;
    Int_t WireID;

    Double_t IntPoint;
    Double_t IntPoint_X;
    Double_t IntPoint_Y;

    Double_t DiffIntPoint;
    Double_t DiffIntPoint_X;
    Double_t DiffIntPoint_Y;
  
    Double_t PathLength;

    Double_t Theta;
    Double_t Phi;

    Double_t BeamEnergy;
    Double_t EnergyLoss;

    Double_t PCZ_Ref;
  }track_obj;
  //-------------------------------------------------------------
  vector<TrackEvent> TrEvent;
  vector<TrackEvent> *ReadTrEvent; 
  //-------------------------------------------------------------
  void ZeroTr_obj(){
    track_obj.TrackType = 0;
    track_obj.HitType = 0;
    
    track_obj.SiEnergy = 0;
    track_obj.SiZ = -100.0;
    track_obj.SiX = 0;
    track_obj.SiY = 0;
    track_obj.SiR = 0;
    track_obj.SiPhi = sqrt(-1);
    track_obj.DetID = -1;

    track_obj.PCEnergy = sqrt(-1);
    track_obj.PCX = 0;
    track_obj.PCZ = sqrt(-1);
    track_obj.PCY = 0;
    track_obj.PCR = 0;
    track_obj.PCPhi = sqrt(-1);
    track_obj.WireID = -1;

    track_obj.IntPoint = 0;
    track_obj.IntPoint_X = 0;
    track_obj.IntPoint_Y = 0;
    track_obj.PathLength = 0;

    track_obj.DiffIntPoint = 0;
    track_obj.DiffIntPoint_X = 0;
    track_obj.DiffIntPoint_Y = 0;

    track_obj.Theta = sqrt(-1);
    track_obj.Phi = sqrt(-1);

    track_obj.BeamEnergy = 0; 
    track_obj.EnergyLoss=0;
  };
  //-------------------------------------------------------------
  void zeroTrack(){

    NTracks = 0;
    NTracks1 = 0;
    NTracks2 = 0;
    NTracks3 = 0;

    TrEvent.clear();
    ZeroTr_obj();
  };
  //-------------------------------------------------------------
  //static bool Si_sort_method(struct Silicon_Event a,struct Silicon_Event b);
  //static bool PC_sort_method(struct PropCounter_Event c,struct PropCounter_Event d);
  static bool Tr_Sisort_method(struct TrackEvent a,struct TrackEvent b);
  static bool Tr_PCsort_method(struct TrackEvent a,struct TrackEvent b);
  //-------------------------------------------------------------
};
//////////////////////////////////////////////////////////////////////////////////////////
