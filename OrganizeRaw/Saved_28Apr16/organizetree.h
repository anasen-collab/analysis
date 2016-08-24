//////////////////////////////////////////////////////////////////////////////////////
#include <TROOT.h>
#include <vector>
using namespace std;

// 12/09/15 Modified for 24 wires in PC
//2016/01/04 Modified for different X, Y, Z for Silicon & PC
//////////////////////////////////////////////////////////////////////////////////////
#define MaxPCHits 24
#define NPCWires  24
#define MaxSiHits 500
#define MaxTracks 100
//////////////////////////////////////////////////////////////////////////////////////
//gROOT->ProcessLine("#include <vector>");

class SiHit {
  // This class is for Si hit
 public:
  SiHit(){};

  Int_t NSiHits;

  struct SortByDetector {
    Int_t DetID;
    Int_t UpMult;
    Int_t DownMult;
    Int_t BackMult;
    Int_t HitType;

    vector<Int_t> UpChNum;
    vector<Int_t> DownChNum;
    vector<Int_t> BackChNum;
    vector<Double_t> EnergyUp_Raw;
    vector<Double_t> EnergyDown_Raw;
    vector<Double_t> EnergyBack_Raw;
    vector<Double_t> EnergyUp_Pulser;
    vector<Double_t> EnergyDown_Pulser;
    vector<Double_t> EnergyBack_Pulser;
    vector<Double_t> EnergyUp_Rel;
    vector<Double_t> EnergyDown_Rel;
    vector<Double_t> EnergyBack_Rel;
    vector<Double_t> EnergyUp_Cal;
    vector<Double_t> EnergyDown_Cal;
    vector<Double_t> EnergyBack_Cal;

    vector<Double_t> TimeUp;
    vector<Double_t> TimeDown;
    vector<Double_t> TimeBack;

  }det_place_holder;

  struct SortByHit {
    Int_t NHitsInDet;
    Int_t DetID;

    Int_t Flag;
    Int_t HitType;
    Double_t FrontChannel, BackChannel;
    Double_t EnergyBack, EnergyFront;
    Double_t Energy, Time;
    Double_t X, Y, Z;
    Double_t XW, YW, ZW;
    Double_t RW, PhiW;
    Double_t RFSubtract;
    Int_t TrackType;

  }hit_place_holder;

  vector<SortByDetector> Detector;
  vector<SortByHit> Hit;

  vector<SortByDetector> *ReadDet;
  vector<SortByHit> *ReadHit;

  //vector<SortByHit>::iterator it;

  void zeroPlaceHolder(){
    det_place_holder.DetID = -1;
    det_place_holder.UpMult = 0;
    det_place_holder.DownMult = 0;
    det_place_holder.BackMult = 0;
    det_place_holder.HitType = 0;

    det_place_holder.UpChNum.clear();
    det_place_holder.DownChNum.clear();
    det_place_holder.BackChNum.clear();
    det_place_holder.EnergyUp_Raw.clear();
    det_place_holder.EnergyDown_Raw.clear();
    det_place_holder.EnergyBack_Raw.clear();
    det_place_holder.EnergyUp_Pulser.clear();
    det_place_holder.EnergyDown_Pulser.clear();
    det_place_holder.EnergyBack_Pulser.clear();
    det_place_holder.EnergyUp_Rel.clear();
    det_place_holder.EnergyDown_Rel.clear();
    det_place_holder.EnergyBack_Rel.clear();
    det_place_holder.EnergyUp_Cal.clear();
    det_place_holder.EnergyDown_Cal.clear();
    det_place_holder.EnergyBack_Cal.clear();

    det_place_holder.TimeUp.clear();
    det_place_holder.TimeDown.clear();
    det_place_holder.TimeBack.clear();

    hit_place_holder.NHitsInDet = 0;
    hit_place_holder.DetID = -1;
    hit_place_holder.Flag = 0;
    hit_place_holder.HitType = 0;
    hit_place_holder.FrontChannel = -1;
    hit_place_holder.BackChannel = -1;
    hit_place_holder.EnergyBack = 0;
    hit_place_holder.EnergyFront = 0;
    hit_place_holder.Energy = 0;
    hit_place_holder.Time = 0;
    hit_place_holder.X = 0;
    hit_place_holder.Y = 0;
    hit_place_holder.Z = 0;
    hit_place_holder.XW = 0;
    hit_place_holder.YW = 0;
    hit_place_holder.ZW = 0;
    hit_place_holder.RW = 0;
    hit_place_holder.PhiW = 0;
    hit_place_holder.RFSubtract = 0;
    hit_place_holder.TrackType = 0;
  };

  void zeroSiHit(){
    NSiHits = 0;

    Detector.clear();
    Hit.clear();
    
    zeroPlaceHolder();
  };
};
//////////////////////////////////////////////////////////////////////////////////////
class PCHit {
  // This class is for Proportional Counter hit
public:
  PCHit(){};

  Int_t NPCHits;

  struct SortByPC {
    Int_t WireID;
    Double_t Down, Up;
    Double_t DownVoltage, UpVoltage;
    Double_t Energy;
    Double_t Z;
    Double_t XW, YW, ZW;
    Double_t RW, PhiW;
    Int_t TrackType;

  }pc_place_holder;

  vector<SortByPC> Hit;
  vector<SortByPC> *ReadHit;

  void zeroPlaceHolder(){
    pc_place_holder.WireID = -1;
    pc_place_holder.Down = 0;
    pc_place_holder.Up = 0;
    pc_place_holder.DownVoltage = 0;
    pc_place_holder.UpVoltage = 0;
    pc_place_holder.Energy = 0;
    pc_place_holder.Z = 0;
    pc_place_holder.XW = 0;
    pc_place_holder.YW = 0;
    pc_place_holder.ZW = 0;
    pc_place_holder.RW = 0;
    pc_place_holder.PhiW = 0;
    pc_place_holder.TrackType = 0;
  };

  void zeroPCHit(){
    NPCHits = 0;
    Hit.clear();
    zeroPlaceHolder();
  };
};

class SETTINGS {
public:
  SETTINGS(){};

  struct set_values {
    Double_t X3Pulser_Offset[24][12];
    Double_t X3Pulser_Slope[24][12];
    Double_t X3RelativeGains[24][12];
    Double_t X3FinalFix[24][12];
    Double_t X3Geometry_Up[24][4][4];
    Double_t X3Geometry_Down[24][4][4];
    Double_t QQQPulser_Offset[4][32];
    Double_t QQQPulser_Slope[4][32];
    Double_t QQQRelativeGains[4][32];
    Double_t QQQFinalFix[4][32];
    Double_t SiAlphaCal[28];

  }set_val;


};
/////////////////////////////////////////////////////////////////////////////////////
class BasicTrack {
 public:
  BasicTrack(){};

  Int_t NTracks;

  struct Track_Event {
    Int_t TrackType;
    Int_t DetID;
    Int_t WireID;
    Double_t SiEnergy;
    Double_t SiZ;
    Double_t SiR;
    Double_t SiPhi;
    Double_t PCEnergy;
    Double_t PCZ;
    Double_t PCR;
    Double_t PCPhi;
    Double_t IntPoint;
    Double_t PathLength;
    Double_t Theta;
  }track_place_holder;

  vector<Track_Event> Event;
  vector<Track_Event> *ReadEvent;


  void zeroPlaceHolder(){
    track_place_holder.TrackType = 0;
    track_place_holder.DetID = -1;
    track_place_holder.WireID = -1;
    track_place_holder.SiEnergy = 0;
    track_place_holder.SiZ = -10;
    track_place_holder.SiR = 0;
    track_place_holder.SiPhi = 0;
    track_place_holder.PCEnergy = 0;
    track_place_holder.PCZ = 0;
    track_place_holder.PCR = 0;
    track_place_holder.PCPhi = 0;
    track_place_holder.IntPoint = 0;
    track_place_holder.PathLength = 0;
    track_place_holder.Theta = 0;
  };

  void zeroBasicTrack(){
    NTracks = 0;
    Event.clear();
    zeroPlaceHolder();
  };

};
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
class Reconstruct {
 public:
  Reconstruct(){};
  Int_t NTracks;
  Int_t NTracks1;
  Int_t NTracks2;
  Int_t NTracks3;

  struct TrackEvent {
    Int_t TrackType;
    Int_t DetID;
    Int_t WireID;
    Double_t SiEnergy;
    Double_t SiZ;
    Double_t SiR;
    Double_t SiPhi;
    Double_t PCEnergy;
    Double_t PCZ;
    Double_t PCR;
    Double_t PCPhi;
    Double_t IntPoint;
    Double_t PathLength;
    Double_t Theta;
    Double_t Phi;
    Double_t BeamEnergy;
  }track_place_holder;

  struct Silicon_Event {
    Int_t TrackType;
    Int_t DetID;
    Double_t SiEnergy;
    Double_t SiZ;
    Double_t SiR;
    Double_t SiPhi;
  }si_place_holder;

  struct PropCounter_Event {
    Int_t TrackType;
    Int_t WireID;
    Double_t PCEnergy;
    Double_t PCZ;
    Double_t PCR;
    Double_t PCPhi;
  }pc_place_holder;

  struct Be8Event {
    Double_t DiffIntPoint;
    Double_t SiEnergy_tot;
    Double_t PCEnergy_tot;
    Double_t Energy_tot;
    Double_t AlphaAngle;
    Double_t BeTheta;
    Double_t BePhi;
    Double_t BeKE;
    Double_t Ex;
    Double_t AvgBeamEnergy;

  }BeEvent;


  vector<TrackEvent> TrEvent;
  vector<TrackEvent> *ReadTrEvent;

  vector<Silicon_Event> SiEvent;
  vector<Silicon_Event> *ReadSiEvent;

  vector<PropCounter_Event> PCEvent;
  vector<PropCounter_Event> *ReadPCEvent;


  void zeroPlaceHolder(){
    track_place_holder.TrackType = 0;
    track_place_holder.DetID = -1;
    track_place_holder.WireID = -1;
    track_place_holder.SiEnergy = 0;
    track_place_holder.SiZ = -10;
    track_place_holder.SiR = 0;
    track_place_holder.SiPhi = 0;
    track_place_holder.PCEnergy = 0;
    track_place_holder.PCZ = 0;
    track_place_holder.PCR = 0;
    track_place_holder.PCPhi = 0;
    track_place_holder.IntPoint = 0;
    track_place_holder.PathLength = 0;
    track_place_holder.Theta = 0;
    track_place_holder.Phi = 0;
    track_place_holder.BeamEnergy = 0;

    si_place_holder.TrackType = 0;
    si_place_holder.DetID = -1;
    si_place_holder.SiEnergy = 0;
    si_place_holder.SiZ = -10;
    si_place_holder.SiR = 0;
    si_place_holder.SiPhi = 0;

    pc_place_holder.TrackType = 0;
    pc_place_holder.WireID = -1;
    pc_place_holder.PCEnergy = 0;
    pc_place_holder.PCZ = 0;
    pc_place_holder.PCR = 0;
    pc_place_holder.PCPhi = 0;

    BeEvent.DiffIntPoint = -10;
    BeEvent.SiEnergy_tot = -10;
    BeEvent.PCEnergy_tot = -10;
    BeEvent.Energy_tot = -10;
    BeEvent.AlphaAngle = -10;
    BeEvent.BeTheta = -10;
    BeEvent.BePhi = -10;
    BeEvent.BeKE = -10;
    BeEvent.Ex = -10;
    BeEvent.AvgBeamEnergy = -10;

  };

  void zeroReconstruct(){
    NTracks = 0;
    NTracks1 = 0;
    NTracks2 = 0;
    NTracks3 = 0;

    TrEvent.clear();
    SiEvent.clear();
    PCEvent.clear();

    zeroPlaceHolder();
  };

};
