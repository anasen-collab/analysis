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

    vector<Int_t> HitType;
    vector<Double_t> EnergyBack, EnergyFront;
    vector<Double_t> Energy, Time;
    vector<Double_t> X, Y, Z;
    vector<Double_t> XW, YW, ZW;
    vector<Double_t> RW, PhiW;
    vector<Double_t> RFSubtract;

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
    hit_place_holder.DetID = 0;

    hit_place_holder.HitType.clear();
    hit_place_holder.EnergyBack.clear();
    hit_place_holder.EnergyFront.clear();
    hit_place_holder.Energy.clear();
    hit_place_holder.Time.clear();
    hit_place_holder.X.clear();
    hit_place_holder.Y.clear();
    hit_place_holder.Z.clear();
    hit_place_holder.XW.clear();
    hit_place_holder.YW.clear();
    hit_place_holder.ZW.clear();
    hit_place_holder.RW.clear();
    hit_place_holder.PhiW.clear();
    hit_place_holder.RFSubtract.clear();
  };

  void zeroSiHit()
  {
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
  Int_t NPCHits,WireID[MaxPCHits];
  Double_t Down[MaxPCHits],Up[MaxPCHits];
  Double_t DownVoltage[MaxPCHits], UpVoltage[MaxPCHits];
  Double_t Energy[MaxPCHits];
  Double_t Zp[MaxPCHits];
  Double_t XWp[MaxPCHits],YWp[MaxPCHits],ZWp[MaxPCHits];
  Double_t RWp[MaxPCHits], PhiWp[MaxPCHits];

  void zeroPCHit(){
    NPCHits = 0;
    for (Int_t i=0; i<MaxPCHits; i++) {
      WireID[i] = -1;
      Down[i] = 0;
      Up[i]   = 0;
      DownVoltage[i] = 0;
      UpVoltage[i] = 0;
      Energy[i] = 0;
      Zp[i] = 0;
      XWp[i] = 0;
      YWp[i] = 0;
      ZWp[i] = 0;
      RWp[i] = 0;
      PhiWp[i] = 0;
    }
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
class Track {
 public:
  Track(){};
  Int_t NTracks;
  Int_t NTracks1;
  Int_t NTracks2;
  Int_t NTracks3;
  Int_t TrackType[MaxTracks];
  Double_t SiZ[MaxTracks];
  Double_t SiR[MaxTracks];
  Double_t SiPhi[MaxTracks];
  Double_t SiEnergy[MaxTracks];
  Double_t PCZ[MaxTracks];
  Double_t PCZ_Ref[MaxTracks];
  Double_t PCR[MaxTracks];
  Double_t PCPhi[MaxTracks];
  Double_t PCEnergy[MaxTracks];
  Double_t PCUpVoltage[MaxTracks];
  Double_t PCDownVoltage[MaxTracks];
  Double_t IntPoint[MaxTracks];
  Double_t Theta[MaxTracks];
  Double_t PathLength[MaxTracks];
  Double_t BeamEnergy[MaxTracks];
  Double_t EnergyLoss[MaxTracks];
  Double_t DiffIntPoint;
  Double_t DiffIntPoint_Proton;
  Double_t SiEnergy_tot;
  Double_t PCEnergy_tot;
  Double_t Energy_tot;
  Double_t Ex;
  Double_t Angle;
  Int_t WireID[MaxTracks];
  Int_t DetID[MaxTracks];
  Int_t ChNum[MaxTracks];
  Int_t HitType[MaxTracks];
  Int_t HitType_Old[MaxTracks];
  Double_t DiffPhi[MaxTracks];
  Double_t DiffZ[MaxTracks];
  Double_t DiffR[MaxTracks];
  Double_t BeTheta;
  Double_t BeKE;
  Double_t AvgBeamEnergy;
  Double_t ProtonEnergy;
  Double_t ProtonEnergy_Rec;
  Double_t ProtonDiffEnergy;
  Double_t ProtonAngle;

  void zeroTrack(){
    NTracks = 0;
    NTracks1 = 0;
    NTracks2 = 0;
    NTracks3 = 0;
    DiffIntPoint = -10;
    DiffIntPoint_Proton = -10;
    SiEnergy_tot = 0;
    PCEnergy_tot = 0;
    Energy_tot = 0;
    Ex = -10;
    Angle = -10;
    BeTheta = -10;
    BeKE = -10;
    AvgBeamEnergy = -10;
    ProtonEnergy = -10;
    ProtonEnergy_Rec = -10;
    ProtonDiffEnergy = -1000;
    ProtonAngle = -1000;

    for (Int_t i=0; i<MaxTracks; i++){
      TrackType[i] = -1;
      SiZ[i] = 0;
      SiR[i] = 0;
      SiPhi[i] = 0;
      SiEnergy[i] = 0;
      PCZ[i] = 0;
      PCZ_Ref[i] = 0;
      PCR[i] = 0;
      PCPhi[i] = 0;
      PCEnergy[i] = 0;
      PCUpVoltage[i] = 0;
      PCDownVoltage[i] = 0;
      IntPoint[i] = 0;
      Theta[i] = 0;
      PathLength[i] = 0;
      BeamEnergy[i] = 0;
      EnergyLoss[i] = 0;
      WireID[i] = -1;
      DetID[i] = -1;
      ChNum[i] = -1;
      HitType[i] = -1;
      HitType_Old[i] = -1;
      DiffPhi[i] = -10;
      DiffZ[i] = -10;
      DiffR[i] = -10;
    }

  };

};
/////////////////////////////////////////////////////////////////////////////////////

