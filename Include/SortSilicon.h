#include <TROOT.h>
#include <TRandom3.h>
#include <TMath.h>
#include <algorithm>
#include <vector>

#include "/home/manasta/Desktop/parker_codes/Include/ChannelMap.h"
#include "/home/manasta/Desktop/parker_codes/Include/organizetree.h"

using namespace std;

class SortSilicon {
 public:
  SortSilicon(){};

  struct data {
    Double_t Energy;
    Double_t Time;
    Double_t Channel;
    Double_t Channel_Up;
    Double_t Channel_Down;
    Double_t Energy_Up;
    Double_t Energy_Down;
  };


  data data_place_holder;
  vector<data> front;
  vector<data> back;
  vector<data> front_compress;
  vector<data> back_compress;

  Double_t QQQR;
  Double_t QQQPhi;

  static bool Esort_method(struct data a,struct data b);
  static bool Csort_method(struct data a,struct data b);

  void ProcessQQQ_1(SiHit *Si, ChannelMap *CMAP, Int_t RFTime);
  void ProcessQQQ_2(SiHit *Si, ChannelMap *CMAP, Int_t RFTime);
  void ProcessSX3_1(SiHit *Si, ChannelMap *CMAP, Int_t RFTime);
  void ProcessSX3_2(SiHit *Si, ChannelMap *CMAP, Int_t RFTime);

};

bool SortSilicon::Esort_method(struct data a,struct data b){
    if(a.Energy>b.Energy)
      return 1;
    return 0;
};

bool SortSilicon::Csort_method(struct data a,struct data b){
  if(a.Channel<b.Channel)
    return 1;
  return 0;
};

void SortSilicon::ProcessQQQ_1(SiHit *Si, ChannelMap *CMAP, Int_t RFTime){

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
	    Si->hit_place_holder.Flag = 3;
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
  Si->hit_place_holder.TrackType = 2;
  Si->hit_place_holder.DetID = Si->det_place_holder.DetID;
  if ( Si->hit_place_holder.Flag == 1){//this means one front is matched with one back
    Si->hit_place_holder.HitType = 11;
  }else if ( Si->hit_place_holder.Flag == 2 ){//this means that the front multiplicity is one less
    Si->hit_place_holder.HitType = back_compress.size()*10 + front_compress.size()-1;
  }else if ( Si->hit_place_holder.Flag == 3 ){//this means that the back multiplicity is one less
    Si->hit_place_holder.HitType = (back_compress.size()-1)*10 + front_compress.size();
  }else{
    Si->hit_place_holder.HitType = back_compress.size()*10 + front_compress.size();
  }
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

};

void SortSilicon::ProcessSX3_1(SiHit *Si, ChannelMap *CMAP, Int_t RFTime){
  TRandom3 *XRandom = new TRandom3();
  TRandom3 *ZRandom = new TRandom3();

  XRandom->SetSeed();
  ZRandom->SetSeed();

  front.clear();
  back.clear();
  front_compress.clear();
  back_compress.clear();

  Double_t up[4][3];
  Double_t down[4][3];
  for ( Int_t i=1; i<4; i++ ){
    for( Int_t j=0; j<3; j++ ){
      up[i][j] = 0;
      down[i][j] = 0;
    }
  }
  for (Int_t i=0; i<4; i++){
    up[i][0] = -1;
  }

  for ( Int_t i=0; i<Si->det_place_holder.UpMult; i++ ){
    up[Si->det_place_holder.UpChNum[i]][0] = (Double_t)Si->det_place_holder.UpChNum[i];
    up[Si->det_place_holder.UpChNum[i]][1] = Si->det_place_holder.EnergyUp_Cal[i];
    up[Si->det_place_holder.UpChNum[i]][2] = Si->det_place_holder.TimeUp[i];
  }
  for ( Int_t i=0; i<Si->det_place_holder.DownMult; i++ ){
    down[Si->det_place_holder.DownChNum[i]][0] = (Double_t)Si->det_place_holder.DownChNum[i];
    down[Si->det_place_holder.DownChNum[i]][1] = Si->det_place_holder.EnergyDown_Cal[i];
    down[Si->det_place_holder.DownChNum[i]][2] = Si->det_place_holder.TimeDown[i];
  }

  //copy contents of det_place_holder into new arrays, get rid of bad energies
  for ( Int_t i=0; i<4; i++ ){
    if ( up[i][1]>0 || down[i][1]>0 ){
      data_place_holder.Channel = up[i][0];
      data_place_holder.Channel_Up = up[i][0];
      data_place_holder.Channel_Down = down[i][0];
      data_place_holder.Energy = up[i][1] + down[i][1];
      data_place_holder.Energy_Up = up[i][1];
      data_place_holder.Energy_Down = down[i][1];
      data_place_holder.Time = up[i][2];
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
  sort( back.begin(), back.end(), Csort_method );//front channels already sorted
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

  Double_t up_mult_compress = 0;
  Double_t down_mult_compress = 0;
  for ( Int_t i=0; i<front_compress.size(); i++ ){

  }
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
  //-----------------------------------------------------------------------------------------
  //determine x and z positions
  Si->hit_place_holder.Z = (3-avg_back_ch + ZRandom->Rndm())*1.875;
  Si->hit_place_holder.X = XRandom->Rndm()+(3-avg_front_ch);
  Double_t xw=0, yw=0, zw=0, rw=0, phiw=0;
  CMAP->GetX3WorldCoordinates(Si->det_place_holder.DetID,Si->hit_place_holder.X,Si->hit_place_holder.Z,xw,yw,zw,rw,phiw);
  //-------------------------------------------------------------
  //fill tree
  Si->hit_place_holder.NHitsInDet = 1;
  Si->hit_place_holder.TrackType = 2;
  Si->hit_place_holder.DetID = Si->det_place_holder.DetID;
  Si->hit_place_holder.HitType = 111;
  Si->hit_place_holder.FrontChannel = avg_front_ch;
  Si->hit_place_holder.BackChannel = avg_back_ch;
  Si->hit_place_holder.EnergyFront = (front_sum);
  Si->hit_place_holder.EnergyBack = (back_sum);
  Si->hit_place_holder.Energy = (back_sum);
  //cout << back_compress.at(0).Time << endl;
  Si->hit_place_holder.Time = (avg_time);
  Si->hit_place_holder.RFSubtract = (Si->hit_place_holder.Time - RFTime);

  Si->hit_place_holder.XW = (xw);
  Si->hit_place_holder.YW = (yw);
  Si->hit_place_holder.ZW = (zw);
  Si->hit_place_holder.RW = (rw);
  Si->hit_place_holder.PhiW = (phiw);

  delete XRandom;
  delete ZRandom;
};

void SortSilicon::ProcessSX3_2(SiHit *Si, ChannelMap *CMAP, Int_t RFTime){
  TRandom3 *XRandom = new TRandom3();
  TRandom3 *ZRandom = new TRandom3();

  XRandom->SetSeed();
  ZRandom->SetSeed();

  Si->hit_place_holder.NHitsInDet = 1;
  Si->hit_place_holder.DetID = Si->det_place_holder.DetID;
  Si->hit_place_holder.HitType = Si->det_place_holder.BackMult*100 + Si->det_place_holder.UpMult*10 + Si->det_place_holder.DownMult;
  Si->hit_place_holder.TrackType = 2;
  if ( Si->det_place_holder.HitType == 111 ){
    if ( Si->det_place_holder.DownChNum[0]==Si->det_place_holder.UpChNum[0] ){
      Si->hit_place_holder.EnergyBack = Si->det_place_holder.EnergyBack_Cal[0];
      Si->hit_place_holder.EnergyFront = Si->det_place_holder.EnergyUp_Cal[0] + Si->det_place_holder.EnergyDown_Cal[0];
      Si->hit_place_holder.Energy = Si->det_place_holder.EnergyBack_Cal[0];
      Si->hit_place_holder.Time = Si->det_place_holder.TimeBack[0];
      Si->hit_place_holder.RFSubtract = Si->hit_place_holder.Time - RFTime;
      Double_t FinalZUp=0,FinalZUpCal=-10, FinalZDown=0, FinalZDownCal=-10;
      FinalZDown = 2*Si->det_place_holder.EnergyDown_Cal[0]/Si->det_place_holder.EnergyBack_Cal[0]-1;
      CMAP->PosCal(Si->hit_place_holder.DetID,Si->det_place_holder.DownChNum[0],Si->det_place_holder.BackChNum[0],FinalZDown,FinalZDownCal);
      FinalZUp = 1-2*Si->det_place_holder.EnergyUp_Cal[0]/Si->det_place_holder.EnergyBack_Cal[0];
      CMAP->PosCal(Si->hit_place_holder.DetID,Si->det_place_holder.UpChNum[0],Si->det_place_holder.BackChNum[0],FinalZUp,FinalZUpCal);
      if (FinalZDownCal > -1 && FinalZUpCal > -1){
	Si->hit_place_holder.Z = FinalZUpCal;
	Si->hit_place_holder.X = XRandom->Rndm()+(3-Si->det_place_holder.UpChNum[0]);
      }else{
	Si->hit_place_holder.HitType += 1000;
      }
      Double_t xw=0, yw=0, zw=0, rw=0, phiw=0;
      CMAP->GetX3WorldCoordinates(Si->hit_place_holder.DetID,Si->hit_place_holder.X,Si->hit_place_holder.Z,xw,yw,zw,rw,phiw);
      Si->hit_place_holder.XW = xw;
      Si->hit_place_holder.YW = yw;
      Si->hit_place_holder.ZW = zw;
      Si->hit_place_holder.RW= rw;
      Si->hit_place_holder.PhiW = phiw;
      MyFill(Form("ZPos%i_%i_%i",Si->det_place_holder.DetID,Si->det_place_holder.UpChNum[0],Si->det_place_holder.BackChNum[0]),800,-1,1,FinalZDown);
    }else{
      Si->hit_place_holder.HitType += 1000;
    }
  }else if ( (Si->det_place_holder.HitType == 110) || (Si->det_place_holder.HitType == 101) ){
    Si->hit_place_holder.EnergyBack = (Si->det_place_holder.EnergyBack_Cal[0]);
    Si->hit_place_holder.Energy = (Si->det_place_holder.EnergyBack_Cal[0]);
    Si->hit_place_holder.Time = (Si->det_place_holder.TimeBack[0]);
    Si->hit_place_holder.RFSubtract = (Si->hit_place_holder.Time - RFTime);
    Double_t FinalZUp=0,FinalZUpCal=-10, FinalZDown=0, FinalZDownCal=-10;
    if ( Si->det_place_holder.HitType == 110 ){
      Si->hit_place_holder.EnergyFront = (Si->det_place_holder.EnergyUp_Cal[0]);
      FinalZUp = 1-2*Si->det_place_holder.EnergyUp_Cal[0]/Si->det_place_holder.EnergyBack_Cal[0];
      CMAP->PosCal(Si->hit_place_holder.DetID,Si->det_place_holder.UpChNum[0],Si->det_place_holder.BackChNum[0],FinalZUp,FinalZUpCal);
      if (FinalZUpCal > -1){
	Si->hit_place_holder.Z = (FinalZUpCal);
	Si->hit_place_holder.X = (XRandom->Rndm()+(3-Si->det_place_holder.UpChNum[0]));
      }else{
	Si->hit_place_holder.HitType += 1000;
      }
     MyFill(Form("ZPos%i_%i_%i",Si->det_place_holder.DetID,Si->det_place_holder.UpChNum[0],Si->det_place_holder.BackChNum[0]),800,-1,1,FinalZUp);
    }else{
      Si->hit_place_holder.EnergyFront = (Si->det_place_holder.EnergyDown_Cal[0]);
      FinalZDown = 2*Si->det_place_holder.EnergyDown_Cal[0]/Si->det_place_holder.EnergyBack_Cal[0]-1;
      CMAP->PosCal(Si->hit_place_holder.DetID,Si->det_place_holder.DownChNum[0],Si->det_place_holder.BackChNum[0],FinalZDown,FinalZDownCal);
      if (FinalZDownCal > -1){
	Si->hit_place_holder.Z = (FinalZDownCal);
	Si->hit_place_holder.X = (XRandom->Rndm()+(3-Si->det_place_holder.DownChNum[0]));
      }else{
	Si->hit_place_holder.HitType += 1000;
      }
     MyFill(Form("ZPos%i_%i_%i",Si->det_place_holder.DetID,Si->det_place_holder.DownChNum[0],Si->det_place_holder.BackChNum[0]),800,-1,1,FinalZDown);
    }
    Double_t xw=0, yw=0, zw=0, rw=0, phiw=0;
    CMAP->GetX3WorldCoordinates(Si->hit_place_holder.DetID,Si->hit_place_holder.X,Si->hit_place_holder.Z,xw,yw,zw,rw,phiw);
    Si->hit_place_holder.XW = (xw);
    Si->hit_place_holder.YW = (yw);
    Si->hit_place_holder.ZW = (zw);
    Si->hit_place_holder.RW = (rw);
    Si->hit_place_holder.PhiW = (phiw);
  }else if ( Si->det_place_holder.HitType == 211 ){
    if ( Si->det_place_holder.DownChNum[0]==Si->det_place_holder.UpChNum[0] ){
      if ( fabs(Si->det_place_holder.BackChNum[1]-Si->det_place_holder.BackChNum[0]) == 1 ){
	Si->hit_place_holder.EnergyBack = (Si->det_place_holder.EnergyBack_Cal[0] + Si->det_place_holder.EnergyBack_Cal[1] );
	Si->hit_place_holder.EnergyFront = (Si->det_place_holder.EnergyUp_Cal[0] + Si->det_place_holder.EnergyDown_Cal[0]);
	Si->hit_place_holder.Energy = (Si->det_place_holder.EnergyBack_Cal[0] + Si->det_place_holder.EnergyBack_Cal[1]);
	Si->hit_place_holder.Time = (Si->det_place_holder.TimeBack[0]);
	Si->hit_place_holder.RFSubtract = (Si->hit_place_holder.Time - RFTime);
	
	Si->hit_place_holder.Z = (1.875*(3-Si->det_place_holder.BackChNum[0]) - 0.1 + 0.2*ZRandom->Rndm());//assuming a mm inbetween strips
	Si->hit_place_holder.X = (XRandom->Rndm()+(3-Si->det_place_holder.DownChNum[0]));
	
	Double_t xw=0, yw=0, zw=0, rw=0, phiw=0;
	CMAP->GetX3WorldCoordinates(Si->hit_place_holder.DetID,Si->hit_place_holder.X,Si->hit_place_holder.Z,xw,yw,zw,rw,phiw);
	Si->hit_place_holder.XW = (xw);
	Si->hit_place_holder.YW = (yw);
	Si->hit_place_holder.ZW = (zw);
	Si->hit_place_holder.RW = (rw);
	Si->hit_place_holder.PhiW = (phiw);
      }else{
	Si->hit_place_holder.HitType += 1000;
      }
    }else{
      Si->hit_place_holder.HitType += 1000;
    }
  }else{
    Si->hit_place_holder.HitType += 1000;
  }
  delete XRandom;
  delete ZRandom;
};

void SortSilicon::ProcessQQQ_2(SiHit *Si, ChannelMap *CMAP, Int_t RFTime){

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

  if ( Si->det_place_holder.HitType == 11 ){
    QQQPhi = (Si->det_place_holder.BackChNum[0]+PhiRandom->Rndm())*StripAngle;
    QQQR = OuterRadius - (Si->det_place_holder.UpChNum[0]+RRandom->Rndm() )*RingPitch;
    Si->hit_place_holder.X = QQQR*TMath::Cos(QQQPhi);
    Si->hit_place_holder.Y = QQQR*TMath::Sin(QQQPhi);
    Si->hit_place_holder.Energy = Si->det_place_holder.EnergyBack_Cal[0];
    Si->hit_place_holder.EnergyBack = Si->det_place_holder.EnergyBack_Cal[0];
    Si->hit_place_holder.EnergyFront = Si->det_place_holder.EnergyUp_Cal[0];
    Si->hit_place_holder.FrontChannel = (Double_t)Si->det_place_holder.UpChNum[0];
    Si->hit_place_holder.BackChannel = (Double_t)Si->det_place_holder.BackChNum[0];
    Si->hit_place_holder.DetID = Si->det_place_holder.DetID;
    Si->hit_place_holder.Time = Si->det_place_holder.TimeBack[0];
    Si->hit_place_holder.RFSubtract = Si->hit_place_holder.Time - RFTime;
    CMAP->GetQQQ3WorldCoordinates(Si->hit_place_holder.DetID, Si->hit_place_holder.X,Si->hit_place_holder.Y,Si->hit_place_holder.XW,Si->hit_place_holder.YW, Si->hit_place_holder.RW, Si->hit_place_holder.PhiW);
    Si->hit_place_holder.ZW = 0;//JJPIV--6Jan2016--explicitly setting w
    Si->hit_place_holder.HitType = 11;
    Si->hit_place_holder.TrackType = 2;
  }else if ( Si->det_place_holder.HitType == 12 ){    
    if ( fabs(Si->det_place_holder.UpChNum[1] - Si->det_place_holder.UpChNum[0]) == 1){//if front channels are adjacent
      Double_t eff_channel = (Double_t)(Si->det_place_holder.UpChNum[1] + Si->det_place_holder.UpChNum[0])/2;
      QQQPhi = (Si->det_place_holder.BackChNum[0]+PhiRandom->Rndm())*StripAngle;
      QQQR = (OuterRadius - 0.5*RingPitch) - eff_channel*RingPitch - (0.5-RRandom->Rndm())*0.02;//assume that we know position to better resolution than channel width, because particle must've hit close to the border to share charge...

      Si->hit_place_holder.X = QQQR*TMath::Cos(QQQPhi);
      Si->hit_place_holder.Y = QQQR*TMath::Sin(QQQPhi);
      Si->hit_place_holder.Energy = Si->det_place_holder.EnergyBack_Cal[0];
      Si->hit_place_holder.EnergyBack = Si->det_place_holder.EnergyBack_Cal[0];
      Si->hit_place_holder.EnergyFront = Si->det_place_holder.EnergyUp_Cal[0] + Si->det_place_holder.EnergyUp_Cal[1];
      Si->hit_place_holder.FrontChannel = eff_channel;
      Si->hit_place_holder.BackChannel = (Double_t)Si->det_place_holder.BackChNum[0];
      Si->hit_place_holder.DetID = Si->det_place_holder.DetID;
      Si->hit_place_holder.Time = Si->det_place_holder.TimeBack[0];
      Si->hit_place_holder.RFSubtract = Si->hit_place_holder.Time - RFTime;
      CMAP->GetQQQ3WorldCoordinates(Si->hit_place_holder.DetID, Si->hit_place_holder.X,Si->hit_place_holder.Y,Si->hit_place_holder.XW,Si->hit_place_holder.YW, Si->hit_place_holder.RW, Si->hit_place_holder.PhiW);
      Si->hit_place_holder.ZW = 0;//JJPIV--6Jan2016--explicitly setting w
      Si->hit_place_holder.HitType = 12;
      Si->hit_place_holder.TrackType = 2;
    }else{
      Si->hit_place_holder.HitType = 10012;
    }
  }else if ( Si->det_place_holder.HitType == 21 ){    
    if ( fabs(Si->det_place_holder.BackChNum[1] - Si->det_place_holder.BackChNum[0]) == 1){//if back channels are adjacent
      Double_t eff_channel = (Double_t)(Si->det_place_holder.BackChNum[1] + Si->det_place_holder.BackChNum[0])/2;
      QQQPhi = (0.5+eff_channel)*StripAngle - (0.5-PhiRandom->Rndm())*0.02 ;
      QQQR = OuterRadius - (Si->det_place_holder.UpChNum[0]+RRandom->Rndm() )*RingPitch;

      Si->hit_place_holder.X = QQQR*TMath::Cos(QQQPhi);
      Si->hit_place_holder.Y = QQQR*TMath::Sin(QQQPhi);
      Si->hit_place_holder.Energy = Si->det_place_holder.EnergyBack_Cal[0] + Si->det_place_holder.EnergyBack_Cal[1];
      Si->hit_place_holder.EnergyBack = Si->det_place_holder.EnergyBack_Cal[0] + Si->det_place_holder.EnergyBack_Cal[1];
      Si->hit_place_holder.EnergyFront = Si->det_place_holder.EnergyUp_Cal[0];
      Si->hit_place_holder.FrontChannel = (Double_t)Si->det_place_holder.UpChNum[0];
      Si->hit_place_holder.BackChannel = eff_channel;
      Si->hit_place_holder.DetID = Si->det_place_holder.DetID;
      Si->hit_place_holder.Time = Si->det_place_holder.TimeBack[0];
      Si->hit_place_holder.RFSubtract = Si->hit_place_holder.Time - RFTime;
      CMAP->GetQQQ3WorldCoordinates(Si->hit_place_holder.DetID, Si->hit_place_holder.X,Si->hit_place_holder.Y,Si->hit_place_holder.XW,Si->hit_place_holder.YW, Si->hit_place_holder.RW, Si->hit_place_holder.PhiW);
      Si->hit_place_holder.ZW = 0;//JJPIV--6Jan2016--explicitly setting w
      Si->hit_place_holder.HitType = 21;
      Si->hit_place_holder.TrackType = 2;
    }else{
      Si->hit_place_holder.HitType = 10021;
    }
  }
  
  delete PhiRandom;
  delete RRandom;
  
};
