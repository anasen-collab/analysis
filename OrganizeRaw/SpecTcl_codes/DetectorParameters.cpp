// Ingo Wiedenhoever 2002,2006
// Eric Diffenderfer 7/2005
#include "DetectorParameters.hpp"

DetectorPar::DetectorPar(const std::string& name,
			 const unsigned int& num_chan):
  num_chan_(num_chan),
  e(name+".e",4096,0.,4096.,"ch",num_chan,0),
  t(name+".t",4096,0.,4096.,"ch",num_chan,0),
  chlist(name+".ch",6,num_chan,0),
  mult(name+".mult",6)
{}

std::string DetectorPar::print()
{
  std::string out;
  std::ostringstream tmp;
  for(int i=0;i<num_chan_;i++){
    if(e[i].isValid())
      tmp << e[i].getName() << " " << e[i] << std::endl;
    if(t[i].isValid())
      tmp << t[i].getName() << " " << t[i] << std::endl;
  }
  out = tmp.str();
  return out;
}

TempList::TempList(const unsigned int no_channels){
  mult=0;
  no_channels_=no_channels;
  chlist = new double[no_channels];
  elist = new double[no_channels];
  tlist = new double[no_channels];
}

TempList::~TempList(){
  delete[] chlist;
  delete[] elist;
  delete[] tlist;
}

 void TempList::InsertHit(const double& e,const double& t,const double& ch)
{ 
  int i,j;
  /* enough space ? */
  if (mult>=no_channels_){
    return;
  }
  /* insert into list sorted by energy */
  for (i=(int)mult-1;i>=0;i--){
    if (e<elist[i])
      break;
  }
  // element i+1 is at the position for ch 
  // so make room for it
  for (j=(int)mult-1;j>i;j--){
    chlist[j+1]=(double)chlist[j];
    elist[j+1]=(double)elist[j];
    tlist[j+1]=(double)tlist[j];
  }
  // and shove it in
  chlist[i+1]=ch;
  elist[i+1]=e;
  tlist[i+1]=t;
  mult ++;
}


HitList::HitList(const std::string& basename,
		 const unsigned int& no_channels,
		 const int& sorted_by_channels):
  no_channels_(no_channels),
  mult(basename+".mult",6),
  chlist(basename+".ch",6,no_channels,0),
  elist(basename+".e",4096,ECAL_LOW,ECAL_HIGH,"MeV",no_channels,0),
  tlist(basename+".t",4096,TCAL_LOW,TCAL_HIGH,"ns",no_channels,0),
  sorted_by_channel(sorted_by_channels)
{}

std::string HitList::print()
{
  std::string out;
  std::ostringstream tmp;
  for(int i=0;i<no_channels_;i++){
    if(chlist[i].isValid())
      tmp << chlist[i].getName() << " " << chlist[i] << std::endl;
    if(elist[i].isValid())
      tmp << elist[i].getName() << " " << elist[i] << std::endl;
    if(tlist[i].isValid())
      tmp << tlist[i].getName() << " " << tlist[i] << std::endl;
  }
  out = tmp.str();
  return out;
}

 void HitList::InsertHit(const double& ch,const double& e,const double& t)
{ 
  int i,j;
  // is the list full already ?
  if (mult>=no_channels_){
    return;
  }
    
  if (sorted_by_channel){ 
    /* insert into list sorted by channel */
    for (i=(int)mult-1;i>=0;i--){
      if (ch>chlist[i])
	break;
    }
  }
  else{ /* sorted by energy, largest e first */
    for (i=(int)mult-1;i>=0;i--){
      if (e<elist[i])
	break;
    }
  }   
  // element i+1 is at the position for ch 
  // so make room for it
  for (j=(int)mult-1;j>i;j--){
    chlist[j+1]=(double)chlist[j];
    elist[j+1]=(double)elist[j];
    tlist[j+1]=(double)tlist[j];
  }
  // and shove it in
  chlist[i+1]=ch;
  elist[i+1]=e;
  tlist[i+1]=t;
  mult += 1.;
}

HitList::~HitList()
{
}

ClusterList::ClusterList(const std::string& basename,
			 const unsigned int& no_channels):
  no_channels_(no_channels),
  mult(basename+".mult",6),
  elist(basename+".e",4096,ECAL_LOW,ECAL_HIGH,"MeV",no_channels,0),
  tlist(basename+".t",4096,TCAL_LOW,TCAL_HIGH,"ns",no_channels,0),
  chflist(basename+".cf",6,no_channels,0),
  chblist(basename+".cb",6,no_channels,0),
  xcrdlist(basename+".xcrd",4096,-300.,300.,"mm",no_channels,0),
  ycrdlist(basename+".ycrd",4096,-300.,300.,"mm",no_channels,0),
  zcrdlist(basename+".zcrd",4096,-300.,300.,"mm",no_channels,0)
{};

ClusterList::~ClusterList()
{
};

std::string ClusterList::print()
{
  std::string out;
  std::ostringstream tmp;
  for(int i=0;i<no_channels_;i++){
    if(chflist[i].isValid())
      tmp << chflist[i].getName() << " " << chflist[i] << std::endl;
    if(chblist[i].isValid())
      tmp << chblist[i].getName() << " " << chblist[i] << std::endl;
    if(elist[i].isValid())
      tmp << elist[i].getName() << " " << elist[i] << std::endl;
    if(tlist[i].isValid())
      tmp << tlist[i].getName() << " " << tlist[i] << std::endl;
  }
  out = tmp.str();
  return out;
}

void ClusterList::InsertHit(double e,double t,double chf,double chb
			    ,double xcrd=0.,double ycrd=0.,double zcrd=0.)
{ int i,j;
/* enough space for another ? */
 if (mult>=no_channels_){
   return;
 } 
 
  /* insert into list sorted by energy, large first */
  for (i=(int)mult-1;i>=0;i--){
    if (e<elist[i])
      break;
  }
  // element i+1 is at the position for ch 
  // so make room for it
  for (j=(int)mult-1;j>i;j--){
    chflist[j+1]=(double)chflist[j];
    chblist[j+1]=(double)chblist[j];
    elist[j+1]=(double)elist[j];
    tlist[j+1]=(double)tlist[j];
    xcrdlist[j+1]=(double)xcrdlist[j];
    ycrdlist[j+1]=(double)ycrdlist[j];
    zcrdlist[j+1]=(double)zcrdlist[j];
  }
  // and shove it in
  chflist[i+1]=chf;
  chblist[i+1]=chb;
  elist[i+1]=e;
  tlist[i+1]=t;
  xcrdlist[j+1]=xcrd;
  ycrdlist[j+1]=ycrd;
  zcrdlist[j+1]=zcrd;
  mult += 1.;
}


////////////////////////////////////////
///////////    S2Detector    ///////////
////////////////////////////////////////

S2DetectorPar::S2DetectorPar(const std::string& basename,
			     const unsigned int& num_front,
			     const unsigned int& num_back,
			     const Vector3& posv):
  num_front_(num_front),
  num_back_(num_back),
  posv_(posv),
  shiftv_(0.,0.,0.),
  rotv_(0.,0.,0.),
  fmult(basename+".fmult",6),
  bmult(basename+".bmult",6),
  patter_e(basename+".patter_e",6),
  patter_t(basename+".patter_t",6),
  e(basename+".e",4096,ECAL_LOW,ECAL_HIGH,"MeV",num_front+num_back,0),
  t(basename+".t",4096,TCAL_LOW,TCAL_HIGH,"ns",num_front+num_back,0),
//  relt(name+".relt",4096,0.,220.,"ns"), //DSG
  efrontmatch(basename+".efrontmatch",4096,ECAL_LOW,ECAL_HIGH,"MeV"),
  ebackmatch(basename+".ebackmatch",4096,ECAL_LOW,ECAL_HIGH,"MeV"),
  efronttotal(basename+".efronttotal",4096,ECAL_LOW,ECAL_HIGH,"MeV"),
  ebacktotal(basename+".ebacktotal",4096,ECAL_LOW,ECAL_HIGH,"MeV"),
  frontmatchstat(basename+".frontmatchstat",256,0,256,""),
  backmatchstat(basename+".backmatchstat",256,0,256,""),
  f_total(basename+".f_total",num_front,1 /* sorted by channels*/ ),
  b_total(basename+".b_total",num_back,1 /* sorted by channels*/ )
{
  ring_pitch_ = (S2OUTERRAD - S2INNERRAD) / static_cast<double>(num_front_);
  delta_phi_ = 360. / static_cast<double>(num_back_);
  match_epsilon=0;
  match_delta=0.1;
  match_maxene=4096;
  match_enefromback=0;
  addback_front=0;
  addback_back=0;
}

 int S2DetectorPar::SetMatchParameters(int match_enefromback,
				       double match_epsilon,
				       double match_delta,
				       double match_maxene,
				       int addback_front,
				       int addback_back)
{
  this->match_enefromback=match_enefromback;
  this->match_epsilon=match_epsilon;
  this->match_delta=match_delta;
  this->match_maxene=match_maxene;
  this->addback_front=addback_front;
  this->addback_back=addback_back;
  return 0;
}


std::string S2DetectorPar::print()
{
  std::string out;
  std::ostringstream tmp;
  for(int i=0;i<num_front_+num_back_;i++){
    if(e[i].isValid())
      tmp << e[i].getName() << " " << e[i] << std::endl;
    if(t[i].isValid())
      tmp << t[i].getName() << " " << t[i] << std::endl;
  }
  tmp << f_total.print() << b_total.print();
  out = tmp.str();
  return out;
}

void S2DetectorPar::reset()
{
  fmult = 0;
  bmult = 0;
  for(int i=0;i<num_front_+num_back_;i++){
    /* Clear the TreeParameters which may not be set by Unpacker */
    e[i].Reset();
    t[i].Reset();
  }
  f_total.reset();
  b_total.reset();
  no_hits = 0;
}

Vector3 S2DetectorPar::chVect(const double& cf,const double& cb){
  //First make a vector to the channel for a detector at the origin
  double s;
  double phi;
  if((cf < static_cast<double>(num_front_)) 
     && (cb < static_cast<double>(num_back_))){
    s = S2INNERRAD + (cf + (((double)rand())/RAND_MAX))*ring_pitch_;
    phi = (static_cast<double>(num_back_) - cb 
	   - (((double)rand())/RAND_MAX))*delta_phi_;
  }
  else{ 
    std::cerr << "InValid channels sent to S2DetectorPar::chVect" << std::endl;
    //    exit(EXIT_FAILURE);
    Vector3 resultv;
    return resultv;
  }

  Vector3 resultv;
  resultv.SetMagThetaPhi(s,M_PI/2.,phi*M_PI/180.);
  
  //then calibrate the vector 
  //add it to the detector position vector
  if(rotv_[0]){
    resultv.RotateX(rotv_[0]*M_PI/180.);
  }
  if(rotv_[1]){
    resultv.RotateY(rotv_[1]*M_PI/180.);
  }
  if(rotv_[2]){
    resultv.RotateZ(rotv_[2]*M_PI/180.);
  }
  resultv = resultv + shiftv_;
  resultv = resultv + posv_;
    
  return resultv;
}

int S2DetectorPar::ReconstructClusters2(HitList& front,
					HitList& back,
					ClusterList& out) 
{ 
  TempList FrontClusters(48),BackClusters(16);
  int nback = static_cast<int>(back.mult);
  if(!nback) return 0;
  int i,j,front_match;
  unsigned int frontstatus=0;
  unsigned int backstatus=0;

  efronttotal=0.;
  ebacktotal=0.;
  // first add-back neighboring hits
  // in the back
  for (i=0;i<back.mult;i++){
    backstatus |= (((unsigned int)back.chlist[i]));
    ebacktotal = ebacktotal+back.elist[i];
    if (addback_back&&
	(i+1<back.mult)&&(back.chlist[i+1]==back.chlist[i]+1)){
      double e=back.elist[i+1]+back.elist[i];
      double t=(back.tlist[i+1]+back.tlist[i])/2;
      double ch=(back.chlist[i+1]+back.chlist[i])/2.0;
      BackClusters.InsertHit(e,t,ch);
      // we analyzed the i+1 hit along with the i, so 
      i++;
    }
    else{
      double e=back.elist[i];
      double t=back.tlist[i];
      double ch=back.chlist[i];
      BackClusters.InsertHit(e,t,ch);
    }
  }
  //in the front
  for (i=0;i<front.mult;i++){
    efronttotal = efronttotal+front.elist[i];
  }
  for (i=0;i<front.mult;i++){
    frontstatus |= (((unsigned int)front.chlist[i]));
    if (addback_front&&
	(i+1<front.mult)&&(front.chlist[i+1]==front.chlist[i]+1)){
      double e=front.elist[i+1]+front.elist[i];
      double t=(front.tlist[i+1]+front.tlist[i])/2;
      double ch=(front.chlist[i+1]+front.chlist[i])/2.0;
      FrontClusters.InsertHit(e,t,ch);
      // we analyzed the i+1 hit along with the i, so 
      i++;
    }
    else{
      double e=front.elist[i];
      double t=front.tlist[i];
      double ch=front.chlist[i];
      FrontClusters.InsertHit(e,t,ch);
    }
  }
  front_match=0;i=0;
  efrontmatch=0.;ebackmatch=0.;
  // Front and Backclusters are sorted by energy, high first
  while ((i<BackClusters.mult)&&(front_match<FrontClusters.mult)){
    double eps = FrontClusters.elist[front_match]-BackClusters.elist[i];
    double delta = match_epsilon+match_delta*BackClusters.elist[i];
    double back_e=BackClusters.elist[i];
    if (back_e>match_maxene || (fabs(eps)<=delta)){
      // we have a match
      double match_e= FrontClusters.elist[front_match];
      double match_ch= FrontClusters.chlist[front_match];
      double cb=BackClusters.chlist[i];
      double back_t=BackClusters.tlist[i];
      Vector3 particlvect = chVect(match_ch,cb);
      double xcrd = particlvect.X();
      double ycrd = particlvect.Y();
      double zcrd = particlvect.Z();
      if (match_enefromback){
	out.InsertHit(back_e,back_t,
		      match_ch,cb,xcrd,ycrd,zcrd);
      }
      else{
	out.InsertHit(match_e,back_t,
		      match_ch,cb,xcrd,ycrd,zcrd);
      }// end if "match_enefromback" 
      {	
	// sum up matched energies front and back
	efrontmatch += match_e;
	ebackmatch += back_e;
      }
    }// end found match
    else if (eps>delta){
      // front energy is larger than back energy+delta , 
      // so we won't ever match 
      // this one in the future, drop it, pick next front, keep the back
      frontstatus|=128; // flags "unmatched" 
      front_match++;
      continue;
    }
    else if (eps< (-delta)){
      // back energy is larger than front energy-delta , 
      // so we won't ever match 
      // this one in the future, drop it, pick next back, keep the front
      backstatus|=128; // flags "unmatched" 
      i++;
      continue;
    }

    front_match++;
    i++;
  }
  frontmatchstat=frontstatus;
  backmatchstat=backstatus;
  return(0);
}


int S2DetectorPar::ReconstructClusters(HitList& front,
				       HitList& back,
				       ClusterList& out) 
{ 
  int nback = static_cast<int>(back.mult);
  double threshold=match_delta;
  if(!nback) return 0;
  
  //front energy and back energy 
  //must be less than eps apart for a match

  //first the simple case of only one back-side hit
  if(nback == 1){
    double eps = match_epsilon+match_delta*back.elist[0];
    double match_ch = -1;
    double match_e = 0.;
    double match_e_back = back.elist[0];
    //look for single channel matches
    for(int i=0;i<front.mult;i++){
      double t_eps = fabs(front.elist[i] - back.elist[0]);
      if(back.elist[0]>match_maxene || (t_eps < eps)){
	eps = t_eps;
	match_ch = front.chlist[i];
	match_e = front.elist[i];
      }
    }
    //now check if double channel has better match
    if(front.mult > 1){
      for (int i=0;i<front.mult;i++){
	int j;
	for (j=i+1;j<front.mult;j++){
	  // look for "out of sequence" channel -> end of cluster
	  if (front.chlist[j]!=front.chlist[i]+(j-i))
	    break;
	}
	//already covered the lone channel case
	if(j==i+1) continue; 
	
	//loop through the cluster looking at pairs of channels
	for(int k=i;k<j;k++){
	  double sume=0.;
	  double sumc=0.;
	  //calculate summed energy and center of mass for the pair
	  for(int m=k;(m<k+2)&&(m<j);m++){
	    sume += front.elist[m];
	    sumc += front.elist[m]*front.chlist[m];
	  }
	  //compare this pair to previous matches
	  double t_eps = fabs(sume - back.elist[0]);
	  if(t_eps < eps){
	    eps = t_eps;
	    match_ch = sumc/sume;
	    match_e = sume;
	  }
	}
      }
    }
    //if there was a match insert the hit
    if(match_ch >= 0.){
      //calculate position of particle hit and put it into the cluster list
      double cf = match_ch;
      double cb = back.chlist[0];
      if((cf < static_cast<double>(num_front_)) 
	   && (cb < static_cast<double>(num_back_))){
	Vector3 particlvect = chVect(cf,cb);
	double xcrd = particlvect.X();
	double ycrd = particlvect.Y();
	double zcrd = particlvect.Z();
	if (match_enefromback){
	  out.InsertHit(match_e_back,back.tlist[0],
			cf,cb,xcrd,ycrd,zcrd);
	}
	else{
	  out.InsertHit(match_e,back.tlist[0],
		      cf,cb,xcrd,ycrd,zcrd);
	}
      }
    }
  }
  else{ //more than one backside hit

    //declare storage containers
    std::vector<std::vector<double> > eps(nback);
    std::vector<double> erow(1,0.);
    std::vector<std::vector<double> > match_e(nback,erow);
    std::vector<std::vector<double> > match_e_back(nback,erow);
    std::vector<double> chrow(1,-1.);
    std::vector<std::vector<double> > match_ch(nback,chrow);

    //loop over the backside hits
    //building a list of acceptable matches for each
    for(int i=0;i<nback;i++){
      eps[i].push_back(match_epsilon+match_delta*back.elist[i]);
      //first look for single channel matches
      for(int j=0;j<front.mult;j++){
	double t_eps = fabs(front.elist[j] - back.elist[i]);
	if(t_eps < eps[i].back()){
	  eps[i].push_back(t_eps);
	  match_ch[i].push_back(front.chlist[j]);
	  match_e[i].push_back(front.elist[j]);
	  match_e_back[i].push_back(back.elist[i]);
	}
      }
      //then check for pairs of channels that match
      if(front.mult > 1){
	for (int n=0;n<front.mult;n++){
	  int j;
	  for (j=n+1;j<front.mult;j++){
	    // look for "out of sequence" channel -> end of cluster
	    if (front.chlist[j]!=front.chlist[n]+(j-n))
	      break;
	  }
	  //already covered the lone channel case
	  if(j==n+1) continue; 
	  
	  //loop through the cluster looking at pairs of channels
	  for(int k=n;k<j;k++){
	    double sume=0.;
	    double sumc=0.;
	    //calculate summed energy and center of mass for the pair
	    for(int m=k;(m<k+2)&&(m<j);m++){
	      sume += front.elist[m];
	      sumc += front.elist[m]*front.chlist[m];
	    }
	    //compare this pair to previous matches
	    //if its a better match push it onto the end
	    double t_eps = fabs(sume - back.elist[i]);
	    if(t_eps < eps[i].back()){
	      eps[i].push_back(t_eps);
	      match_ch[i].push_back(sumc/sume);
	      match_e[i].push_back(front.elist[j]);
	      match_e_back[i].push_back(back.elist[i]);
	    }
	  }
	}
      }
    }
    
    //Now look for conflicts
    //resolve them with lowest eps value
    for(int i=0;i<nback;i++){
      for(int j=0;j<nback;j++){
	if((j != i) && (match_ch[i].size() > 1) 
	   && (match_ch[j].size() > 1)){
	  double ch1 = match_ch[i].back();
	  double ch2 = match_ch[j].back();
	  if(1. > fabs(ch1 - ch2)){
	    //we have conflicting matches
	    //remove the worst one
	    if(eps[i].back() <= eps[j].back()){
	      eps[j].pop_back();
	      match_ch[j].pop_back();
	      match_e[j].pop_back();
	      match_e_back[j].pop_back();
	    }
	    else{
	      eps[i].pop_back();
	      match_ch[i].pop_back();
	      match_e[i].pop_back();
	      match_e_back[i].pop_back();
	    }
	  }
	}

      }
    }

    //now put the hits into the cluster list
    for(int i=0;i<nback;i++){
      if(match_ch[i].size() > 1){
	//calculate position of particle hit and put it into the cluster list
	double cf = match_ch[i].back();
	double cb = back.chlist[i];
	if((cf < static_cast<double>(num_front_)) 
	   && (cb < static_cast<double>(num_back_))){
	  Vector3 particlvect = chVect(cf,cb);
	  double xcrd = particlvect.X();
	  double ycrd = particlvect.Y();
	  double zcrd = particlvect.Z();
	  if (match_enefromback){
	    out.InsertHit(match_e_back[i].back(),back.tlist[i],
			  cf,cb,xcrd,ycrd,zcrd);
	  }
	  else{
	    out.InsertHit(match_e[i].back(),back.tlist[i],
			  cf,cb,xcrd,ycrd,zcrd);
	    
	  }
	}
      }
    }
  }
  return static_cast<int>(out.mult);

}

////////////////////////////////////////
//////// Beging WDetector Changes //////
////////////////////////////////////////


WDetectorPar::WDetectorPar(const std::string& basename,
			     const unsigned int& num_front,
			     const unsigned int& num_back,
			     const Vector3& posv):
  num_front_(num_front),
  num_back_(num_back),
  posv_(posv),
  shiftv_(0.,0.,0.),
  rotv_(0.,0.,0.),
  fmult(basename+".fmult",6),
  bmult(basename+".bmult",6),
  patter_e(basename+".patter_e",6),
  patter_t(basename+".patter_t",6),
  e(basename+".e",4096,ECAL_LOW,ECAL_HIGH,"MeV",num_front+num_back,0),
  t(basename+".t",4096,TCAL_LOW,TCAL_HIGH,"ns",num_front+num_back,0),
  efrontmatch(basename+".efrontmatch",4096,ECAL_LOW,ECAL_HIGH,"MeV"),
  ebackmatch(basename+".ebackmatch",4096,ECAL_LOW,ECAL_HIGH,"MeV"),
  efronttotal(basename+".efronttotal",4096,ECAL_LOW,ECAL_HIGH,"MeV"),
  ebacktotal(basename+".ebacktotal",4096,ECAL_LOW,ECAL_HIGH,"MeV"),
  frontmatchstat(basename+".frontmatchstat",256,0,256,""),
  backmatchstat(basename+".backmatchstat",256,0,256,""),
  f_total(basename+".f_total",num_front,1 /* sorted by channels*/ ),
  b_total(basename+".b_total",num_back,1 /* sorted by channels*/ )
{
  delta_x_ = (WIDTH_X) / static_cast<double>(num_front_);
  delta_y_ = (WIDTH_Y) / static_cast<double>(num_back_);
  match_epsilon=0;
  match_delta=0.1;
  match_maxene=4096;
  match_enefromback=0;
  addback_front=0;
  addback_back=0;
}

 int WDetectorPar::SetMatchParameters(int match_enefromback,
				       double match_epsilon,
				       double match_delta,
				       double match_maxene,
				       int addback_front,
				       int addback_back)
{
  this->match_enefromback=match_enefromback;
  this->match_epsilon=match_epsilon;
  this->match_delta=match_delta;
  this->match_maxene=match_maxene;
  this->addback_front=addback_front;
  this->addback_back=addback_back;
  return 0;
}


std::string WDetectorPar::print()
{
  std::string out;
  std::ostringstream tmp;
  for(int i=0;i<num_front_+num_back_;i++){
    if(e[i].isValid())
      tmp << e[i].getName() << " " << e[i] << std::endl;
    if(t[i].isValid())
      tmp << t[i].getName() << " " << t[i] << std::endl;
  }
  tmp << f_total.print() << b_total.print();
  out = tmp.str();
  return out;
}

void WDetectorPar::reset()
{
  fmult = 0;
  bmult = 0;
  for(int i=0;i<num_front_+num_back_;i++){
    /* Clear the TreeParameters which may not be set by Unpacker */
    e[i].Reset();
    t[i].Reset();
  }
  f_total.reset();
  b_total.reset();
  no_hits = 0;
}

Vector3 WDetectorPar::chVect(const double& cf,const double& cb){
  //First make a vector to the channel for a detector at the origin
  //Variables: wtheta = angle of det. relative to beam axis
  //           wphi = angle of det. relative to perpendicular of axis
  double x;
  double y;
  double s;
  double wtheta, wphi, d;
  //double dettheta, detphi, detzpos;
  //double phipos, thetapos, zpos;
  if((cf < static_cast<double>(num_front_)) 
     && (cb < static_cast<double>(num_back_))){
    x = - 25 + ((((double)rand())/RAND_MAX)+cb)*delta_x_;
    y = - 25 + ((((double)rand())/RAND_MAX)+cf)*delta_y_; //delta_y_ * (7.5 - (((double)rand())/RAND_MAX)-cf);
    s = sqrt(x*x + y*y);
    d = 660; //zpos (should match Resolut.cpp);
    wtheta = atan(s/d); //angle of chamber table;
    wphi = atan(y/x); //-85*3.14195/180;
    }
  else{ 
    std::cerr << "InValid channels sent to WDetectorPar::chVect" << std::endl;
    //    exit(EXIT_FAILURE);
    Vector3 resultv;
    return resultv;
  }

  Vector3 resultv;
  //  resultv.SetXYZ(x*cos((M_PI/180)*(wtheta))+d*cos((M_PI/180)*(90-wtheta)),y,d*(cos((M_PI/180)*(wtheta))));
  resultv.SetXYZ(x,y,0.0);
  //not sure on +/- of last component of x and z
  //resultv.SetXYZ(x*cos(wphi),y,x*sin(wphi));
  //resultv.SetXYZ(x,y,0);

  //then calibrate the vector 
  //add it to the detector position vector
  //if(rotv_[0]){
  //  resultv.RotateX(rotv_[0]*M_PI/180.);
  //}
  //if(rotv_[1]){
  //  resultv.RotateY(rotv_[1]*M_PI/180.);
  //}
  //if(rotv_[2]){
  //  resultv.RotateZ(rotv_[2]*M_PI/180.);
  //}
  resultv = resultv + shiftv_;
  resultv = resultv + posv_;

  return resultv;
}

int WDetectorPar::ReconstructClusters2(HitList& front,
					HitList& back,
					ClusterList& out) 
{ 
  TempList FrontClusters(16),BackClusters(16);
  int nback = static_cast<int>(back.mult);
  if(!nback) return 0;
  int i,j,front_match;
  unsigned int frontstatus=0;
  unsigned int backstatus=0;

  efronttotal=0.;
  ebacktotal=0.;
  // first add-back neighboring hits
  // in the back
  for (i=0;i<back.mult;i++){
    backstatus |= (((unsigned int)back.chlist[i]));
    ebacktotal = ebacktotal+back.elist[i];
    if (addback_back&&
	(i+1<back.mult)&&(back.chlist[i+1]==back.chlist[i]+1)){
      double e=back.elist[i+1]+back.elist[i];
      double t=(back.tlist[i+1]+back.tlist[i])/2;
      double ch=(back.chlist[i+1]+back.chlist[i])/2.0;
      BackClusters.InsertHit(e,t,ch);
      // we analyzed the i+1 hit along with the i, so 
      i++;
    }
    else{
      double e=back.elist[i];
      double t=back.tlist[i];
      double ch=back.chlist[i];
      BackClusters.InsertHit(e,t,ch);
    }
  }
  //in the front
  for (i=0;i<front.mult;i++){
    efronttotal = efronttotal+front.elist[i];
  }
  for (i=0;i<front.mult;i++){
    frontstatus |= (((unsigned int)front.chlist[i]));
    if (addback_front&&
	(i+1<front.mult)&&(front.chlist[i+1]==front.chlist[i]+1)){
      double e=front.elist[i+1]+front.elist[i];
      double t=(front.tlist[i+1]+front.tlist[i])/2;
      double ch=(front.chlist[i+1]+front.chlist[i])/2.0;
      FrontClusters.InsertHit(e,t,ch);
      // we analyzed the i+1 hit along with the i, so 
      i++;
    }
    else{
      double e=front.elist[i];
      double t=front.tlist[i];
      double ch=front.chlist[i];
      FrontClusters.InsertHit(e,t,ch);
    }
  }
  front_match=0;i=0;
  efrontmatch=0.;ebackmatch=0.;
  // Front and Backclusters are sorted by energy, high first
  while ((i<BackClusters.mult)&&(front_match<FrontClusters.mult)){
    double eps = FrontClusters.elist[front_match]-BackClusters.elist[i];
    double delta = match_epsilon+match_delta*BackClusters.elist[i];
    double back_e=BackClusters.elist[i];
    if (back_e>match_maxene || (fabs(eps)<=delta)){
      // we have a match
      double match_e= FrontClusters.elist[front_match];
      double match_ch= FrontClusters.chlist[front_match];
      double cb=BackClusters.chlist[i];
      double back_t=BackClusters.tlist[i];
      Vector3 particlvect = chVect(match_ch,cb);
      double xcrd = particlvect.X();
      double ycrd = particlvect.Y();
      double zcrd = particlvect.Z();
      if (match_enefromback){
	out.InsertHit(back_e,back_t,
		      match_ch,cb,xcrd,ycrd,zcrd);
      }
      else{
	out.InsertHit(match_e,back_t,
		      match_ch,cb,xcrd,ycrd,zcrd);
      }// end if "match_enefromback" 
      {	
	// sum up matched energies front and back
	efrontmatch += match_e;
	ebackmatch += back_e;
      }
    }// end found match
    else if (eps>delta){
      // front energy is larger than back energy+delta , 
      // so we won't ever match 
      // this one in the future, drop it, pick next front, keep the back
      frontstatus|=128; // flags "unmatched" 
      front_match++;
      continue;
    }
    else if (eps< (-delta)){
      // back energy is larger than front energy-delta , 
      // so we won't ever match 
      // this one in the future, drop it, pick next back, keep the front
      backstatus|=128; // flags "unmatched" 
      i++;
      continue;
    }

    front_match++;
    i++;
  }
  frontmatchstat=frontstatus;
  backmatchstat=backstatus;
  return(0);
}


int WDetectorPar::ReconstructClusters(HitList& front,
				       HitList& back,
				       ClusterList& out) 
{ 
  int nback = static_cast<int>(back.mult);
  double threshold=match_delta;
  if(!nback) return 0;
  
  //front energy and back energy 
  //must be less than eps apart for a match

  //first the simple case of only one back-side hit
  if(nback == 1){
    double eps = match_epsilon+match_delta*back.elist[0];
    double match_ch = -1;
    double match_e = 0.;
    double match_e_back = back.elist[0];
    //look for single channel matches
    for(int i=0;i<front.mult;i++){
      double t_eps = fabs(front.elist[i] - back.elist[0]);
      if(back.elist[0]>match_maxene || (t_eps < eps)){
	eps = t_eps;
	match_ch = front.chlist[i];
	match_e = front.elist[i];
      }
    }
    //now check if double channel has better match
    if(front.mult > 1){
      for (int i=0;i<front.mult;i++){
	int j;
	for (j=i+1;j<front.mult;j++){
	  // look for "out of sequence" channel -> end of cluster
	  if (front.chlist[j]!=front.chlist[i]+(j-i))
	    break;
	}
	//already covered the lone channel case
	if(j==i+1) continue; 
	
	//loop through the cluster looking at pairs of channels
	for(int k=i;k<j;k++){
	  double sume=0.;
	  double sumc=0.;
	  //calculate summed energy and center of mass for the pair
	  for(int m=k;(m<k+2)&&(m<j);m++){
	    sume += front.elist[m];
	    sumc += front.elist[m]*front.chlist[m];
	  }
	  //compare this pair to previous matches
	  double t_eps = fabs(sume - back.elist[0]);
	  if(t_eps < eps){
	    eps = t_eps;
	    match_ch = sumc/sume;
	    match_e = sume;
	  }
	}
      }
    }
    //if there was a match insert the hit
    if(match_ch >= 0.){
      //calculate position of particle hit and put it into the cluster list
      double cf = match_ch;
      double cb = back.chlist[0];
      if((cf < static_cast<double>(num_front_)) 
	   && (cb < static_cast<double>(num_back_))){
	Vector3 particlvect = chVect(cf,cb);
	double xcrd = particlvect.X();
	double ycrd = particlvect.Y();
	double zcrd = particlvect.Z();
	if (match_enefromback){
	  out.InsertHit(match_e_back,back.tlist[0],
			cf,cb,xcrd,ycrd,zcrd);
	}
	else{
	  out.InsertHit(match_e,back.tlist[0],
		      cf,cb,xcrd,ycrd,zcrd);
	}
      }
    }
  }
  else{ //more than one backside hit

    //declare storage containers
    std::vector<std::vector<double> > eps(nback);
    std::vector<double> erow(1,0.);
    std::vector<std::vector<double> > match_e(nback,erow);
    std::vector<std::vector<double> > match_e_back(nback,erow);
    std::vector<double> chrow(1,-1.);
    std::vector<std::vector<double> > match_ch(nback,chrow);

    //loop over the backside hits
    //building a list of acceptable matches for each
    for(int i=0;i<nback;i++){
      eps[i].push_back(match_epsilon+match_delta*back.elist[i]);
      //first look for single channel matches
      for(int j=0;j<front.mult;j++){
	double t_eps = fabs(front.elist[j] - back.elist[i]);
	if(t_eps < eps[i].back()){
	  eps[i].push_back(t_eps);
	  match_ch[i].push_back(front.chlist[j]);
	  match_e[i].push_back(front.elist[j]);
	  match_e_back[i].push_back(back.elist[i]);
	}
      }
      //then check for pairs of channels that match
      if(front.mult > 1){
	for (int n=0;n<front.mult;n++){
	  int j;
	  for (j=n+1;j<front.mult;j++){
	    // look for "out of sequence" channel -> end of cluster
	    if (front.chlist[j]!=front.chlist[n]+(j-n))
	      break;
	  }
	  //already covered the lone channel case
	  if(j==n+1) continue; 
	  
	  //loop through the cluster looking at pairs of channels
	  for(int k=n;k<j;k++){
	    double sume=0.;
	    double sumc=0.;
	    //calculate summed energy and center of mass for the pair
	    for(int m=k;(m<k+2)&&(m<j);m++){
	      sume += front.elist[m];
	      sumc += front.elist[m]*front.chlist[m];
	    }
	    //compare this pair to previous matches
	    //if its a better match push it onto the end
	    double t_eps = fabs(sume - back.elist[i]);
	    if(t_eps < eps[i].back()){
	      eps[i].push_back(t_eps);
	      match_ch[i].push_back(sumc/sume);
	      match_e[i].push_back(front.elist[j]);
	      match_e_back[i].push_back(back.elist[i]);
	    }
	  }
	}
      }
    }
    
    //Now look for conflicts
    //resolve them with lowest eps value
    for(int i=0;i<nback;i++){
      for(int j=0;j<nback;j++){
	if((j != i) && (match_ch[i].size() > 1) 
	   && (match_ch[j].size() > 1)){
	  double ch1 = match_ch[i].back();
	  double ch2 = match_ch[j].back();
	  if(1. > fabs(ch1 - ch2)){
	    //we have conflicting matches
	    //remove the worst one
	    if(eps[i].back() <= eps[j].back()){
	      eps[j].pop_back();
	      match_ch[j].pop_back();
	      match_e[j].pop_back();
	      match_e_back[j].pop_back();
	    }
	    else{
	      eps[i].pop_back();
	      match_ch[i].pop_back();
	      match_e[i].pop_back();
	      match_e_back[i].pop_back();
	    }
	  }
	}

      }
    }

    //now put the hits into the cluster list
    for(int i=0;i<nback;i++){
      if(match_ch[i].size() > 1){
	//calculate position of particle hit and put it into the cluster list
	double cf = match_ch[i].back();
	double cb = back.chlist[i];
	if((cf < static_cast<double>(num_front_)) 
	   && (cb < static_cast<double>(num_back_))){
	  Vector3 particlvect = chVect(cf,cb);
	  double xcrd = particlvect.X();
	  double ycrd = particlvect.Y();
	  double zcrd = particlvect.Z();
	  if (match_enefromback){
	    out.InsertHit(match_e_back[i].back(),back.tlist[i],
			  cf,cb,xcrd,ycrd,zcrd);
	  }
	  else{
	    out.InsertHit(match_e[i].back(),back.tlist[i],
			  cf,cb,xcrd,ycrd,zcrd);
	    
	  }
	}
      }
    }
  }
  return static_cast<int>(out.mult);
}

////////////////////////////////////////
///////// END WDETECTOR CHANGES ////////
////////////////////////////////////////

ICDetectorPar::ICDetectorPar(std::string name,unsigned int no_channels):
  e(name+".e",4096,ECAL_LOW,ECAL_HIGH,"MeV",no_channels,0),
  t(name+".t",4096,TCAL_LOW,TCAL_HIGH,"ns",no_channels,0),
  patter_e(name+".patter_e",6),
  patter_t(name+".patter_t",6),
  mult(name+".mult",6),
  total(name+".total",no_channels,0 /* sorted by energy */)
{
}

std::string ICDetectorPar::print()
{
  std::string out;
  std::ostringstream tmp;
  for(int i=0;i<e.size();i++){
    if(e[i].isValid())
      tmp << e[i].getName() << " " << e[i] << std::endl;
    if(t[i].isValid())
      tmp << t[i].getName() << " " << t[i] << std::endl;
  }
  tmp << total.print();
  out = tmp.str();
  return out;
}

void ICDetectorPar::reset()
{
  mult = 0;
  total.reset();
}

SIDetectorPar::SIDetectorPar(std::string name,unsigned int no_channels):
  e(name+".e",4096,ECAL_LOW,ECAL_HIGH,"MeV",no_channels,0),
  t(name+".t",4096,TCAL_LOW,TCAL_HIGH,"ns",no_channels,0),
  x(name+".x",4096,0,4096,"mm"),
  y(name+".y",4096,0,4096,"mm"),
  patter_e(name+".patter_e",6),
  patter_t(name+".patter_t",6),
  mult(name+".mult",6),
//  relt(name+".relt",4096,0.,220.,"ns"),
  total(name+".total",no_channels,0 /* sorted by energy*/ )
{
}

std::string SIDetectorPar::print()
{
  std::string out;
  std::ostringstream tmp;
  for(int i=0;i<e.size();i++){
    if(e[i].isValid())
      tmp << e[i].getName() << " " << e[i] << std::endl;
    if(t[i].isValid())
      tmp << t[i].getName() << " " << t[i] << std::endl;
  }
  if(x.isValid())
    tmp << x.getName() << " " << x << std::endl;
  if(x.isValid())
    tmp << y.getName() << " " << y << std::endl;
  out = tmp.str();
  tmp << total.print();
  out = tmp.str();
  return out;
}

void SIDetectorPar::reset()
{
  mult = 0;
  total.reset();
}

MCPDetectorPar::MCPDetectorPar(std::string name,unsigned int no_channels):
  e(name+".e",4096,ECAL_LOW,ECAL_HIGH,"MeV",no_channels,0),
  t(name+".t",4096,TCAL_LOW,TCAL_HIGH,"ns",no_channels,0),
  total(name+".total",4096,0,4096,"MeV"),
  x(name+".x",4096,0,4096,"mm"),
  y(name+".y",4096,0,4096,"mm"),
  patter_e(name+".patter_e",6),
  patter_t(name+".patter_t",6),
  mult(name+".mult",6)
{
}

std::string MCPDetectorPar::print()
{
  std::string out;
  std::ostringstream tmp;
  for(int i=0;i<e.size();i++){
    if(e[i].isValid())
      tmp << e[i].getName() << " " << e[i] << std::endl;
    if(t[i].isValid())
      tmp << t[i].getName() << " " << t[i] << std::endl;
  }
  tmp << total.getName() << " " << total << std::endl;
  if(x.isValid())
    tmp << x.getName() << " " << x << std::endl;
  if(x.isValid())
    tmp << y.getName() << " " << y << std::endl;
  out = tmp.str();
  out = tmp.str();
  return out;
}

void MCPDetectorPar::reset()
{
  mult = 0;
}

ResneutDetectorPar::ResneutDetectorPar(std::string name,unsigned int no_channels):
  e(name+".e",4096,ECAL_LOW,ECAL_HIGH,"MeV",no_channels,0),
  t(name+".t",4096,TCAL_LOW,TCAL_HIGH,"ns",no_channels,0),
  total(name+".total",4096,0,4096,"MeV"),
  x(name+".x",4096,0,4096,"mm"),
  y(name+".y",4096,0,4096,"mm"),
  patter_e(name+".patter_e",6),
  patter_t(name+".patter_t",6),
  mult(name+".mult",6)
{
}

std::string ResneutDetectorPar::print()
{
  std::string out;
  std::ostringstream tmp;
  for(int i=0;i<e.size();i++){
    if(e[i].isValid())
      tmp << e[i].getName() << " " << e[i] << std::endl;
    if(t[i].isValid())
      tmp << t[i].getName() << " " << t[i] << std::endl;
  }
  tmp << total.getName() << " " << total << std::endl;
  if(x.isValid())
    tmp << x.getName() << " " << x << std::endl;
  if(x.isValid())
    tmp << y.getName() << " " << y << std::endl;
  out = tmp.str();
  out = tmp.str();
  return out;
}

void ResneutDetectorPar::reset()
{
  mult = 0;
}


RawPar::RawPar(std::string name,
	       unsigned int si_a_chan,
	       unsigned int si_b_chan,
	       unsigned int si_c_chan,
	       unsigned int si_d_chan,
	       unsigned int si_x_chan,
	       unsigned int s2_stop_chan,
	       unsigned int ion_chan,
	       unsigned int rftime_chan):
  si_a("sia",si_a_chan),
  si_b("sib",si_b_chan),
  si_c("sic",si_c_chan),
  si_d("sid",si_d_chan),
  si_x("six",si_x_chan),
  s2_stop("s2_stop",s2_stop_chan),
  ion("ion",ion_chan),
  rftime("rftime",rftime_chan),
  mcp1("mcp1",5),
  mcp2("mcp2",5),
  resneut1("resneut1",3),
  resneut2("resneut2",3),
  resneut3("resneut3",3),
  resneut4("resneut4",3),
  ic_stop("ic_stop",1),
  triggerbits("trigger.bits",256,0,256,"bits")
{}

std::string RawPar::print()
{
  std::string out;
  std::ostringstream tmp;
  tmp << si_a.print() << si_b.print() << si_c.print() << si_d.print() ;
  tmp << ion.print() << rftime.print() << s2_stop.print() << mcp1.print() << mcp2.print();
  tmp << ic_stop.print();
  out = tmp.str();
  return out;
}

void RawPar::reset()
{
  si_a.reset();
  si_b.reset();
  si_c.reset();
  si_d.reset();
  si_x.reset();
  s2_stop.reset();
  ion.reset();
  rftime.reset();
  mcp1.reset();
  mcp2.reset();
  resneut1.reset();
  resneut2.reset();
  resneut3.reset();
  resneut4.reset();
  ic_stop.reset();
}

CalPar::CalPar(std::string name,
	       unsigned int si_a_chanf,unsigned int si_a_chanb,Vector3 aposv,
	       unsigned int si_b_chanf,unsigned int si_b_chanb,Vector3 bposv,
	       unsigned int si_c_chanf,unsigned int si_c_chanb,Vector3 cposv,
	       unsigned int si_d_chanf,unsigned int si_d_chanb,Vector3 dposv,
	       unsigned int si_x_chanf,unsigned int si_x_chanb,Vector3 xposv,
	       unsigned int s2_stop_chan,unsigned int ion_chan,unsigned int rftime_chan):
  si_a(name+".si_a",si_a_chanf,si_a_chanb,aposv),
  si_b(name+".si_b",si_b_chanf,si_b_chanb,bposv),
  si_c(name+".si_c",si_c_chanf,si_c_chanb,cposv),
  si_d(name+".si_d",si_d_chanf,si_d_chanb,dposv),
  si_x(name+".si_x",si_x_chanf,si_x_chanb,xposv),
  s2_stop(name+".s2_stop",s2_stop_chan),
  ion(name+".ion",ion_chan),
  rftime(name+".rftime",rftime_chan),
  mcp1(name+".mcp1",5),
  mcp2(name+".mcp2",5),
  resneut1(name+".resneut1",3),
  resneut2(name+".resneut2",3),
  resneut3(name+".resneut3",3),
  resneut4(name+".resneut4",3),
  ic_stop(name+".ic_stop",1),
  clusters_a(name+".clusters_a",si_a_chanf+si_a_chanb),
  clusters_b(name+".clusters_b",si_b_chanf+si_b_chanb),
  clusters_c(name+".clusters_c",si_c_chanf+si_c_chanb),
  clusters_d(name+".clusters_d",si_d_chanf+si_d_chanb),
  clusters_x(name+".clusters_x",si_x_chanf+si_x_chanb),
  esum(name+".esum",4096,ECAL_LOW,3*ECAL_HIGH,"MeV"),
  abc_bmult(name+".abc_bmult",6),
  ab_bmult(name+".ab_bmult",6),
  ac_bmult(name+".ac_bmult",6),
  bc_bmult(name+".bc_bmult",6),
  a_unmatch(name+".a_unmatch",6),
  b_unmatch(name+".b_unmatch",6),
  c_unmatch(name+".c_unmatch",6),
  d_unmatch(name+".d_unmatch",6)
{
  // set non-default restructions for si_c:
  // energy from back, 2 MeV+2%E match, addback on front only 
  si_c.SetMatchParameters(1,4096,1,10,0,0);
  si_a.SetMatchParameters(1,0.0,0.1,1024,0,0);
  si_b.SetMatchParameters(1,0.0,0.1,1024,0,0);
  si_d.SetMatchParameters(1,0.0,0.1,1024,0,0);
  si_x.SetMatchParameters(1,0.0,0.1,1024,0,0);
}

void CalPar::reset()
{
  si_a.reset();
  si_b.reset();
  si_c.reset();
  si_d.reset();
  si_x.reset();
  s2_stop.reset();
  ion.reset();
  rftime.reset();
  mcp1.reset();
  mcp2.reset();
  resneut1.reset();
  resneut2.reset();
  resneut3.reset();
  resneut4.reset();
  
  clusters_a.reset();
  clusters_b.reset();
  clusters_c.reset();
  clusters_d.reset();
  clusters_x.reset();

  abc_bmult = 0;
  ab_bmult = 0;
  ac_bmult = 0;
  bc_bmult = 0;
}

std::string CalPar::print()
{
  std::string out;
  std::ostringstream tmp;
  tmp << si_a.print() << si_b.print() << si_c.print() << si_d.print()<< s2_stop.print();
  tmp << ion.print() << rftime.print();
  tmp << clusters_a.print() << clusters_b.print() << clusters_c.print()<< clusters_d.print();
  out = tmp.str();
  return out;
}
