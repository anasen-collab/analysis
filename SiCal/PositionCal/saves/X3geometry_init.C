// create calibration file with trivial parameters
#include <fstream>

void X3geometry_init(void) {
  using namespace std;

  ofstream outfile;
  outfile.open("X3geometry_init.dat");
  outfile<<"DetNum\tFrontCh\tBackCh\tmin\tmax\n";

  Double_t up=0;//default, back
  Double_t down=0;//default, back
  
  for (Int_t DetNum=4; DetNum<28; DetNum++) {
    for (Int_t FrNum=0; FrNum<4; FrNum++) {
      for (Int_t BkNum=0; BkNum<4; BkNum++) {
	up=-1+0.5*BkNum;
	down=up+0.5;
	outfile << DetNum << "\t" << FrNum << "\t" << BkNum << "\t" << up << "\t" << down << endl;
    }
  }
  }
  outfile.close();
}
