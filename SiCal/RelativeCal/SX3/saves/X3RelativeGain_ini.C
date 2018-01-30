// create calibration file with trivial parameters
#include <fstream>

void X3RelativeGain_ini(void) {
  using namespace std;

  ofstream outfile;
  outfile.open("X3RelativeGains_init.dat");
  outfile<<"DetNum\tChNum\tGain\n";

  for (Int_t DetNum=4; DetNum<28; DetNum++) {
    Double_t gain=1;//default, back
    for (Int_t ChNum=0; ChNum<12; ChNum++) {
      if(ChNum>7)//Step1:up
	gain=1;
      if(ChNum>3&&ChNum<8)
	gain=1;//Step2:down
      outfile << DetNum << "\t" << ChNum << "\t" << gain << endl;
    }
  }
  outfile.close();
}
