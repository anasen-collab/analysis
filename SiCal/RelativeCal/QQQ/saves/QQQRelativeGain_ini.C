// create calibration file with trivial parameters
#include <fstream>

void QQQRelativeGain_ini(void) {
  using namespace std;

  ofstream outfile;
  outfile.open("QQQRelativeGains_ini.dat");
  outfile<<"DetNum\tChNum\tGain\n";

  for (Int_t DetNum=0; DetNum<4; DetNum++) {
    Double_t gain=1;
    for (Int_t ChNum=0; ChNum<32; ChNum++) {
      if(ChNum>16)
	gain=1;//front
      outfile << DetNum << "\t" << ChNum << "\t" << gain << endl;
    }
  }
  outfile.close();
}
