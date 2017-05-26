// create calibration file with trivial parameters
#include <fstream>

void X3FinalFix_ini(void) {
  using namespace std;

  ofstream outfile;
  outfile.open("X3FinalFix_init.dat");
  outfile<<"DetNum\tChNum\tOffset\n";

  for (Int_t DetNum=4; DetNum<28; DetNum++) {
    Double_t offset=0;//default, back
    for (Int_t ChNum=0; ChNum<12; ChNum++) {
      if(ChNum>7)//Step1:up
	offset=0;
      if(ChNum>3&&ChNum<8)
	offset=0;//Step2:down
      outfile << DetNum << "\t" << ChNum << "\t" << offset << endl;
    }
  }
  outfile.close();
}
