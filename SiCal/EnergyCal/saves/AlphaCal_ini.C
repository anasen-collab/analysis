// create calibration file with trivial parameters
#include <fstream>
void AlphaCal_ini (void) {
  Double_t gain = 1;
  ofstream outfile; 
  outfile.open("AlphaCal_init.dat");
  outfile << "DetNum\tDummy\tGain\n";
  for (Int_t DetNum=0; DetNum<28; DetNum++) {
    if(DetNum>3)
      gain=1;//R1
    if(DetNum>15)
      gain=1;//R2
    outfile << DetNum << "\t"<< -1 << "\t" << gain <<endl;
  }
}
