// create calibration file with trivial parameters
#include <fstream>
void SiPulser_ini (void) {
  Double_t zeroshift = 0;
  Double_t vperch = 1;
  ofstream outfile; 
  outfile.open("Sipulser_ini.dat");
  outfile << "MBID\tCBID\tChan\tOffset\tVolts/Chan"<<endl;
  for (Int_t MBID=1; MBID<3; MBID++) {  
    for (Int_t CBID=1; CBID<15; CBID++) { 
      for (Int_t ChNum=0; ChNum<16; ChNum++) {
	if(MBID==2&&CBID>12)
	  continue;
	outfile << MBID << "\t" << CBID << "\t" << ChNum << "\t"<< zeroshift << "\t" << vperch <<endl;
      }
    }
  }
}
