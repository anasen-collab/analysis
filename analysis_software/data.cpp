/* Program: data.cpp
 * Description: Runs Main on the files listed in runs.lst.
 * See readme.md for general instructions.
 * Developed by J. Lighthall Oct 2017
 */

//C and C++ libraries
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cstdlib>

using namespace std;

int main() {
  int list[99];
  const char* fname="runs.lst"; //name of file with list of run numbers
  
  ifstream infile(fname);
  printf("Reading list file \"%s\"\n",fname);
  printf(" The following runs will be processed\n");
  int i=0;
  while (infile >> list[i]) {
    printf("  %d\n",list[i]);
    i++;
  }

  char ans;
  bool del=false;
  cout << " Do you want to delete the raw .root files after conversion? (y/n)" << endl << "  ";
  cin >> ans;
  if (ans=='y')
    del=true;
      
  for(int j=0;j<i;j++) {
    string str0 = "./Main "; //command
    string str1 = "/data0/lighthall/root/raw/"; //location of raw .root files

    string str;
    stringstream num;
    num << list[j];
    str=num.str();

    string str2 = "run";//file name
    str2+=str;
    str2+=".root ";
    
    string str3 = "/data0/lighthall/root/main/"; //output directory
    
    string str4 =  str0+str1+str2+str3+str2;
    //cout << str4 <<endl;
    system(str4.c_str()); 

    if(del) {
      str4= "rm -v ";
      str4+=str1+str2;
      system(str4.c_str()); //use this line to delete the files after conversion
    }
  }
}
