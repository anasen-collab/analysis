# Analyzer
// Goal: To Analyze the Tracking and Reconstruction of the 7Be+d-> p + alpha + alpha Events  
// & to Analyze the Elastic Scattering of Deuterons to Calibrate ANASEN for 7Be+d Experiments..
// & for the other (d,p),(d,alpha)..etc..ANASEN experiments with Gas volume target
//
// Uses Lookup tables instead of doing integration multiple times for Energyloss, 
// Final Energy, Initial Energy & Distance calculation.

## Versions

* Analyzer_Cal
* Analyzer_PCHeavyHit
  Included is the section of the code to check/plot for the heavy hit in the PC. This part will be  relevant when  the overflow is in your data, so make your Max PC Threshold be 4097 or higher. The plot puts the tracked signal in the 12th bin. Let me know if any questions!
* Analyzer_ES
  Here is the code for the elastic scattering from gas that basically reconstructs the beam energy (7Be in this case) from the elastically scattered gas (deuteron) detected in the Silicon. You need to change your beam and gas accordingly. Let me know if there are any questions.
  Goal: To Analyze the Elastic Scattering of Deuterons to Calibrate ANASEN for 7Be+d Experiments..
 & for the other (d,p),(d,alpha)..etc..ANASEN experiments with Gas volume target
## Input

the file has 3 inputs and one output

the raw file, the calibated file and a cut file.

the second input, `DataList.txt` gives the location of a root file. 


## Usage
To create a dictionary:
````
rootcint -f tr_dict.cxx -c tree_structure.h LinkDef.h
rootcint -f tr_dict.cxx -c ../Include/tree_structure.h LinkDef.h
````
Usage: 

````
g++ -o Analyzer tr_dict.cxx LookUp.cpp Analyzer.cpp `root-config --cflags --glibs`
g++ -o Analyzer_ES tr_dict.cxx LookUp.cpp Analyzer_ES.cpp `root-config --cflags --glibs`
````

excecution
````
./Analyzer DataListCal.txt 282_3_4Cal5Analyzer20161102.root cut/He4.root 
./Analyzer_ES DataListCal.txt 2430Cal5Analyzer20170303.root cut/D2.root 
````
