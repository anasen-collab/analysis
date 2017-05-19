# Analyzer

## Versions

* Analyzer_Cal
* Analyzer_PCHeavyHit
  Included is the section of the code to check/plot for the heavy hit in the PC. This part will be  relevant when  the overflow is in your data, so make your Max PC Threshold be 4097 or higher. The plot puts the tracked signal in the 12th bin. Let me know if any questions!
* Analyzer_ES
  Here is the code for the elastic scattering from gas that basically reconstructs the beam energy (7Be in this case) from the elastically scattered gas (deuteron) detected in the Silicon. You need to change your beam and gas accordingly. Let me know if there are any questions.
## Input

the file has 3 inputs and one output

the raw file, the calibated file and a cut file.

the second input, `DataList.txt` gives the location of a root file. 
