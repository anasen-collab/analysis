script to do the PC UP vs Down calibration. This needs to be done right after the PC pulser calibration. I have only attached the changes in the files.

## Instructions
 
First you first have to create 2D histos in the main for up vs down in Main

run on file processed through Main.

 Usage: `root -l PC_UpDown_RelCal.C`

then use the script to calculate the slope and offset. 
A histogram for each of the 24 wires is plotted. 
use the cut tool to trace over the line

Then load the file with slope and offset in the main (add these changes files  to your Main.cpp & ChannelMap.h).


