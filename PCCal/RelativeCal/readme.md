#PC Relative calibration
The relative calibration for the PC is accomplished in two steps.
This needs to be done right after the PC pulser calibration.
First, the Up and Down signals are gain matched to each other.

Second, the sum of the gain-matched up and down signals, ie the energy, is gain matched.
Mutual gain matching is performed in lieu of energy calibration.

## Instructions
First you first have to create 2D histos in the main for up vs down in Main

run on file processed through Main.

 Usage: `root -l PC_UpDown_RelCal.C`

then use the script to calculate the slope and offset. 
A histogram for each of the 24 wires is plotted. 
use the cut tool to trace over the line

## Input file
Use a negative polarity pulser run to gain match the up and down signals.
By definition, after pulser calibration and gain matching, all pulser signals will be mutually gain matched.

## Output files
The output file has the following collumns
`Wire	Slope	Shift`
Each slope should be nearly equal to one and each offset should be nearly equal to zero.
The ouput file should have 25 lines, including the header line.

## Next steps
Then load the file with slope and offset in the main (add these changes files  to your Main.cpp & ChannelMap.h).
