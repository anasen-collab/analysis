# Super-X3 Position Calibration
Relative position calibration of the SX3 silicon detectors.
The SX3 detectors correspond to detectors number 4-27.
Each detector has 4 front channels which are position sensitive and 4 back channels.
The detector strips on the front and back are oriented in a perpendicular fashion, forming 16 pseudo pixels.
Each of the 384 pseudo pixels is calibrated using the known position range of each pseudo pixel.

## Instructions
### General Usage
Use the following command to run the program `root -l GeomtryCal.C+`
There are a few important setting which mist be considered.
* First, the max bin content ratio used to determine the edge of the detectors.
The default value is 3.0 and seems to work well.
* Re-binning the histograms should be used as necessary, depending on statistics.
* The starting bin of the loop should be considered to skip junk bins, such as z=-1.0.

### Output files (.dat files)
The output file (e.g.`X3geometry.dat`) has the following columns:
Detector number, Front channel number, Back channel number, upstream edge, downstream edge.
The first line of `.dat` file is a header line.
Each of the 24 SX3 detectors has 16 pseudo pixels.
Each calibration file should have 385 lines, including the header line.

The `.dat` files  included in the repository are an example. The run-to-run changes in the `.dat` files are excuded by `.gitignore`. To force the updated files to be saved to the repository, use the command `git add -f file.dat`.

For reference, trivial calibration files can be generated with the file `X3geometry_ini.C` in the [saves](saves) directory. These programs are  run using the command `root -l X3geometry_ini.C`

### Input files
Data with high statistics is needed. Each pseudo pixel (1/16th of each detector) needs to have sufficient statistics to locate the edges of the strips.

### Histograms
enable the generation of the 
`TString hname=Form("SX3ZposCal_%i_%i_%i",DetNum,FrontChNum,BackChNum);`
histograms using `#define ZPosCal` in Main.ccp.

### Gains
The position on the detector is calculated in Silicon_Cluster.h.

The calculated parameters are applied by the `ChannelMap` method `PosCal`
````
FinalZPosCal = (EdgeDCal-EdgeUCal)/(EdgeDown[DNum-4][StripNum][BChNum]-EdgeUp[DNum-4][StripNum][BChNum])
    *(FinalZPos-EdgeDown[DNum-4][StripNum][BChNum])+EdgeDCal;
````
Here EdgeDCal and EdgeUCal are the known locations of the strip edges and EdgeDown and  EdgeUp are the edges measured in GeomtryCal.C.
FinalZPos is defined over the range -1 to +1.
The variable FinalZPosCal is defined over the range 0 to 7.5cm and gives the physical position in cm on the silicon strip.
### Next steps

## Fiting Method
a half max threshold is set in the code. A typical value is the maximum peak height divided by 3.
