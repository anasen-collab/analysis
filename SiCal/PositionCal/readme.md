# Super-X3 Position Calibration
Relative position calibration of the SX3 silicon detectors.
The SX3 detectors correspond to detectors number 4-27.
Each detector has 4 front channels which are position sensitive and 4 back channels.
The detector strips on the front and back are oriented in a perpendicular fashion, forming 16 pseudo pixels.
Each of the 384 pseudo pixels is calibrated using the known position range of each pseudo pixel.

## Instructions
### General Usage
`root -l GeomtryCal.C+`


### Output files (`.dat` files)
The output file (e.g.`X3geometry.dat`) has the following columns:
Detector number, Front channel number, Back channel number, upstream edge, downstream edge.
The first line of `.dat` file is a header line.
Each of the 24 SX3 detectors has 16 pseudo pixels.
Each calibration file should have 385 lines, including the header line.

The `.dat` files  included in the repository are an example. The run-to-run changes in the `.dat` files are excuded by `.gitignore`. To force the updated files to be saved to the repository, use the command `git add -f file.dat`.

For reference, trivial calibration files can be generated with the file `X3geometry_ini.C` in the [saves](saves) directory. These programs are  run using the command `root -l X3geometry_ini.C`

### Input files

### Histograms
`TString hname=Form("SX3ZposCal_%i_%i_%i",DetNum,FrontChNum,BackChNum);`
turn on using `#define ZPosCal` in Main.ccp.

### Gains
### Next steps
## Fiting Methods 
