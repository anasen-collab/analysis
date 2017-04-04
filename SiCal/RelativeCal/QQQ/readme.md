# QQQ Relative Calibration
Relative calibration of Si gains
////Relative calibration of Si gains for QQQ
////Essentially the same progam as that for the SX3
//// Edited by : John Parker , 2016Jan22
//   Developed by : Jon Lighthall, 2016.12
## General Usage
### Instructions
////root -l SiRelativeGains_Step1.C+
1. First, run `Main.cpp` with a given `.dat` file all 1's
   * Input `organize.root` and the `.dat` file into this code
   * This code will output another .dat file that you should define
   * `root -l SiRelativeGains_Step1.C++`
2. Second, run Main.cpp again with the new .dat file
   * Input the new organize.root and new .dat file into Step2.C
   * `root -l SiRelativeGains_Step2.C++`

### .dat files
Output file (e.g.`_Slope1.dat`) has the following columns:
Detector number, Front channel, Slope
The first line of dat files is a dummy line.
The QQQ detectors correspond to detectors number 0-3. Each detector has 32 channels (0-31).

## Step 1
Loop over front channels. channels 16-31 in the `.dat` file will be filled in.

### Methods 
1. //Method 1 - calculates slope of points wihtin pre-defined cut using TGraph
2. //Method 2 - calculates slope of points wihtin user-defined cut using TGraph
  // This method works very similar to method 1, except that instead of a predefined cut, the user must
  // draw their own graphical cut. To do this, in the canvas: View->Toolbar and click on the scissors on the top
  // right hand corner. Then, select the region around your data by clicking. Double click to close the cut.
  // A best fit line should appear through your data
3.
4. //Method 4 - automated cut generation based on TProfile slope
5.
6.   //Method 6 - automated cut generation based on TProfile slope; cone-shaped

## Step 2
Loop over back channels
channels 0-15 in the `.dat` file will be filled in.