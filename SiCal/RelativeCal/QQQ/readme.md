# QQQ Relative Calibration
Relative calibration of Si gains for QQQ. Essentially the same progam as that for the SX3 (or vice-versa).
Edited by : John Parker , 2016Jan22
Developed by : Jon Lighthall, 2016.12
## General Usage
The data for the QQQ relaive gains calibration should have high statistics over a range of energies; such as no-target in-gas runs.
### Instructions
1. First, run `Main.cpp` with a given `.dat` file all 1's
   * Edit `Main.cpp` to initialize the channel map with a `.dat` file with all 1's
   * Compile and run `Main.cpp` to generate a new `.root` file
   * Input the new `.root` and the `.dat` file into this code. Make sure you are using the same `.dat` file that was reference by `Main.cpp`.
   * Run `root -l SiRelativeGains_Step1.C+`
   * This code will output a new `.dat` file with updated gains
2. Second, run Main.cpp again with the new .dat file
   * Edit `Main.cpp` to initialize the channel map with the new `.dat` file
   * Compile and run `Main.cpp` to generate a new `.root` file
   * Input the new `.root` and the `.dat` file from Step 1 into this code. Make sure you are using the same `.dat` file that was reference by `Main.cpp`.
   * Run  `root -l SiRelativeGains_Step2.C+`

### Output files (`.dat` files)
Output file (e.g.`QQQRelativeGains_Step1.dat`) has the following columns:
Detector number, Front channel, Slope
The first line of dat files is a dummy line.
The QQQ detectors correspond to detectors number 0-3. Each detector has 32 channels (0-31).

## Step 1
Loop over front channels. Channels 16-31 in the `.dat` file will be filled in. For each detector (`DetNum=0-3`), each front channel (`FrontChNum=0-15`) is gain-matched with a particular back channel (typically `BackChNum=0`).

### Methods 
1. Method 1 - calculates slope of points wihtin pre-defined cut using TGraph
2. Method 2 - calculates slope of points wihtin user-defined cut using TGraph
   This method works very similar to method 1, except that instead of a predefined cut, the user must
   draw their own graphical cut. To do this, in the canvas: View->Toolbar and click on the scissors on the top
   right hand corner. Then, select the region around your data by clicking. Double click to close the cut.
   A best fit line should appear through your data
3. Expanded Method 1 - three slopes are calculated:		
   1) slope from un-gated fit of 2D histogram		
   2) slope from gated fit of TGraph (same as Method 1)		
   3) slope from un-gated fit of profile of 2D histogram (uses the average y-position for each x-position)
4. Automated cut generation based on TProfile slope. A parrallelpiped cut is generated.
5. Expansion of Method 4. Used for Det 2. Uses upper and lower limits for the lower bound of the gate; i.e., the gate "backs in" towards the origin over the steps.
6. Automated cut generation based on TProfile slope; cone-shaped cut. Calculates slopes within 0.2% of those of Mehtod 4 with the advantage of beign able to better detect the slope for detectors with multiple loci.

## Step 2
Loop over back channels. Channels 0-15 in the `.dat` file will be filled in. For each detector (`DetNum=0-3`), each back channel (`BackChNum=1-15`) is gain-matched with a particular back channel (typically `FrontChNum=0`).
