# Main.cpp
 This file creates a tree with vectors, objects and members
## Compiling 
Before compiling, the dictionary is created with the follwoing command:
```rootcint -f Main_dict.cxx -c tree_structure.h LinkDef.h```
 To code is then compiled with the command:
```g++ -o Main Main_dict.cxx Main.cpp `root-config --cflags --glibs` -O3```
These two commands can be run using the makefile.
```make```
## Execution
To run the code: 
```
./Main input.root output.root
```

## Files 
auto-generated files can be removed using the command `make clean`. 

The auto-generated files `Main_dict.cxx` and `Main_dict.h` may be included in the repository as an example. These filese, in addition to the executable file `Main`, are
excuded by `.gitignore`. The force a save of updated files use the command `git add -f foo.bar`.

## Options
* Beam diagnostics
   * `#define MCP_RF_Cut` To select the component of the beam, mostly for radio-active beams, disable while you work with calibration data & enable while you do data analysis
   * `#define IC_hists`
* Calibration 
   * `#define Pulser_ReRun`redefine this for cal
   * `#define Hist_for_Cal` 
   * `#define Hist_after_Cal` Select the Histograms for Calibration or for a Check.
   * `#define ZPosCal `

## ROOT
After compiling, the output `.root` files may be viewed in root. Doing so will yield class warnings unless the folling line is added to your `rootlogon.C` file.

````
gROOT->ProcessLine(".L /home/<user>/anasen/analysis_software/tree_structure.h");
````
