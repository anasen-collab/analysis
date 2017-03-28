# Main.cpp
 This file creates a tree with vectors, objects and members
## Execution
To do so, we have to create a dictionary
 To create a dictionary:
`rootcint -f Main_dict.cxx -c tree_structure.h LinkDef.h`
 To compile the code:
`g++ -o Main Main_dict.cxx Main.cpp `root-config --cflags --glibs` -O3`
These two commands can be run using the makefile.
`make`
To run the code:
 `./Main input.root output.root`
  for example:
  `./Main ../../../../../data0/manasta/evt2root_files_new/24Mg/runXXX.root ../../../../../data0/manasta/OrganizeRaw_files/24Mg/runXXX.root`

## Options
To select the component of the beam, mostly for Radio-active beams, disable while you work with Calibration data & enable while you do data analysis

`#define MCP_RF_Cut`