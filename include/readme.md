# Shared files ('/include')

This directory contains the header files which are shared by programs in different stages of the conversion, calibration, and analysis of ANSEN data.

## detclass.h
This file defines classes for detector readout.
### Used by
* evt2root_NSCL11.C
* Organize.cpp (deprecated)
* Main.cpp
## tree_structure.h
Defines the classes 'SiHit', 'PCHit', and 'Track' for use in calibration and analysis. This file has been developed from 'treestructure.h' and 'organizetree.h'.
### Used by
* Main.cpp
* Analyser.cpp
## LinkDef.h
