Main: Main_dict.cxx Main.cpp
	g++ -o Main Main_dict.cxx Main.cpp `root-config --cflags --glibs` -O3

Main_dict.cxx: tree_structure.h LinkDef.h
	rootcint -f Main_dict.cxx -c tree_structure.h LinkDef.h
