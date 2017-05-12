B: tr_dict.cxx Main_ dict.cxx Main.cpp
	@echo compiling Analyzer code...
	g++ -o B tr_dict.cxx LookUp.cpp Analyzer.cpp `root-config --cflags --glibs`

tr_dict.cxx: tree_structure.h LinkDef.h
	@echo generating tracking dictionary...
	rootcint -f tr_dict.cxx -c tree_structure.h LinkDef.h

clean:
	@echo removing B files...
	rm B tr_dict.cxx
