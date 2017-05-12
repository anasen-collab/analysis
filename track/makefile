Analyzer: tr_dict.cxx LookUp.cpp Analyzer_Cal.cpp
	@echo compiling Analyzer code...
	g++ -o Analyzer tr_dict.cxx LookUp.cpp Analyzer_Cal.cpp `root-config --cflags --glibs`

tr_dict.cxx: tree_structure.h LinkDef.h
	@echo generating tracking dictionary...
	rootcint -f tr_dict.cxx -c tree_structure.h LinkDef.h

clean:
	@echo removing Analyzer files...
	rm Analyzer tr_dict.cxx tr_dict.h
