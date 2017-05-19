all: Analyzer Analyzer_Cal Analyzer_ES Analyzer_Maria

Analyzer: tr_dict.cxx LookUp.cpp Analyzer.cpp
	@echo compiling Analyzer code...
	g++ -o Analyzer tr_dict.cxx LookUp.cpp Analyzer.cpp `root-config --cflags --glibs`

Analyzer_Cal: tr_dict.cxx LookUp.cpp Analyzer_Cal.cpp
	@echo compiling Analyzer code...
	g++ -o Analyzer_Cal tr_dict.cxx LookUp.cpp Analyzer_Cal.cpp `root-config --cflags --glibs`

Analyzer_ES: tr_dict.cxx LookUp.cpp Analyzer_ES.cpp
	@echo compiling Analyzer code...
	g++ -o Analyzer_ES tr_dict.cxx LookUp.cpp Analyzer_ES.cpp `root-config --cflags --glibs`

Analyzer_Maria: tr_dict.cxx LookUp.cpp Analyzer_Maria.cpp
	@echo compiling Analyzer code...
	g++ -o Analyzer_Maria tr_dict.cxx LookUp.cpp Analyzer_Maria.cpp `root-config --cflags --glibs`

tr_dict.cxx: tree_structure.h LinkDef.h
	@echo generating tracking dictionary...
	rootcint -f tr_dict.cxx -c tree_structure.h LinkDef.h

clean:
	@echo removing Analyzer files...
	rm Analyzer Analyzer_Cal Analyzer_ES Analyzer_Maria tr_dict.cxx tr_dict.h
