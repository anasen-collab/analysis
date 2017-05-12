Main: Main_dict.cxx Main.cpp
	@echo compiling Main code...
	g++ -o Main Main_dict.cxx Main.cpp `root-config --cflags --glibs` -O3

Main_dict.cxx: tree_structure.h LinkDef.h
	@echo generating Main dictionary...
	rootcint -f Main_dict.cxx -c tree_structure.h LinkDef.h

clean:
	@echo removing Main files...
	rm Main Main_dict.cxx Main_dict.h
