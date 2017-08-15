Organize: organize_dictionary.cxx Organize.cpp
	@echo compiling Organize code...
	g++ -o Organize organize_dictionary.cxx Organize.cpp `root-config --cflags --glibs`	

organize_dictionary.cxx: ../include/organizetree.h LinkDef.h
	@echo generating organize dictionary...
	rootcint -f organize_dictionary.cxx -c ../include/organizetree.h LinkDef.h

clean:
	@echo removing Organize files...
	rm Organzie organize_dictionary.cxx organize_dictionary.h
