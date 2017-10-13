Main: Main_dict.cxx Main.cpp ChannelMap.h Silicon_Cluster.h
	@echo compiling Main code...
	g++ -o Main Main_dict.cxx Main.cpp `root-config --cflags --glibs` -O3

Main_dict.cxx: ../include/tree_structure.h ../include/LinkDef.h
	@echo generating Main dictionary...
	rootcint -f Main_dict.cxx -c ../include/tree_structure.h ../includeLinkDef.h

clean:
	@echo removing Main files...
	rm Main Main_dict.cxx Main_dict.h
