NETCDF = -lnetcdf_c++4 -lnetcdf
DEP = describe.cpp vecop.cpp
APP = main.cpp
EXECUTABLE = program

install: $(APP) $(DEP)
	g++ -o $(EXECUTABLE) $(APP) $(DEP) $(NETCDF) -std=c++11 -I.
