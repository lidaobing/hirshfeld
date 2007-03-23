AM_CXXFLAGS = -Wall -Wextra
INCLUDES = -DPACKAGE_DATA_DIR=\""$(pkgdatadir)"\"
bin_PROGRAMS = hirshfeld

commmonsources = atom.cpp atom.hpp atomdata.cpp atomdata.hpp contraction.cpp \
		 contraction.hpp hirshfeld.cpp hirshfeld.hpp lebedev.cpp \
		 lebedev.hpp molecule.cpp molecule.hpp molorbit.hpp slater.cpp \
		 slater.hpp port.h
		 
hirshfeld_SOURCES = main.cpp $(commmonsources)
hirshfeld_LDADD = @LIBOBJS@

convertdir = $(datadir)/hirshfeld
convert_PROGRAMS = convert
convert_SOURCES = convert.cpp $(commmonsources)
convert_LDADD = @LIBOBJS@
