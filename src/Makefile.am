AM_CXXFLAGS = -Wall -Wextra

INCLUDES = -I$(top_srcdir)/lib
LDADD = $(top_builddir)/lib/libgnu.a

bin_PROGRAMS = hirshfeld

commmonsources = atom.cpp atom.hpp atomdata.cpp atomdata.hpp contraction.cpp \
		 contraction.hpp hirshfeld.cpp hirshfeld.hpp lebedev.cpp \
		 lebedev.hpp molecule.cpp molecule.hpp molorbit.hpp slater.cpp \
		 slater.hpp base_directory.cpp
		 
hirshfeld_SOURCES = main.cpp $(commmonsources)

convertdir = $(datadir)/hirshfeld
convert_PROGRAMS = convert
convert_SOURCES = convert.cpp $(commmonsources)
