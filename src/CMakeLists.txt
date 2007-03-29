ADD_LIBRARY(mylib STATIC
	atom.cpp
	atomdata.cpp
	atomdata.hpp
	atom.hpp
	contraction.cpp
	contraction.hpp
	hirshfeld.cpp
	hirshfeld.hpp
	lebedev.cpp
	lebedev.hpp
	molecule.cpp
	molecule.hpp
	molorbit.hpp
	slater.cpp
	slater.hpp
)

ADD_EXECUTABLE(hirshfeld main.cpp) 
TARGET_LINK_LIBRARIES(hirshfeld mylib)

ADD_EXECUTABLE(convert convert.cpp) 
TARGET_LINK_LIBRARIES(convert mylib)

INSTALL(PROGRAMS hirshfeld 
	DESTINATION bin)


INSTALL(PROGRAMS convert
	DESTINATION lib/hirshfeld)
