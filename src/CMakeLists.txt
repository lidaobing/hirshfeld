include_directories("${PROJECT_BINARY_DIR}")
include_directories("${PROJECT_SOURCE_DIR}/lib")

add_executable (hirshfeld
	atom.cpp
	atomdata.cpp
	base_directory.cpp
	contraction.cpp
	debug.cpp
	hirshfeld.cpp
	lebedev.cpp
	main.cpp
	molecule.cpp
	slater.cpp
	)

install (TARGETS hirshfeld DESTINATION bin)

add_executable (convert
	atom.cpp
	atomdata.cpp
	base_directory.cpp
	contraction.cpp
	convert.cpp
	debug.cpp
	hirshfeld.cpp
	lebedev.cpp
	molecule.cpp
	slater.cpp
	)

install (TARGETS convert DESTINATION ${CMAKE_INSTALL_LIBDIR}/${PACKAGE_NAME})