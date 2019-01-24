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
target_include_directories(hirshfeld PRIVATE ${PROJECT_BINARY_DIR})
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
target_include_directories(convert PRIVATE ${PROJECT_BINARY_DIR})
install (TARGETS convert DESTINATION ${CMAKE_INSTALL_LIBDIR}/${PROJECT_NAME})
