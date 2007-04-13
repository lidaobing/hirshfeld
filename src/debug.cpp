#include "config.h"

#include "debug.hpp"
#include <iostream>

using namespace std;

#ifdef DEBUG

ostream& debug = clog;

#else

DummyOstream theDummyOstream;
DummyOstream& debug = theDummyOstream;

#endif

