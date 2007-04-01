#include "config.h"

#include "debug.hpp"
#include <iostream>

using namespace std;

void debug(const string& str) {
  clog << str << '\n';
}
