#include "atom.hpp"
#include <iomanip>
using namespace std;

namespace chemistry {

ostream& operator<<(ostream& os, const Atom& atom) {
	os << setw(7) << atom.atomicnumber << "      "
	   << setw(16) << showpoint << setprecision(8) << fixed << atom.x
	   << setw(16) << atom.y
	   << setw(16) << atom.z
	   << '\n';
	return os;
}
}

