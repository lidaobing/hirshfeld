
// All units is under a.u.
#ifndef ATOM_H
#define ATOM_H
#include <iosfwd>

namespace chemistry {
	class Atom {
	public:
		int atomicnumber;  //Atomic number
		double x;
		double y;
		double z;
		friend std::ostream& operator<<(std::ostream& os, const Atom& atom);
	};
	std::ostream& operator<<(std::ostream& os, const Atom& atom);
};

#endif
