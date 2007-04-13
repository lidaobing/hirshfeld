#include "config.h"
#include "slater.hpp"

static const double default_slater_radius = 1.35; // unit: Angstrom

//J.C. Slater J. Chem. Phys. 41. 3199(1964)
double Slater_radius(int atmtype) {
	const static double r[] = {
		-1.0,
		0.25,   //H
		-1.0,
		1.45,   //Li
		1.05,
		0.85,
		0.70,
		0.65,
		0.60,
		0.50,
		-1.0,
		1.80,	//Na
		1.50,
		1.25,
		1.10,
		1.00,
		1.00,
		1.00,
		-1.0,
		2.20,
		1.80,
		1.60,
		1.40,
		1.35,
		1.40,
		1.40,
		1.40,
		1.35,
		1.35,
		1.35,
		1.35,
		1.30,
		1.25,
		1.15,
		1.15,
		1.15};   //Br
	int maxnum = sizeof(r)/sizeof(r[0]);
	const double Angstrom2au = 1.88972;
	if(atmtype > 0 
           and atmtype < maxnum
           and r[atmtype] > 0) {
		return r[atmtype] * Angstrom2au;
	} else {
          return default_slater_radius * Angstrom2au;
	}
}

#ifdef TEST
int main() {
	int i;
	for(i = 1; i < 50; i++) {
		printf("%i\t%f\n", i, Slater_radius(i));
	}
}
#endif
