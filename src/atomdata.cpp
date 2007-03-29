#include "atomdata.hpp"
#include "slater.hpp"
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <sstream>

using namespace std;
using namespace chemistry;

Atomdata::Atomdata() {
}

Atomdata::Atomdata(int atomicnumber)
	:m_atomicnumber(atomicnumber)
{
  ostringstream fname;
  fname << "/usr/share/hirshfeld/" << atomicnumber << ".data";
  ifstream in(fname.str().c_str());
  if(!in) {
    cerr << "Can't open file " << fname << '\n';
    exit(1);
  }
  for(int i = 0; i < 100; i++) {
    double tmp1;
    double tmp2;
    in >> tmp1 >> tmp2;
    m_r.push_back(tmp1);
    m_density.push_back(tmp2);
  }
  m_rm = Slater_radius(atomicnumber);
}

double Atomdata::rm() const {
	return m_rm;
}

double Atomdata::density(double r) const{
	int idx = 0;
	int size = m_r.size();
	while(idx < size && m_r[idx] < r)
		idx++;
	if(idx == 0)
		return m_density[0];
	if(idx == size)
		return m_density[idx - 1];
	double x1 = m_r[idx-1];
	double x2 = m_r[idx];
	double y1 = m_density[idx-1];
	double y2 = m_density[idx];

	return y1 + (y2 - y1) * (r - x1) / (x2 - x1);
}

int Atomdata::number() const {
	return m_atomicnumber;
}
