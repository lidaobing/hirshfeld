#include "atomdata.hpp"
#include "slater.hpp"
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>

using namespace std;
using namespace chemistry;

Atomdata::Atomdata() {
}

Atomdata::Atomdata(int atomicnumber)
	:m_atomicnumber(atomicnumber)
{
	char fname[100];
	snprintf(fname, 100, PACKAGE_DATA_DIR "/%i.data", atomicnumber);
	ifstream in(fname, ios::in);
	if(!in) {
		cerr << "Can't open file " << fname << endl;
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

Atomdata::Atomdata(const Atomdata& A)
	: m_r(A.m_r),
	  m_density(A.m_density),
	  m_rm(A.m_rm),
	  m_atomicnumber(A.m_atomicnumber)
{
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
