// -*- coding: utf-8 -*-
// $Id: hirshfeld.cpp 18 2006-07-11 12:35:43Z nichloas $
#include "config.h"
#include "hirshfeld.hpp"
#include <cstdio>
#include <cmath>
#include <iomanip>
#include <cassert>
#include <fstream>
#include "lebedev.hpp"
#include "atomdata.hpp"
#include "atom.hpp"

using namespace chemistry;
using namespace std;

Hirshfeld::Hirshfeld(istream& is)
  : m_dirty(false)
{
  mol.read(is);

  if(not mol) {
    m_error_str = "format error in fchk file";
    m_dirty = true;
    return;
  }

  Lebedev_Laikov_sphere(110, lebedev_x, lebedev_y, lebedev_z, lebedev_w);

  if(not initatoms()) {
    m_error_str = "init atom data failed";
    m_dirty = true;
    return;
  }
}

Hirshfeld::~Hirshfeld() {
  for(unsigned i = 0; i < atoms.size(); i++) {
    delete atoms[i];
  }
}

void Hirshfeld::run(ostream& os) {
  if(m_dirty) return;
	int atmnum = mol.atmnum();
	os << "No.\tAtomic\telctron\t\tcharge\n";
	for(int i = 0; i < atmnum; i++) {
		activeatom = i;
		double electron = gauss_chebyshev_integrate(30);
		int atomicnumber = mol.atom(i).atomicnumber;
		os << i+1 << "\t"
                   << atomicnumber << "\t"
                   << fixed << setw(10) << electron << "\t"
                   << showpos << atomicnumber - electron << '\n';
		os.unsetf(ios::showpos);
	}
}

// FIXME: for each atom-type, we only need one Atomdata
bool Hirshfeld::initatoms() {
  int atomnum = mol.atmnum();
  for(int i = 0; i < atomnum; i++) {
    int number = mol.atom(i).atomicnumber;
    Atomdata * tmp = new Atomdata(number);
    if(not *tmp) {
      delete tmp;
      return false;
    }
    atoms.push_back(tmp);
  }
  return true;
}

// void Hirshfeld::test() const {
// 	double r = 2.0;
// 	for(int i = 0; i < 110; i++) {
// 		double x = r * lebedev_x[i];
// 		double y = r * lebedev_y[i];
// 		double z = r * lebedev_z[i];
// 		cout << mol->density(x,y,z) << endl;
// 	}
// }


double Hirshfeld::density(double r) const {
	double result = 0.0;
	double x0 = mol.atom(activeatom).x;
	double y0 = mol.atom(activeatom).y;
	double z0 = mol.atom(activeatom).z;

	int atmnum = atoms.size();
	double * atmdensity = new double[atmnum];

	for(int i = 0; i < 110; i++) {
		double x = x0 + r * lebedev_x[i];
		double y = y0 + r * lebedev_y[i];
		double z = z0 + r * lebedev_z[i];
		double w = lebedev_w[i];
		double moldensity = mol.density(x,y,z);
		if(moldensity < 1e-20)
			return 0.0;
		double sumatmdensity = 0.0;
		for(int j = 0; j < atmnum; j++) {
			const Atom & tmpatom = mol.atom(j);
			double dx = x - tmpatom.x;
			double dy = y - tmpatom.y;
			double dz = z - tmpatom.z;
			double r = sqrt(dx * dx
			              + dy * dy
				      + dz * dz);
			atmdensity[j] = atomdensity(j, r);
			sumatmdensity += atmdensity[j];
		}
//		result += w * atmdensity[activeatom];


		result += w * atmdensity[activeatom] / sumatmdensity * moldensity;

	}
	delete []atmdensity;
	return result;
}

double Hirshfeld::sphereint(double x) const{
	double rm = radiusmedium(activeatom);
	double r = rm * (1+x) / (1-x);
	double dens = density(r);
	double result = dens * 4 * M_PI * r * r * 2 * rm /(x-1)/(x-1);
//	cout << rm << "\t" << r << "\t" << dens << endl;
	return result;
}


// double Hirshfeld::atomdensity(int n, double r) const{
// 	if(1 == mol->m_atoms[n].atmtype) {
// 		return atoms[n]->density(r,0,0);
// 	}
// 	double result = 0.0;
// 	for(int i = 0; i < 110; i++) {
// 		double x = r * lebedev_x[i];
// 		double y = r * lebedev_y[i];
// 		double z = r * lebedev_z[i];
// 		double w = lebedev_w[i];
// 		result += w * atoms[n]->density(x,y,z);
// 	}
// 	return result;
// }



//《数值分析》 陈昌明，厦门大学出版社， P109
double Hirshfeld::gauss_chebyshev_integrate(int n) {
	double result = 0.0;
	double coef = M_PI / (n + 1.0);
	for(int i=0; i <= n; i++) {    // confirm "<="
		double x = cos((2*i+1.0)/(2*n+2.0)*M_PI);
		result += coef * sphereint(x) * sqrt(1-x*x);
//		cout << i << "\t" << x << "\t" <<result << endl;
	}
	return result;
}

double Hirshfeld::atomdensity(int atmidx, double r) const {
	return atoms[atmidx]->density(r);
}

double Hirshfeld::radiusmedium(int atmidx) const {
	return atoms[atmidx]->rm();
}
