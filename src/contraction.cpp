#include "config.h"
#include "contraction.hpp"
#include <vector>
#include <cassert>
#include <cmath>
#include <iostream>
#include <cstdlib>
#include "atom.hpp"

using chemistry::Contraction;
using namespace std;

static double prim_calcN(double alpha, int majortype, int minortype);
static double prim_calcN_D5(double alpha, int minortype);
static double prim_calcN_gaussian(double alpha, int l, int m, int n);
static int dfactor(int n);

// dfactor(n) = (2n - 1)!!

Contraction::Contraction(Atom & atom, int majortype, int minortype)
	: m_atom(atom),
	  m_majortype(majortype),
	  m_minortype(minortype),
	  m_x(atom.x),
	  m_y(atom.y),
	  m_z(atom.z)
{
//	cout << m_x << m_y << m_z << endl;
}

double Contraction::calc(double x, double y, double z) const {

	double result = 0.0;
	x -= m_x;
	y -= m_y;
	z -= m_z;
	double r2 = x*x + y*y + z*z;
	int primnum = m_prims.size();
	for(int i = 0; i < primnum; i++) {
		result += m_prims[i].coef * m_prims[i].N * exp(-m_prims[i].alpha * r2);
	}


	switch(m_majortype) {
		case S:
			break;
		case P:
			switch(m_minortype) {
				case P_x:
					result *= x;
					break;
				case P_y:
					result *= y;
					break;
				case P_z:
					result *= z;
					break;
				default:
					assert(false);
			}
			break;
		case D:
			switch(m_minortype) {
				case D_xx:
					result *= x*x;
					break;
				case D_yy:
					result *= y*y;
					break;
				case D_zz:
					result *= z*z;
					break;
				case D_xy:
					result *= x*y;
					break;
				case D_xz:
					result *= x*z;
					break;
				case D_yz:
					result *= y*z;
					break;
				default:
					assert(false);
			}
			break;
		case D5:
			switch(m_minortype) {
				case D5_z2:
					result *= (3*z*z - r2);
					break;
				case D5_xz:
					result *= x*z;
					break;
				case D5_yz:
					result *= y*z;
					break;
				case D5_x2y2:
					result *= (x*x - y*y);
					break;
				case D5_xy:
					result *= x*y;
					break;
				default:
					assert(false);
			}
			break;
		default:
			assert(false);

	}
return result;
}


void Contraction::add(double alpha, double coef) {
	Prim prim;
	prim.alpha = alpha;
	prim.N = prim_calcN(alpha, m_majortype, m_minortype);
	prim.coef = coef;
	m_prims.push_back(prim);
}

void Contraction::calcN() {
//	m_N = 1.0;
	assert(m_prims.size() != 0);
	if(1 == m_prims.size()) {
		m_N = 1.0;
		return;
	}
	if(m_majortype == D5) {
		printf("This program can't deal with 5D shell contraction!\n");
		exit(1);
	}


	double integrator = 0;

	for(int i = 0; i < int(m_prims.size()); i++) {
		integrator += m_prims[i].coef * m_prims[i].coef;
		for(int j = i+1; j < int(m_prims.size()); j++) {
			integrator += int2prim(i,j) * 2;
		}
	}

	m_N = 1 / sqrt(integrator);
}

double Contraction::int2prim(int index1, int index2) const {
	int l;
	int m;
	int n;
	switch(m_majortype) {
		case S:
			l = 0;
			m = 0;
			n = 0;
			break;
		case P:
			l = 1;
			m = 0;
			n = 0;
			break;
		case D:
			if(m_minortype == D_xx || m_minortype == D_yy || m_minortype == D_zz) {
				l = 2;
				m = 0;
				n = 0;
			} else {
				l = 1;
				m = 1;
				n = 0;
			}
			break;
		default:
			assert(false);
	}

	const Prim& prim1 = m_prims[index1];
	const Prim& prim2 = m_prims[index2];

	double N2 = pow(((prim1.alpha + prim2.alpha) / M_PI), 1.5) *
		pow(2 * (prim1.alpha + prim2.alpha), l+m+n) /
		(dfactor(l) * dfactor(m) * dfactor(n));
	return prim1.coef * prim2.coef * prim1.N * prim2.N / N2;
}

ostream& chemistry::operator<<(ostream& os, const Contraction& con) {
	os << con.m_atom.atomicnumber << "\t"
	   << con.m_majortype << "\t"
	   << con.m_minortype << "\t"
	   << con.m_prims.size() << "\t"
	   << fixed << con.m_N
	   << '\n';
	return os;
}

static int dfactor(int n) {
	assert(n >= 0);
	int result = 1;
	for(int i = 3; i <= 2*n-1; i += 2) {
		result *= i;
	}
	return result;
}

inline double prim_calcN_gaussian(double alpha, int l, int m, int n) {
	return pow((2 * alpha / M_PI), 0.75) *
		sqrt( pow(4 * alpha, l+m+n) /
			(dfactor(l) * dfactor(m) * dfactor(n))
			);
}

double prim_calcN(double alpha, int majortype, int minortype) {
	if(majortype == Contraction::S)
		return prim_calcN_gaussian(alpha, 0,0,0);
	if(majortype == Contraction::P)
		return prim_calcN_gaussian(alpha, 1,0,0);
	if(majortype == Contraction::D)
	{
		if(minortype == Contraction::D_xx ||
		   minortype == Contraction::D_yy ||
		   minortype == Contraction::D_zz)
			return prim_calcN_gaussian(alpha, 2,0,0);
		else
			return prim_calcN_gaussian(alpha, 1,1,0);
	}
	if(majortype == Contraction::D5)
		return prim_calcN_D5(alpha, minortype);
	assert(false && "F shell underconstruction!!");

}

double prim_calcN_D5(double alpha, int minortype) {
	if(minortype == Contraction::D5_xz ||
	   minortype == Contraction::D5_yz ||
       minortype == Contraction::D5_xy)
		return prim_calcN_gaussian(alpha, 1,1,0);

	//d_{3z^2-r^2}
	if(minortype == Contraction::D5_z2) {
		return prim_calcN_gaussian(alpha, 2,0,0) * 0.5;
	} else {
	//d_{x^2-y^2}
		return prim_calcN_gaussian(alpha, 2,0,0) * sqrt(3.0) * 0.5;
	}
}
