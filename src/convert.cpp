// $Id: convert.cpp 19 2007-01-09 16:35:00Z nichloas $

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include "dirname.h"
#include "molecule.hpp"
#include "lebedev.hpp"
#include "slater.hpp"
#include "atom.hpp"

using namespace std;
using namespace chemistry;

static const int LEBEDEVNUM = 110;
static const int RADIUSNUM = 100;
static double lebedev_x[LEBEDEVNUM];
static double lebedev_y[LEBEDEVNUM];
static double lebedev_z[LEBEDEVNUM];
static double lebedev_w[LEBEDEVNUM];

static void convert(const Molecule& mol);
static double density(const Molecule& mol, double r);
static void initlebedev();

int main(int argc, char * argv[]) {
  if(argc != 2) {
    cerr << "Usage: " << basename(argv[0]) << " *.fchk\n";
    exit(1);
  }
  ifstream ifs(argv[1]);
  if(not ifs) {
    cerr << "can't open file: " << argv[1] << "\n";
    exit(1);
  }
  Molecule mol(ifs);
  ifs.close();
  
  if(mol.atmnum() != 1) {
    cerr << "Only for atom's fchk file!\n";
    exit(1);
  }
  convert(mol);
  return 0;
}

void convert(const Molecule & mol) {
  initlebedev();
  int atomicno = mol.atom(0).atomicnumber;
  double rm = Slater_radius(atomicno);
  for(int i = 0; i < RADIUSNUM; i++) {
    double x = cos((i+0.5)/RADIUSNUM*M_PI);
    double r = rm * (1 - x) / (1 + x);
    cout << scientific << r << "\t" << density(mol, r) << '\n';
  }
}

double density(const Molecule& mol, double r) {
  double result = 0.0;
  for(int i = 0; i < LEBEDEVNUM; i++) {
    result += lebedev_w[i] *
    mol.density( r * lebedev_x[i], r * lebedev_y[i], r * lebedev_z[i]);
  }
  return result;
}

inline void initlebedev() {
  Lebedev_Laikov_sphere (LEBEDEVNUM, lebedev_x, lebedev_y, lebedev_z, lebedev_w);
}
