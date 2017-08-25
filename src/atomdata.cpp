#include "config.h"
#include "atomdata.hpp"
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <sstream>
#include <cassert>
#include <algorithm>
#include "slater.hpp"
#include "base_directory.hpp"
#include "debug.hpp"

using namespace std;
using namespace chemistry;

Atomdata::Atomdata()
  : m_dirty(false)
{
}

// FIXME: don't assume the data have 100 lines
// FIXME: don't assume data is sorted
// FIXME: need check all data is valid (i.e. all data should be positive)
Atomdata::Atomdata(int atomicnumber)
  : m_atomicnumber(atomicnumber),
    m_dirty(false)
{
  ostringstream resource;
  resource << "hirshfeld/" << atomicnumber << ".data";

  ifstream in;
  char fname[32];

  snprintf(fname, 32, "%d.data", atomicnumber);
  in.open(fname);
  if(not in) {
    string fname = load_first_data(resource.str());
    if(fname.empty()) {
      debug << string("open ") + resource.str() + " failed.\n";
      m_dirty = true;
      return;
    }
    in.open(fname.c_str());
  }
  
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
    debug << tmp1 << '\t' << tmp2 << '\n';
  }
  m_rm = Slater_radius(atomicnumber);
}

double Atomdata::rm() const {
	return m_rm;
}

// FIXME: use bisect
double Atomdata::density(double r) const{
  assert(*this);
  int idx = 0;
  int size = m_r.size();
  while(idx < size && m_r[idx] < r)
    idx++;
  double result;
  if(idx == 0) {
    result = m_density[0];
  } else if(idx == size) {
    result = m_density[idx - 1];
  } else {
    double x1 = m_r[idx-1];
    double x2 = m_r[idx];
    double y1 = m_density[idx-1];
    double y2 = m_density[idx];
    result = y1 + (y2 - y1) * (r - x1) / (x2 - x1);
  }
  debug << "Atomdata::density: " << number() << '\t' << r << '\t' << idx << '\t' << result << '\n';
  return result;
}

int Atomdata::number() const {
	return m_atomicnumber;
}
