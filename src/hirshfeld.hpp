#ifndef HIRSHFELD_H
#define HIRSHFELD_H

#include <vector>
#include <iosfwd>
#include <string>

#include "molecule.hpp"

namespace chemistry {
  class Atomdata;
}

class Hirshfeld{
  // noncopyable
  Hirshfeld(const Hirshfeld&);
  Hirshfeld& operator=(const Hirshfeld&);
public:
  explicit Hirshfeld(std::istream& is);

  operator void*() const {return m_dirty?0:(void *)(1);}
  const std::string& error_str() const {return m_error_str;}

  void run(std::ostream& os);
  ~Hirshfeld();
private:
  bool m_dirty;
  std::string m_error_str;
  
  chemistry::Molecule mol;
  std::vector<chemistry::Atomdata *> atoms;
  int activeatom;

  double gauss_chebyshev_integrate(int n);
//	double atomdensity(int n, double r) const;
  double lebedev_x[110];
  double lebedev_y[110]; 
  double lebedev_z[110];
  double lebedev_w[110];
  double density(double r) const;
  double sphereint(double r) const;
  bool initatoms();
  double atomdensity(int atomicnumber, double r) const;
  double radiusmedium(int atmidx) const;
};

#endif
