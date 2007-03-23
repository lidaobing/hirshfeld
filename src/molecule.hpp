#ifndef MOLECULE_H
#define MOLECULE_H

#include <vector>
#include <cstdio>
#include <iosfwd>

namespace chemistry {
  class Atom;
  class Contraction; 
  class Molecule {
  public:
    Molecule(const char * fname);
    Molecule();
    ~Molecule();
    void read(std::FILE * f);
    void read(const char * filename);
    double density(double x, double y, double z) const;
    friend std::ostream&
    operator<<(std::ostream& os, const Molecule& mol);
    int atmnum() const;
    const chemistry::Atom& atom(int idx) const;
  private:
    bool closeshell;
    std::vector<chemistry::Atom> m_atoms;
    std::vector<chemistry::Contraction *> conts;
    std::vector<std::vector<double> > density_matrix;
    void check();
    void init();
    void initatoms();
    void initconts();
    void initdensitymatrix();
    void init_S(int shell_idx);
    void init_P(int shell_idx);
    void init_SP(int shell_idx);
    void init_D(int shell_idx);
    void init_D5(int shell_idx);
  };
}
#endif
