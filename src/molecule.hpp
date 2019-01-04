#ifndef MOLECULE_H
#define MOLECULE_H

#include <vector>
#include <cstdio>
#include <iosfwd>

namespace chemistry {
  class Atom;
  class Contraction; 
  class Molecule {
    // noncopyable
    Molecule(const Molecule&);
    Molecule& operator=(const Molecule&);
  public:
    Molecule();
    explicit Molecule(std::istream& is);
    ~Molecule();

    operator void*() const {return m_dirty?0:(void *)(1);}

    void read(std::istream& is);
    
    double density(double x, double y, double z) const;
    friend std::ostream&
    ::operator<<(std::ostream& os, const Molecule& mol);
    int atmnum() const;
    const chemistry::Atom& atom(int idx) const;
  private:
    bool m_dirty;
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
