#include "molecule.hpp"
#include <cassert>
#include <cstring>
#include <iostream>
#include "contraction.hpp"
#include "atom.hpp"

using namespace std;
using namespace chemistry;

static std::vector<int> atoms;
static std::vector<double> coords;

//static int alpha_electrons;
//static int beta_electrons;

static std::vector<int> shell_type;
static std::vector<int> npps; //Number of primitives per shell
static std::vector<int> shell_map;

static int pri_idx = 0;
static std::vector<double> pri_exp; //Primitive exponents
static std::vector<double> con_coe; //Contraction coefficients
static std::vector<double> p_con_coe; //P(S=P) Contraction coefficients
//static std::vector<double> alpha_mo_coe; //Alpha MO coefficients
//static std::vector<double> beta_mo_coe;  //Beta MO coefficients
static std::vector<double> total_SCF_density;  //Total SCF Density

static void readintvector(FILE* fp, vector<int>& vi, const char * keyword);
static void readdoublevector(FILE* fp, vector<double>& vi, const char * keyword);
//static void readint(FILE* fp, int & i, const char * keyword);
static bool hasspshell();


Molecule::Molecule() {
}


Molecule::Molecule(char * fname) {
  read(fname);
}

Molecule::~Molecule() {
  for(int i = 0; i < int(conts.size()); i++) {
    delete conts[i];
  }
}

int Molecule::atmnum() const {
  return m_atoms.size();
}

const Atom& Molecule::atom(int idx) const {
  return m_atoms[idx];
}

void Molecule::read(char * fname) {
  FILE * fp;
  if((fp = fopen(fname, "r")) != NULL) {
    read(fp);
    fclose(fp);
  } else {
    cerr << "Can't open file " << fname << endl;
    exit(1);
  }
}

void Molecule::read(FILE * fp) {
  char tmpstr[100];

  //read title(A72)
  fgets(tmpstr, 80, fp);
  //read type, method, basis
  fgets(tmpstr, 80, fp);
  if('R' == tmpstr[10]) {
    closeshell = true;
  } else if('U' == tmpstr[10]) {
    closeshell = false;
  } else {
    cerr << "Can't judge the Molecule is close shell or open shell\n" << endl;
    exit(1);
  }

  //	readint(fp, alpha_electrons, "Number of alpha electrons");
  //	readint(fp, beta_electrons, "Number of beta electrons");

  readintvector(fp, atoms, "Atomic numbers");
  readdoublevector(fp, coords,"Current cartesian coordinates"); // in fchk file, units is a.u.
  readintvector(fp, shell_type, "Shell types");
  readintvector(fp, npps, "Number of primitives per shell");
  readintvector(fp, shell_map, "Shell to atom map");
  readdoublevector(fp, pri_exp, "Primitive exponents");
  readdoublevector(fp, con_coe, "Contraction coefficients");
  if(hasspshell()) readdoublevector(fp, p_con_coe, "P(S=P) Contraction coefficients");
  //	readdoublevector(fp, alpha_mo_coe, "Alpha MO coefficients");
  //	if(!closeshell) {
  //		readdoublevector(fp, beta_mo_coe, "Beta MO coefficients");
  //	}
  readdoublevector(fp, total_SCF_density, "Total SCF Density");

  check();
  init();
}

void Molecule::check() {
}

void Molecule::init() {
  initatoms();
  initconts();
  initdensitymatrix();
}

void Molecule::initatoms() {
  int atmnum = atoms.size();
  for(int i = 0; i < atmnum; i++) {
    Atom tmpatm;
    tmpatm.atomicnumber = atoms[i];
    tmpatm.x = coords[3*i];
    tmpatm.y = coords[3*i+1];
    tmpatm.z = coords[3*i+2];
    m_atoms.push_back(tmpatm);
  }
}

void Molecule::initconts() {
  for(int i = 0; i < int(shell_type.size()); i++) {
    switch(shell_type[i]) {
    case Contraction::S:
      init_S(i);
      break;
    case Contraction::P:
      init_P(i);
      break;
    case Contraction::SP:
      init_SP(i);
      break;
    case Contraction::D:
      init_D(i);
      break;
    case Contraction::D5:
      init_D5(i);
      break;
    default:
      assert(false);
    }
  }
}


void Molecule::initdensitymatrix() {
  int contnum = conts.size();
  int density_idx = 0;

  for(int i = 0; i < contnum; i++) {
    density_matrix.push_back(vector<double>());
    for(int j = 0; j <= i; j++) {
      density_matrix[i].push_back(total_SCF_density[density_idx++]);
    }
  }
  for(int i = 0; i < contnum; i++) {
    for(int j = i+1; j < contnum; j++) {
      density_matrix[i].push_back(density_matrix[j][i]);
    }
  }
}

void Molecule::init_S(int shell_idx) {
  int atmnum = shell_map[shell_idx] - 1;

  Contraction * cont = new Contraction(m_atoms[atmnum], Contraction::S, Contraction::S_0);

  for(int j = 0; j < npps[shell_idx]; j++) {
    cont->add(pri_exp[pri_idx], con_coe[pri_idx]);
    pri_idx++;
  }
  cont->calcN();
  conts.push_back(cont);
}

void Molecule::init_P(int shell_idx) {
  int atmnum = shell_map[shell_idx] - 1;

  Contraction * cont[3];
  cont[0] = new Contraction(m_atoms[atmnum], Contraction::P, Contraction::P_x);
  cont[1] = new Contraction(m_atoms[atmnum], Contraction::P, Contraction::P_y);
  cont[2] = new Contraction(m_atoms[atmnum], Contraction::P, Contraction::P_z);

  for(int j = 0; j < npps[shell_idx]; j++) {
    cont[0]->add(pri_exp[pri_idx], con_coe[pri_idx]);
    cont[1]->add(pri_exp[pri_idx], con_coe[pri_idx]);
    cont[2]->add(pri_exp[pri_idx], con_coe[pri_idx]);
    pri_idx++;
  }
  for(int i = 0; i < 3; i++) {
    cont[i]->calcN();
    conts.push_back(cont[i]);
  }
}

void Molecule::init_SP(int shell_idx) {
  int atmnum = shell_map[shell_idx] - 1;

  Contraction * cont[4];
  cont[0] = new Contraction(m_atoms[atmnum], Contraction::S, Contraction::S_0);
  cont[1] = new Contraction(m_atoms[atmnum], Contraction::P, Contraction::P_x);
  cont[2] = new Contraction(m_atoms[atmnum], Contraction::P, Contraction::P_y);
  cont[3] = new Contraction(m_atoms[atmnum], Contraction::P, Contraction::P_z);

  for(int j = 0; j < npps[shell_idx]; j++) {
    cont[0]->add(pri_exp[pri_idx], con_coe[pri_idx]);
    cont[1]->add(pri_exp[pri_idx], p_con_coe[pri_idx]);
    cont[2]->add(pri_exp[pri_idx], p_con_coe[pri_idx]);
    cont[3]->add(pri_exp[pri_idx], p_con_coe[pri_idx]);
    pri_idx++;
  }
  for(int i = 0; i < 4; i++) {
    cont[i]->calcN();
    conts.push_back(cont[i]);
  }
}

void Molecule::init_D(int shell_idx) {
  int atmnum = shell_map[shell_idx] - 1;

  Contraction * cont[6];
  cont[0] = new Contraction(m_atoms[atmnum], Contraction::D, Contraction::D_xx);
  cont[1] = new Contraction(m_atoms[atmnum], Contraction::D, Contraction::D_yy);
  cont[2] = new Contraction(m_atoms[atmnum], Contraction::D, Contraction::D_zz);
  cont[3] = new Contraction(m_atoms[atmnum], Contraction::D, Contraction::D_xy);
  cont[4] = new Contraction(m_atoms[atmnum], Contraction::D, Contraction::D_xz);
  cont[5] = new Contraction(m_atoms[atmnum], Contraction::D, Contraction::D_yz);

  for(int j = 0; j < npps[shell_idx]; j++) {
    for(int i = 0; i < 6; i++) {
      cont[i]->add(pri_exp[pri_idx], con_coe[pri_idx]);
    }
    pri_idx++;
  }
  for(int i = 0; i < 6; i++) {
    cont[i]->calcN();
    conts.push_back(cont[i]);
  }
}

void Molecule::init_D5(int shell_idx)
{
  int atmnum = shell_map[shell_idx] - 1;

  Contraction * cont[5];
  cont[0] = new Contraction(m_atoms[atmnum], Contraction::D5, Contraction::D5_z2);
  cont[1] = new Contraction(m_atoms[atmnum], Contraction::D5, Contraction::D5_xz);
  cont[2] = new Contraction(m_atoms[atmnum], Contraction::D5, Contraction::D5_yz);
  cont[3] = new Contraction(m_atoms[atmnum], Contraction::D5, Contraction::D5_x2y2);
  cont[4] = new Contraction(m_atoms[atmnum], Contraction::D5, Contraction::D5_xy);

  for(int j = 0; j < npps[shell_idx]; j++)
  {
    for(int i = 0; i < 5; i++)
    {
      cont[i]->add(pri_exp[pri_idx], con_coe[pri_idx]);
    }
    pri_idx++;
  }

  for(int i = 0; i < 5; i++)
  {
    cont[i]->calcN();
    conts.push_back(cont[i]);
  }
}

double Molecule::density(double x, double y, double z) const {
  int contnum = conts.size();
  double * contcalc = new double[contnum];
  for(int i = 0; i < contnum; i++) {
    contcalc[i] = conts[i]->calc(x,y,z);
  }
  double result = 0.0;
  for(int i = 0; i <contnum; i++) {
    for(int j = 0; j < i; j++) {
      result += 2 * density_matrix[i][j] * contcalc[i] * contcalc[j];
    }
    result += density_matrix[i][i] * contcalc[i] * contcalc[i];
  }
  delete []contcalc;
  return result;
}

ostream& chemistry::operator<<(ostream& os, const Molecule& mol) {
  os << "Atomic Number          x               y               z" << endl;
  int atmnum = mol.m_atoms.size();
  for(int i = 0; i < atmnum; i++) {
    os << mol.m_atoms[i];
  }
  os << "atom\tmajort\tminort\tprims\tm_N" << endl;
  int contnum = mol.conts.size();
  for(int i = 0; i < contnum; i++) {
    os << *(mol.conts[i]);
  }
  return os;
}


static void readintvector(FILE* fp, vector<int>& vi, const char * keyword) {
  vi.clear();
  char tmpstr[100];
  do {
    fgets(tmpstr, 80, fp);
  } while((strstr(tmpstr, keyword) == NULL) && !feof(fp));
  if(feof(fp))
    printf(keyword);
  int num = atoi(&tmpstr[49]);
  for(int i = 0; i < num; i++) {
    int tmpint;
    fscanf(fp, "%12i", &tmpint);
    vi.push_back(tmpint);
  }
}

static void readdoublevector(FILE* fp, vector<double>& vi, const char * keyword) {
  vi.clear();
  char tmpstr[100];
  do {
    fgets(tmpstr, 80, fp);
  } while(strstr(tmpstr, keyword) == NULL);
  int num = atoi(&tmpstr[49]);
  for(int i = 0; i < num; i++) {
    double tmpdouble;
    fscanf(fp, "%le", &tmpdouble);
    vi.push_back(tmpdouble);
  }
}

// static void readint(FILE* fp, int & i, const char * keyword) {
// 	char tmpstr[100];
// 	do {
// 		fgets(tmpstr, 80, fp);
// 	} while(strstr(tmpstr, keyword) == NULL);
// 	i = atoi(&tmpstr[49]);
// }

static bool hasspshell() {
  const int SP = -1;
  int shellnum = shell_type.size();
  for(int i = 0; i < shellnum; i++) {
    if(SP == shell_type[i])
      return true;
  }
  return false;
}

#ifdef TEST
void printintvector(const vector<int>& vi, char* str) {
  printf("%s:%i\n", str, vi.size());
  for(unsigned int i = 0; i < vi.size(); i++) {
    printf("%12i", vi[i]);
    if(i%6 == 5) printf("\n");
  }
  printf("\n");
}

void printdoublevector(const vector<double>& vd, char* str) {
  printf("%s:%i\n", str, vd.size());
  for(unsigned int i = 0; i < vd.size(); i++) {
    printf("%16.8e", vd[i]);
    if(i%5 == 4) printf("\n");
  }
  printf("\n");
}



void print() {
  printintvector(atoms, "atoms");
  printdoublevector(coords, "coords");
  printintvector(shell_type, "shell_type");
  printintvector(npps, "npps");
  printintvector(shell_map, "shell_map");
  printdoublevector(pri_exp, "pri_exp");
  printdoublevector(con_coe, "con_coe");
  printdoublevector(p_con_coe, "p_con_coe");
  printdoublevector(alpha_mo_coe, "alpha_mo_coe");

}


int main() {
  FILE * fp;
  fp = fopen("HCN.fchk", "r");
  Molecule mol;
  mol.read(fp);
  fclose(fp);
  print();

}

#endif
