#include <iostream>
#include <fstream>
#include "hirshfeld.hpp"
#include "port.h"
using namespace std;

int main(int argc, char* argv[]) {
  if(argc != 2 && argc != 3) {
    cout << "Hirshfeld version " VERSION ", by LI Daobing <lidaobing@gmail.com>\n";
    cout << "Usage: " << basename(argv[0]) << " fchkfile|- [outputfile]\n";
    exit(1);
  }

  istream* p_is;
  ostream* p_os;

  ifstream ifs;

  if(argv[1] == string("-")) {
    p_is = &cin;
  } else {
    ifs.open(argv[1]);
    if(not ifs) {
      cerr << "can't open input file: " << argv[1] << "\n";
      exit(1);
    }
    p_is = &ifs;
  }

  Hirshfeld hir(*p_is);
  if(not hir) {
    cout << "format error in fchk file: " << argv[1] << "\n";
    exit(1);
  }

  ofstream ofs;
  if(argc == 2) {
    p_os = &cout;
  } else {
    ofs.open(argv[2]);
    if(!ofs) {
      cerr << "Can't open output file: " << argv[2] << '\n';
      exit(1);
    }
    p_os = &ofs;
  }
  hir.run(*p_os);

  ifs.close();
  ofs.close();
}
