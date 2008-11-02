#include "config.h"

#include <iostream>
#include <fstream>
#include <cstdlib>
#include "dirname.h"
#include "long-options.h"
#include "hirshfeld.hpp"
using namespace std;

void usage(int) {
  cout << "calculate hirshfeld charge from gaussian's fchk file.\n"
          "Usage: hirshfeld <fchkfile|-> [outputfile]\n"
          "       hirshfeld --help\n"
          "       hirshfeld --version\n"
          "\n";
}

int main(int argc, char* argv[]) {
  parse_long_options(argc,
                     argv,
                     "hirshfeld",
                     PACKAGE_NAME,
                     VERSION,
                     usage,
                     "LI Daobing <lidaobing@gmail.com>");

  if(argc != 2 && argc != 3) {
    usage(0);
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
    cout << hir.error_str() << "\n";
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
