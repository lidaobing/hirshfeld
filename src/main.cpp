#include <iostream>
#include <fstream>
#include "hirshfeld.hpp"
#include "port.h"
using namespace std;

int main(int argc, char* argv[]) {
  if(argc != 2 && argc != 3) {
    cout << "Hirshfeld version " VERSION ", by LI Daobing <lidaobing@gmail.com>\n";
    cout << "Usage: " << basename(argv[0]) << " fchkfile [outputfile]\n";
    exit(1);
  }

  if(argc == 3) {
    ofstream out(argv[2]);
    if(!out) {
      cout << "Can't open output file: " << argv[2] << endl;
      exit(1);
    }
    streambuf * oldcout = cout.rdbuf(out.rdbuf());

    Hirshfeld hir(argv[1]);
    hir.run(cout);
    //	hir.test();

    cout.rdbuf(oldcout);
    out.close();
  } else {
    Hirshfeld hir(argv[1]);
    hir.run(cout);
  }


}
