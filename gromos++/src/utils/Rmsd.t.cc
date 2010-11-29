// utils_Rmsd.t.cc
#include <cassert>
#include <cstdlib>
#include "Rmsd.h"
#include "../fit/RotationalFit.h"
#include "../fit/Reference.h"
#include "../gio/InG96.h"
#include "../gcore/System.h"
#include "../gio/InTopology.h"
#include <string>
#include <iostream>

using namespace gcore;
using namespace gio;
using namespace utils;
using namespace fit;

using namespace std;

int main(int argc, char *argv[]) {
  if (argc != 4) {
    cerr << "Usage: " + string(argv[0]) + " <Topology> <RefG96inates> <G96inates>\n";
    exit(1);
  }
  string top = argv[1];
  string refFile = argv[2];
  string file = argv[3];
  // read Simulation data

  InTopology it(top);

  System refsys(it.system());

  InG96 ic(refFile);
  ic >> refsys;
  ic.close();

  System sys(refsys);

  Reference ref(&refsys);

  ref.addClass(0, "CA");
  ref.addClass(0, "N");
  ref.addClass(0, "C");

  RotationalFit rf(&ref);
  Rmsd rmsd(&ref);

  ic.open(file);

  while (!ic.eof()) {
    ic >> sys;
    rf.fit(&sys);
    double d = rmsd.rmsd(sys);
    cout << d << endl;
  }
  ic.close();
  return 0;
}
