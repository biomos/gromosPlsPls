// fit_PositionUtils.t.cc

#include <cassert>
#include <cassert>
#include <cstdlib>
#include "PositionUtils.h"
#include "../gio/InG96.h"
#include "../gcore/System.h"
#include "../gio/InTopology.h"
#include "../gio/OutG96.h"
#include "../gmath/Vec.h"
#include <string>
#include <iostream>
#include <vector>

using namespace gcore;
using namespace gio;
using namespace fit;
using namespace gmath;
using namespace std;

int debug_level = 0;

int main(int argc, char *argv[]) {
  if (argc != 3) {
    cerr << "Usage: " + string(argv[0]) + " <Topology> <G96inates>\n";
    exit(1);
  }
  string top = argv[1];
  string file = argv[2];
  // read Simulation data

  InTopology it(top);
  System sys(it.system());

  InG96 ic;
  ic.open(file);


  OutG96 oc;
  oc.open(cout);
  oc.writeTitle(ic.title());

  while (!ic.eof()) {
    ic >> sys;
    Vec com = PositionUtils::com(sys);
    cout << "Center of Mass: " << com[0] << ' ' << com[1]
            << ' ' << com[2] << endl;
    Vec cog = PositionUtils::cog(sys);
    cout << "Center of Geometry: " << cog[0] << ' ' << cog[1]
            << ' ' << cog[2] << endl;
    PositionUtils::translate(&sys, -com);
    cout << "Mol(0) fitted to center of mass:\n";
    oc << sys;

  }
  return 0;
}
