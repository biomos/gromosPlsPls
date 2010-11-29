#include <iostream>
#include <cassert>
#include <map>
#include <cstdlib>
#include "InTopology.h"
#include "InG96.h"
#include "OutGromacs.h"
#include "../gcore/GromosForceField.h"
#include "../gcore/System.h"

using namespace std;

using namespace gcore;
using namespace gio;

int main(int argc, char *argv[]) {
  if (argc != 2 && argc != 3) {
    cerr << "Usage: " + string(argv[0]) + " <Topology>   [<coordinate file>]\n";
    exit(1);
  }
  try {
    InTopology it(argv[1]);
    System sys(it.system());
    InG96 ic;
    if (argc == 3) {
      ic.open(argv[2]);
      ic.select("ALL");
      ic >> sys;
    }

    GromosForceField gff(it.forceField());
    OutGromacs ot(cout);
    ot.setTitle(it.title());
    ot.write(sys, gff);

    return 0;
  }  catch (gromos::Exception e) {
    cerr << e.what() << endl;
    return 1;
  }
}

