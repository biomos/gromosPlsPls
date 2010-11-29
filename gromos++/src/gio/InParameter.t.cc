#include <cassert>
#include <map>
#include <iostream>
#include <cstdlib>
#include "InParameter.h"
#include "OutTopology.h"
#include "../gcore/GromosForceField.h"
#include "../gcore/System.h"

using namespace std;
using namespace gcore;
using namespace gio;

int main(int argc, char *argv[]) {
  if (argc != 2) {
    cerr << "Usage: " + string(argv[0]) + " <Parameterfile>\n";
    exit(1);
  }
  try {
    InParameter ip(argv[1]);
    GromosForceField gff(ip.forceField());
    System sys;

    OutTopology ot(cout);
    ot.setTitle(ip.title());
    ot.write(sys, gff);

    return 0;
  }  catch (gromos::Exception e) {
    cerr << e.what() << endl;
    return 1;
  }
}
