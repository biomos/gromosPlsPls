#include <cassert>
#include <set>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include "../gio/InTopology.h"
#include "../gmath/Vec.h"
#include "AtomSpecifier.h"
#include "PropertyContainer.h"
#include "../gcore/System.h"
#include "../gcore/MoleculeTopology.h"
#include "../gcore/Molecule.h"
#include "../bound/Boundary.h"
#include "../bound/RectBox.h"

using namespace gcore;
using namespace gio;
using namespace utils;

using namespace std;

int main(int argc, char *argv[]) {
  if (argc != 3) {
    cerr << "Usage: " + string(argv[0]) + " <Topology> <propertyspecifier>\n";
    exit(1);
  }
  try {
    InTopology it(argv[1]);
    System sys(it.system());
    string s = argv[2];

    bound::Boundary *pbc;
    pbc = new bound::RectBox(&sys);

    PropertyContainer bs(sys, NULL);
    bs.addSpecifier(argv[2]);

    cout << bs << endl;

    return 0;
  }  catch (gromos::Exception e) {
    cerr << e.what() << endl;
    return 1;
  }
}
