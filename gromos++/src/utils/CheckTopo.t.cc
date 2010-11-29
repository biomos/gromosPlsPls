#include <set>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include "CheckTopo.h"
#include "../gio/InTopology.h"
#include "../gcore/System.h"
#include "../gcore/LJException.h"
#include "../gcore/MoleculeTopology.h"
#include "../gcore/Molecule.h"
#include "../gcore/Solvent.h"
#include "../gio/InG96.h"

using namespace gcore;
using namespace gio;
using namespace utils;

using namespace std;

int main(int argc, char *argv[]) {
  if (argc != 2) {
    cerr << "Usage: " + string(argv[0]) + " <Topology>\n";
    exit(1);
  }
  try {
    InTopology it(argv[1]);
    System sys(it.system());

    CheckTopo ct(sys.mol(0).topology());
    ct.checkBonds();
    ct.checkAngles();
    ct.checkImpropers();
    ct.checkChargeGroups();
    ct.checkExclusions();

    cout << ct.numErrors() << endl;
    for (int i = 0; i < ct.numErrors(); i++) {
      cout << ct.error(i) << endl;
    }
    return 0;
  }  catch (gromos::Exception e) {
    cerr << e.what() << endl;
    return 1;
  }
}
