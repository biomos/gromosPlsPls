#include <cassert>
#include <cstdlib>
#include <iostream>
#include <map>
#include "InTopology.h"
#include "OutTopology.h"
#include "../gcore/GromosForceField.h"
#include "../gcore/System.h"
#include "../gcore/LJException.h"
#include "../gcore/Molecule.h"
#include "../gcore/MoleculeTopology.h"

using namespace gcore;
using namespace gio;

using namespace std;

int main(int argc, char *argv[]) {
  if (argc != 2) {
    std::cerr << "Usage: " + std::string(argv[0]) + " <Topology>\n";
    exit(1);
  }
  try {
    cout << "create an it" << endl;

    InTopology it(argv[1]);

    cout << "done with it" << endl;
    cout << "create a system" << endl;

    System sys(it.system());
    cout << "done with sys" << endl;
    cout << sys.mol(0).topology().numRes() << endl;

    GromosForceField gff(it.forceField());
    std::cout << sys.numMolecules() << endl;
    OutTopology ot(std::cout);
    ot.setTitle(it.title());
    ot.write(sys, gff);

    return 0;
  }  catch (gromos::Exception e) {
    cerr << e.what() << endl;
    return 1;
  }
}
