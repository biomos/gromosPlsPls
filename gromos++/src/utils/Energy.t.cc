#include <cassert>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include "../gio/InTopology.h"
#include "../gio/InG96.h"
#include "Energy.h"
#include "AtomSpecifier.h"
#include "PropertyContainer.h"
#include "../gcore/System.h"
#include "../gcore/GromosForceField.h"
#include "../bound/Boundary.h"
#include "../bound/RectBox.h"
#include "../gcore/MoleculeTopology.h"
#include "../gcore/Molecule.h"
#include "../gmath/Vec.h"

using namespace gcore;
using namespace gio;
using namespace utils;
using namespace gmath;


using namespace std;

int main(int argc, char *argv[]) {
  if (argc != 3) {
    cerr << "Usage: " + string(argv[0]) + " <Topology> <coordinates>\n";
    exit(1);
  }
  try {
    InTopology it(argv[1]);
    System sys(it.system());
    GromosForceField gff(it.forceField());

    InG96 ic(argv[2]);
    ic.select("ALL");
    ic >> sys;
    ic.close();

    bound::Boundary *pbc;
    pbc = new bound::RectBox(&sys);

    /*
    string s="1:20";
    string t="a%1:1,2,3";
    cout << s << " " << t << endl;
    AtomSpecifier as(sys, s);
    PropertyContainer pc(sys, pbc);
    pc.addSpecifier(t);

    Energy en(sys, gff, *pbc);
    en.setAtoms(as);
    en.setProperties(pc);
    en.setCutOff(1.4);
    en.setRF(62.0, 0.0);
    en.calc();
    cout << en.vdw_m(0) << "\t" << en.el_m(0) << "\t"
         << en.vdw_s(0) << "\t" << en.el_s(0) << endl;
    cout << en.vdw(0) << "\t" << en.el(0) << endl;
    cout << en.vdw() << "\t" << en.el() << endl;
    cout << en.nb() << endl;
    cout << en.cov(0) << endl;
    cout << en.tot() << endl;
     */
    return 0;
  }  catch (gromos::Exception e) {
    cerr << e.what() << endl;
    return 1;
  }
}

