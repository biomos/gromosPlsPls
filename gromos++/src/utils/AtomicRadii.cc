#include <cassert>
#include <iostream>
#include <sstream>
#include <cstdio>
#include <string>
#include <cmath>
#include <set>
#include <map>

#include "../gcore/System.h"
#include "../gcore/Molecule.h"
#include "../gcore/LJException.h"
#include "../gcore/MoleculeTopology.h"
#include "../gcore/Solvent.h"
#include "../gcore/SolventTopology.h"
#include "../gcore/AtomTopology.h"
#include "../gcore/LJType.h"
#include "../gcore/AtomPair.h"
#include "../gcore/GromosForceField.h"

#include "AtomicRadii.h"

using namespace gcore;

void utils::compute_atomic_radii_vdw(int probe_iac, double probe_radius, gcore::System & sys, const gcore::GromosForceField & gff) {
  static const double small = 1.0E-20;
  for (int m = 0; m < sys.numMolecules(); ++m) {
    for (int a = 0; a < sys.mol(m).topology().numAtoms(); ++a) {
      int atom_iac = sys.mol(m).topology().atom(a).iac();
      LJType lj(gff.ljType(AtomPair(atom_iac, probe_iac)));
      if (lj.c6() >= small) {
        sys.mol(m).topology().atom(a).setradius((exp(log(2.0 * lj.c12() / lj.c6()) / 6.0)) - probe_radius);
      } else {
        sys.mol(m).topology().atom(a).setradius(0.0);
      }
    }
  }
  for(int s = 0; s < sys.numSolvents(); ++s) {
    for(int a = 0; a < sys.sol(s).topology().numAtoms(); ++a) {
      int atom_iac = sys.sol(s).topology().atom(a).iac();
      LJType lj(gff.ljType(AtomPair(atom_iac, probe_iac)));
      if (lj.c6() >= small) {
        sys.sol(s).topology().atom(a).setradius((exp(log(2.0 * lj.c12() / lj.c6()) / 6.0)) - probe_radius);
      } else {
        sys.sol(s).topology().atom(a).setradius(0.0);
      }
    }
  }
}
