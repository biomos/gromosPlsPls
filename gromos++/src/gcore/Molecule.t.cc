#include <cassert>
#include <set>

#include "LJException.h"
#include "Molecule.h"
#include "MoleculeTopology.h"
#include "AtomTopology.h"
#include "Bond.h"
#include "Angle.h"
#include "../gmath/Vec.h"

#include <iostream>

using namespace std;

using namespace gcore;
using namespace std;
using namespace gmath;

ostream & operator<<(ostream &o, Vec &v) {
  o << '(' << v[0] << ',' << v[1] << ',' << v[2] << ')';
  return o;
}

int main() {
  MoleculeTopology mt;
  AtomTopology at;
  at.setIac(3);
  mt.addAtom(at);
  at.setIac(4);
  mt.addAtom(at);

  Bond b(1, 2);
  mt.addBond(b);
  b = Bond(2, 3);
  mt.addBond(b);

  Angle ang(1, 2, 3);
  mt.addAngle(ang);
  ang = Angle(4, 5, 6);
  mt.addAngle(ang);

  mt.setResName(4, "ASP");

  Molecule mol(mt);

  cout << "Number of atoms: " << mol.numAtoms() << endl;

  mol.pos(0) = Vec(1, 2, 3);
  mol.pos(1) = Vec(4, 5, 6);

  cout << "Pos0: " << mol.pos(0) << endl;
  cout << "Pos1: " << mol.pos(1) << endl;

  cout << "Bonds: ";
  BondIterator bi(mol.topology());
  for (; bi; ++bi)
    cout << '(' << bi()[0] << ',' << bi()[1] << ") ";
  cout << endl;

  return 0;
}




