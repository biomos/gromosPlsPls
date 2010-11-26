#include <cassert>

#include "Solvent.h"
#include "SolventTopology.h"
#include "AtomTopology.h"
#include "Constraint.h"
#include "../gmath/Vec.h"
#include <iostream>

using namespace gcore;
using namespace std;
using namespace gmath;

using namespace std;

ostream & operator<<(ostream &o, Vec &v) {
  o << '(' << v[0] << ',' << v[1] << ',' << v[2] << ')';
  return o;
}

int main() {
  SolventTopology mt;
  AtomTopology at;
  at.setIac(3);
  mt.addAtom(at);
  at.setIac(4);
  mt.addAtom(at);

  Constraint b(1, 2);
  b.setDist(0.123);

  mt.addConstraint(b);
  b = Constraint(2, 3);
  mt.addConstraint(b);

  mt.setSolvName("SPC");

  Solvent solv(mt);

  cout << "Number of atoms in topo: " << solv.topology().numAtoms() << endl;
  cout << "Number of solvent coordinates: " << solv.numPos() << endl;

  solv.addPos(Vec(1, 2, 3));
  solv.addPos(Vec(4, 5, 6));
  solv.addPos(Vec(1, 2, 3));
  solv.addPos(Vec(4, 5, 6));
  solv.pos(2) = Vec(7, 8, 9);

  cout << "Number of solvent coordinates: " << solv.numPos() << endl;

  for (int i = 0; i < solv.numPos(); i++)
    cout << "Pos: " << i + 1 << " " << solv.pos(i) << endl;


  cout << "Constraints: ";
  ConstraintIterator bi(solv.topology());
  for (; bi; ++bi)
    cout << '(' << bi()[0] << ',' << bi()[1] << ") " << bi().dist();
  cout << endl;
  cout << "Number of solvent coordinates: " << solv.numPos() << endl;

  cout << "Number of solvent molecules: "
          << solv.numPos() / solv.topology().numAtoms() << endl;

  cout << "Setting number of Solv-Coords to 7" << endl;
  solv.setNumPos(7);
  cout << "Number of solvent Coords: ";
  cout << solv.numPos() << endl;

  return 0;
}
