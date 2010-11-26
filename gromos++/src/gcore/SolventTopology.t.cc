#include "SolventTopology.h"
#include "AtomTopology.h"
#include "Constraint.h"
#include <iostream>

using namespace gcore;
using namespace std;

int main() {
  SolventTopology mt;
  AtomTopology at;
  at.setIac(3);
  mt.addAtom(at);
  at.setIac(4);
  mt.addAtom(at);

  cout << "AtomTopology IAC: ";
  for (int i = 0; i < mt.numAtoms(); ++i)
    cout << mt.atom(i).iac() << ' ';
  cout << endl;

  Constraint b(1, 2);
  mt.addConstraint(b);
  b = Constraint(2, 3);
  b.setDist(0.1);

  mt.addConstraint(b);

  cout << "Constraints: ";
  ConstraintIterator bi(mt);
  for (; bi; ++bi)
    cout << '(' << bi()[0] << ',' << bi()[1] << ") " << bi().dist();
  cout << endl;

  mt.setSolvName("SPC");
  cout << "Solvname: ";
  cout << mt.solvName();
  cout << endl;

  return 0;
}
