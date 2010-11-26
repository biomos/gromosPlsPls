#include <set>
#include <vector>
#include "LJException.h"
#include "Exclusion.h"
#include "MoleculeTopology.h"
#include "BbSolute.h"
#include "AtomTopology.h"
#include "Bond.h"
#include "Angle.h"
#include <iostream>

using namespace gcore;
using namespace std;

int main() {
  BbSolute mt;
  AtomTopology at;
  Exclusion e;
  e.insert(2);
  mt.addPexcl(e);
  at.setIac(3);
  mt.addAtom(at);
  at.setIac(4);
  mt.addAtom(at);

  cout << "AtomTopology IAC: ";
  for (int i = 0; i < mt.numAtoms(); ++i)
    cout << mt.atom(i).iac() << ' ';
  cout << endl;
  for (int i = 0; i < mt.numPexcl(); ++i) {

    cout << mt.pexcl(i).size() << " :";
    for (int j = 0; j < mt.pexcl(i).size(); ++j)
      cout << mt.pexcl(i).atom(j) << " ";

    cout << endl;
  }

  Bond b(1, 2);
  mt.addBond(b);
  b = Bond(2, 3);
  mt.addBond(b);

  cout << "Bonds: ";
  BondIterator bi(mt);
  for (; bi; ++bi)
    cout << '(' << bi()[0] << ',' << bi()[1] << ") ";
  cout << endl;

  Angle ang(1, 2, 3);
  mt.addAngle(ang);
  ang = Angle(4, 5, 6);
  mt.addAngle(ang);

  cout << "Angles: ";
  AngleIterator ai(mt);
  for (; ai; ++ai)
    cout << '(' << ai()[0] << ' ' << ai()[1] << ' ' << ai()[2] << ") ";
  cout << endl;

  mt.setResName("ASP");
  cout << "Resnames: ";
  for (int i = 0; i < 5; ++i)
    cout << i << ' ' << mt.resName() << ' ';
  cout << endl;

  return 0;
}
