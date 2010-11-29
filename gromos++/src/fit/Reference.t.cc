// fit_Reference.t.cc

#include <cassert>
#include <cstdlib>
#include "Reference.h"
#include "../gcore/System.h"
#include "../gcore/Molecule.h"
#include "../gio/InTopology.h"
#include <string>
#include <vector>
#include <iostream>
#include <algorithm>

using namespace gcore;
using namespace gio;
using namespace fit;

using namespace std;

ostream & operator<<(ostream &os, const vector<int>& v) {
  for (vector<int>::const_iterator i = v.begin(); i != v.end(); ++i)
    os << *i << " ";
  return os;
}

ostream & operator<<(ostream &os, const Reference &w) {
  for (int i = 0; i < w.sys().numMolecules(); ++i)
    for (int j = 0; j < w.sys().mol(i).numAtoms(); ++j)
      os << w.weight(i, j) << ' ';
  return os;
}

int main(int argc, char *argv[]) {
  if (argc != 2) {
    cerr << "Usage: " + string(argv[0]) + " <Topology>\n";
    exit(1);
  }
  string top = argv[1];
  // read Simulation data

  InTopology it(top);
  System bla(it.system());
  System *refSys = new System(it.system());
  Reference w(refSys);


  //get atom positions from name

  int mole = 0;
  const string t = "CA";
  vector<int> polist;
  w.makePosList(bla, mole, t, polist);
  int tmp = polist.size();

  cout << "Number of " << t << " atoms: " << tmp << endl;
  cout << "Respective Positions:" << endl;
  std::cout << polist << endl;

  const string tt = "N";
  w.makePosList(bla, mole, tt, polist);
  tmp = polist.size() - tmp;

  cout << "Number of " << tt << " atoms: " << tmp << endl;
  cout << "Respective Positions:" << endl;
  cout << "1st element " << polist[0] << endl;
  std::cout << polist << endl;

  //sort the positions

  std::sort(polist.begin(), polist.end());
  cout << "Sorted positions:" << endl;
  std::cout << polist << endl;

  // do it for keyword 'ALL'

  polist.erase(polist.begin(), polist.end());
  cout << polist.size() << polist << endl;
  const string ttt = "ALL";
  w.makePosList(bla, mole, ttt, polist);
  tmp = polist.size();

  cout << "Number of " << ttt << " atoms: " << tmp << endl;
  cout << "Respective Positions:" << endl;
  std::cout << polist << endl;



  cout << w << endl;

  w.addClass(0, "CA");
  cout << "All CA:\n";
  cout << w << endl;
  w.addClass(0, "N");
  cout << "All CA and N:\n";
  cout << w << endl;
  w.addAtom(0, 0);
  cout << "All CA and N and atom 0:\n";
  cout << w << endl;

  if (it.system().numMolecules() > 1) {
    w.addClass(1, "FE");
    cout << "Added a FE from mol 1:\n"
            << w << endl;
  }

  w.setWeight(0, 1, 0.2);
  w.normalise();
  cout << "Added 0.2 to atom 1 and normalised...\n";
  cout << w << endl;

  delete refSys;
  return 0;
}
