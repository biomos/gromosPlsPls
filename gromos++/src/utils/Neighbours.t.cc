// utils_Neighbours.t.cc

#include <cassert>
#include <cstdlib>
#include "Neighbours.h"
#include "../gcore/System.h"
#include "../gio/InTopology.h"
#include <string>
#include <iostream>

using namespace gcore;
using namespace gio;
using namespace utils;

using namespace std;

int main(int argc, char *argv[]) {
  if (argc != 4) {
    cerr << "Usage: " << argv[0] << " <topology> <mol nr.> <atom nr.>" << endl;
    exit(1);
  }
  string top = argv[1];
  int mol = atoi(argv[2]) - 1;
  int nr = atoi(argv[3]) - 1;

  InTopology it(top);
  System sys = it.system();

  Neighbours neigh(sys, mol, nr);
  cout << "Neighbours of atom " << nr + 1 << " of molecule " << mol + 1 << ": ";
  for (Neighbours::iterator iter = neigh.begin(); iter != neigh.end(); ++iter)
    cout << *iter + 1 << ' ';
  cout << endl;
  return 0;
}
