#include <cassert>
#include <cstdlib>
#include "VirtualAtom.cc"
#include "../gcore/System.h"
#include "../gio/InTopology.h"
#include "../gio/InG96.h"
#include <string>
#include <iostream>

using namespace gcore;
using namespace gio;
using namespace utils;

using namespace std;

int debug_level = 0;

ostream & operator<<(ostream &os, const gmath::Vec &v) {
  os << v[0] << ' ' << v[1] << ' ' << v[2];
  return os;
}

int main(int argc, char *argv[]) {
  if (argc != 6) {
    cerr << "Usage: " << argv[0] << " <topology> <coordinates> <mol nr.> <atom nr.> <type>" << endl;
    exit(1);
  }
  string top = argv[1];
  string coord = argv[2];

  int mol = atoi(argv[3]) - 1;
  int atom = atoi(argv[4]) - 1;
  VirtualAtom::virtual_type type = VirtualAtom::virtual_type(atoi(argv[5]));
  try {

    InTopology it(top);
    System sys = it.system();

    InG96 ic(coord);
    ic >> sys;

    VirtualAtom virt0(sys, mol, atom, type, 0);

    cout << "Virtual atom created with configuration\n";
    for (int i = 0; i < 4; ++i)
      cout << virt0.conf().atom(i) << ' ';
    cout << virt0.type() << " " << 0 << endl;
    cout << "Position calculated as: ";
    cout << virt0.pos() << endl;

    VirtualAtom virt1(sys, mol, atom, type, 1);

    cout << "Virtual atom created with configuration\n";
    for (int i = 0; i < 4; ++i)
      cout << virt1.conf().atom(i) << ' ';
    cout << virt1.type() << " " << 1 << endl;

    cout << "Position calculated as: ";
    cout << virt1.pos() << endl;

  } catch (gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }

  return 0;
}
