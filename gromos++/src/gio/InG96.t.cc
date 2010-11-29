// gio_InG96.t.cc

#include <cassert>
#include <cstdlib>

#include "InG96.h"
#include "../gcore/System.h"
#include "InTopology.h"
#include "OutG96S.h"
#include <string>
#include <iostream>

using namespace std;
using namespace gcore;
using namespace gio;

using namespace std;

int main(int argc, char *argv[]) {
  if (argc != 3) {
    cerr << "Usage: " + string(argv[0]) + " <Topology> <Filename>\n";
    exit(1);
  }
  string top = argv[1];
  string file = argv[2];
  // read Simulation data

  InTopology it(top);
  System sys(it.system());
  cout << "read topology and created system" << endl;

  InG96 ic;
  ic.open(file);
  ic.select("ALL");
  OutG96S oc;

  oc.open(cout);
  oc.select("ALL");
  oc.writeTitle(ic.title());

  //while(!ic.eof()){
  ic >> sys;
  cout << "read in something " << endl;

  oc << sys;
  //}
  return 0;
}


