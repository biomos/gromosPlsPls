// fit_TranslationalFit.t.cc

#include <cassert>
#include <cstdlib>
#include "TranslationalFit.h"
#include "PositionUtils.h"
#include "Reference.h"
#include "../gio/InG96.h"
#include "../gcore/System.h"
#include "../gio/InTopology.h"
#include "../gio/OutG96.h"
#include "../gmath/Vec.h"
#include <string>
#include <iostream>

using namespace gcore;
using namespace gio;
using namespace fit;
using namespace gmath;

using namespace std;

ostream & operator<<(ostream &os, const Vec &v) {
  os << v[0] << ' ' << v[1] << ' ' << v[2];
  return os;
}

int main(int argc, char *argv[]) {
  if (argc != 4) {
    cerr << "Usage: " + string(argv[0]) + " <Topology> <RefG96inates> <G96inates>\n";
    exit(1);
  }
  string top = argv[1];
  string refFile = argv[2];
  string file = argv[3];
  // read Simulation data

  InTopology it(top);

  System refsys(it.system());

  InG96 ic(refFile);
  ic >> refsys;
  ic.close();

  System sys(refsys);

  Reference ref(&refsys);

  ref.addClass(0, "CA");
  ref.addClass(0, "N");
  ref.addClass(0, "C");

  TranslationalFit tf(&ref, fit::cog);

  //  cout << "COM: " << tf.com() << endl;
  // cout << "COG: " << tf.cog() << endl;

  ic.open(file);

  OutG96 oc;
  oc.open(cout);
  oc.writeTitle(ic.title());

  cout << "Fitting to COG\n";
  //  while(!ic.eof()){
  //    ic >> sys;
  //    tf.fitToCog(&sys);
  //    oc << sys;
  //  }
  ic.close();
  oc.close();
  return 0;
}
