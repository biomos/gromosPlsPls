#include <cassert>
#include <set>
#include <iostream>
#include <string>
#include <map>
#include <cstdlib>
#include "../gcore/BuildingBlock.h"
#include "../gcore/AtomTopology.h"
#include "../gcore/Bond.h"
#include "../gcore/Angle.h"
#include "../gcore/Dihedral.h"
#include "../gcore/Improper.h"
#include "../gio/InBuildingBlock.h"
#include "FfExpert.h"

using namespace gcore;
using namespace gio;
using namespace utils;

using namespace std;

int debug_level = 0;

int main(int argc, char *argv[]) {
  if (argc != 2) {
    cerr << "Usage: " + string(argv[0]) + " <BuildingBlock>\n";
    exit(1);
  }
  try {
    InBuildingBlock ibb(argv[1]);

    FfExpert exp(ibb.building());

    cout << "succes" << endl;
    vector<FfExpert::counter> v;
    exp.name2iac("C", v);

    cout << "we have " << v.size() << " atoms with C" << endl;
    for (unsigned int i = 0; i < v.size(); i++) {
      cout << v[i].type << " " << v[i].occurence << endl;
    }
    exp.iac2mass(12, v);

    cout << "we have " << v.size() << " masses for iac=12" << endl;
    for (unsigned int i = 0; i < v.size(); i++) {
      cout << v[i].type << " " << v[i].occurence << endl;
    }
    int b = 2;

    exp.iac2charge(b, v);
    cout << "we have " << v.size() << " charges for iac=" << b << endl;
    for (unsigned int i = 0; i < v.size(); i++) {
      cout << exp.charge(v[i].type) << " " << v[i].occurence << endl;
    }


  }  catch (gromos::Exception e) {
    cerr << e.what() << endl;
    return 1;
  }
}
