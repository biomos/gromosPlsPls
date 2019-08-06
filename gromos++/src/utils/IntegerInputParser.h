#include <cassert>
#include <vector>
#include <iomanip>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <set>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/fit/Reference.h"
#include "../src/gio/InG96.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/Box.h"
#include "../src/gio/Ginstream.h"
#include "../src/gio/InTopology.h"
#include "../src/gio/OutG96.h"
#include "../src/gio/OutG96S.h"
#include "../src/bound/Boundary.h"
#include "../src/gmath/Vec.h"

using namespace std;
using namespace gcore;
using namespace gio;
using namespace bound;
using namespace args;
using namespace utils;

namespace utils
{
  class IntegerInputParser: public set<int>
  {
  public:
    void addSpecifier(string const s, int maxnum);
  protected:
    void parse(string const s, int maxnum);
  };
}
