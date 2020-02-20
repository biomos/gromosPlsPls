#include <cassert>
#include <vector>
#include <iomanip>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <set>

#include "../args/Arguments.h"
#include "../args/BoundaryParser.h"
#include "../args/GatherParser.h"
#include "../fit/Reference.h"
#include "../gio/InG96.h"
#include "../gcore/System.h"
#include "../gcore/Molecule.h"
#include "../gcore/Box.h"
#include "../gio/Ginstream.h"
#include "../gio/InTopology.h"
#include "../gio/OutG96.h"
#include "../gio/OutG96S.h"
#include "../bound/Boundary.h"
#include "../gmath/Vec.h"

using namespace std;
using namespace gcore;
using namespace gio;
using namespace bound;
using namespace args;
using namespace utils;

// A class that can take comma-separated identifiers, ranges, or a combination (e.g. 1,4,7-9)
// to parse user input in Integer format for various uses (e.g. selection of energy groups or similar)


namespace utils
{
  class IntegerInputParser: public set<int>
  {
  public:
   // Method to add the Integer numbers in a string (e.g. from user input)
   // maxnum is the highest number allowed (smallest is hardcoded to zero)
    void addSpecifier(string const s, int maxnum);
  protected:
    void parse(string const s, int maxnum);
  };
}
