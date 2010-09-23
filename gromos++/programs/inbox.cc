/**
 * @file inbox.cc
 * Put atoms into box according to pbc.
 */

/**
 * @page contrib Contrib Program Documentation
 *
 * @anchor inbox
 * @section inbox put atoms into the box
 * @author @ref ns
 * @date 11-03-2009
 *
 * Program inbox puts the atoms into the positive computational box
 * according to the periodic boundary conditions. It can be used to visualize
 * the computational box in a crystal simulation. The connectivity and gathering
 * of chargegroups is ignored, thus the chargegroups (and solvent molecules)
 * won't be gathered after application of this program.
 *
 * Using the \@atoms argument one can specify the atoms which are put into the
 * box. All other atoms are not affected by the program. By default all atoms
 * are put into the box.
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;boundary type&gt; </td></tr>
 * <tr><td> \@traj</td><td>&lt;a trajectory&gt; </td></tr>
 * <tr><td>[\@atoms</td><td>&lt;@ref AtomSpecifier "atoms"&gt;]</td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
 inbox
    @topo             ex.top
    @pbc              r
    @atoms            a:a
    @traj             ex.tr
 @endverbatim
 *
 * <hr>
 */

#include <cassert>
#include <cmath>
#include <map>
#include <map>
#include <map>
#include <map>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <set>
#include <algorithm>
#include <iterator>

#include "../src/gcore/Box.h"
#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/gio/InG96.h"
#include "../src/gio/OutG96.h"
#include "../src/gcore/System.h"
#include "../src/gio/InTopology.h"
#include "../src/bound/Boundary.h"
#include "../src/gmath/Vec.h"
#include "../src/utils/AtomSpecifier.h"

using namespace gcore;
using namespace gio;
using namespace bound;
using namespace args;
using namespace gmath;
using namespace utils;
using namespace std;

int main(int argc, char **argv) {
  Argument_List knowns;
  knowns << "topo" << "pbc" << "atoms" << "traj";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo           <molecular topology file>\n";
  usage += "\t@pbc            <boundary type> [<gathermethod>]\n";
  usage += "\t[@atoms         <atoms>]\n";
  usage += "\t@traj           <trajectory files>\n";


  try {
    Arguments args(argc, argv, knowns, usage);

    InTopology it(args["topo"]);
    System sys(it.system());

    Boundary *pbc = BoundaryParser::boundary(sys, args);

    AtomSpecifier atoms(sys);
    {
      Arguments::const_iterator to = args.upper_bound("atoms"),
              iter = args.lower_bound("atoms");
      if (iter == to) {
        // add all atoms
        atoms.addSpecifier("a:a");
        atoms.addSpecifier("s:a");
      } else {
        for (; iter != to; iter++)
          atoms.addSpecifier(iter->second);
      }
    }

    // loop over all trajectories
    InG96 ic;
    Arguments::const_iterator iter = args.lower_bound("traj"),
            to = args.upper_bound("traj");

    if (iter == to)
      throw gromos::Exception(argv[0], "no trajetctory given.");

    // prepare the output
    OutG96 oc(cout);
    oc.writeTitle("Trajectory put in the box");

    for (;iter != to; ++iter) {
      // open file
      ic.open((iter->second).c_str());

      // loop over single trajectory
      while (!ic.eof()) {
        ic.select("ALL");
        ic >> sys;

        if (!sys.hasBox)
          throw gromos::Exception(argv[0], "no box in frame.");

        // get the centre of the box
        const Vec centre = (sys.box().K() * 0.5) + (sys.box().L() * 0.5) +
                (sys.box().M() * 0.5);

        // put atom into positive box
        for(int i = 0; i < atoms.size(); ++i) {
          atoms.pos(i) = atoms.pos(i) - pbc->nearestImage(atoms.pos(i), centre, sys.box()) + centre;
        }

        // write out the frame
        oc << sys;
      } // while frames
      ic.close();
    } // for traj

  } catch (const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}

