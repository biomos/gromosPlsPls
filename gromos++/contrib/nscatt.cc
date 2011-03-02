/**
 * @file nscatt.cc
 */

/**
 * @page programs Program Documentation
 *
 * @anchor nscatt
 * @section nscatt short description
 * @author A. Eichenberger
 * @date March 2011
 *
 * Description...
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * ...
 * <tr><td> \@traj</td><td>&lt;trajectory files&gt; </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
 nscatt
    @topo   ex.top
    @traj   ex.tr
 @endverbatim
 *
 * <hr>
 */

#include <cassert>
#include <iostream>
#include <cstdlib>
#include <sstream>

#include "../src/args/Arguments.h"
#include "../src/gio/InTopology.h"
#include "../src/gcore/System.h"
#include "../src/utils/NeutronScattering.h"


using namespace std;
using namespace args;
using namespace gcore;
using namespace gio;
using namespace utils;

int main(int argc, char **argv) {

  Argument_List knowns;
  knowns << "topo" << "atoms" << "grid" << "traj";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo   <molecular topology file>\n";
  usage += "\t@centre   <AtomSpecifier: atoms to be considered as centre atoms>\n";
  usage += "\t@with     <AtomSpecifier: atoms to be considered as with atoms>\n";
  usage += "\t[@grid    <number of data points>]\n";


  try {
    Arguments args(argc, argv, knowns, usage);

    // read topology
    args.check("topo", 1);
    InTopology it(args["topo"]);
    System sys(it.system());

    if(args.count("traj") < 1) {
      throw gromos::Exception(argv[0], "no trajectory file(s) specified");
    }
    NS ns(&sys, args.lower_bound("traj"), args.upper_bound("traj"));

    // set the centre and with atoms
    if(args.count("atoms") < 1) {
      throw gromos::Exception(argv[0], "no atoms defined (@atoms)");
    }
    {
      Arguments::const_iterator iter = args.lower_bound("atoms");
      Arguments::const_iterator to = args.upper_bound("atoms");
      for (; iter != to; iter++) {
        ns.addAtoms(iter->second.c_str());
      }
    }

    // this is the sequence it HAS to be done:
    //   1) get the number of combinations (also resetting all the vector lengths
    //   2) set the system (for all subvectors too
    //   3) set the atoms to the AtomSpecifiers (for all subvectors)
    cerr << "Getting the atom combinations...";
    ns.getCombinations();
    cerr << "done\n";
    cerr << "Setting the system variable...";
    ns.setSystem(&sys);
    cerr << "done\n";
    cerr << "Setting the trajectories...";
    ns.setTrajectories(args.lower_bound("traj"), args.upper_bound("traj"));
    cerr << "done\n";
    cerr << "Setting the RDF atoms...";
    ns.setRDFatoms();
    cerr << "done\n";
    cerr << "Getting the weights of the partial structure facors...";
    ns.getWeights();
    cerr << "done\n";

    // set the grid number to the specified integer number, if there is an @ grid flag
    if(args.count("grid") > 0) {
      stringstream ss;
      ss << args["grid"];
      int grid;
      ss >> grid;
      if(ss.fail() || ss.bad()) {
        stringstream msg;
        msg << "could not convert " << args["grid"] << " into an integer number";
        throw gromos::Exception(argv[0], msg.str());
      }
      ns.setGrid(grid);
    }

    cerr << "Checking the NeutronScattering class...";
    ns.check();
    cerr << "done\n";

    ns.print(cerr);
    cerr << "Calculating the RDFs...";
    ns.calcRDFsInterAll();
    cerr << "done\n";
    ns.printRDFs(cout);

  } catch (const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}

