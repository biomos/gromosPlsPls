/**
 * @file sad.cc
 * Calculates the solute averaged distance
 */

/**
 * @page programs Program Documentation
 *
 * @anchor sad
 * @section sad Calculates time serie of solute averaged distances
 * @author @ref nb
 * @date 22. 11. 2004
 *
 * Program sad calculates the solute averaged distance for fg and cg solvent.
 * If the force constant and cutoff radius is supplied, the energy is
 * calculated as well.
 * 
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;boundary type&gt; [&lt;gathermethod&gt;] </td></tr>
 * <tr><td> [\@time</td><td>&lt;@ref utils::Time "time and dt"&gt;] </td></tr>
 * <tr><td> \@solute</td><td>&lt;@ref AtomSpecifier "solute atoms"&gt; </td></tr>
 * <tr><td> \@fgsolv</td><td> &lt;@ref AtomSpecifier "fg solvent atoms"&gt;</td></tr>
 * <tr><td> \@cgsolv</td><td> &lt;@ref AtomSpecifier "cg solvent atoms"&gt;</td></tr>
 * <tr><td> [\@energy</td><td>&lt;fg force constant&gt; &lt;fg force cut off&gt;  
 *                            &lt;cg force constant&gt; &lt;cg force cut off&gt;]</td></tr>
 * <tr><td> \@traj</td><td>&lt;trajectory files&gt; </td></tr>
 * </table>
 * 
 * <hr>
 */
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cassert>
#include <string>

#include "../src/args/Arguments.h"
#include "../src/gio/InG96.h"
#include "../src/gcore/System.h"
#include "../src/gio/InTopology.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/utils/groTime.h"

#include "../src/utils/SoluteAverageDistance.h"

int main(int argc, char **argv) {

  args::Argument_List knowns;
  knowns << "topo" << "pbc" << "time" << "solute" << "traj"
          << "fgsolv" << "cgsolv" << "energy";

  std::string usage = "# " + std::string(argv[0]) + "\n";
  usage += "\t@topo       <molecular topology file>\n";
  usage += "\t@pbc        <boundary type> [<gathermethod>]\n";
  usage += "\t[@time       <time and dt>]\n";
  usage += "\t@solute     <solute atoms>\n";
  usage += "\t@fgsolv     <fg solvent atoms>\n";
  usage += "\t@cgsolv     <cg solvent atoms>\n";
  usage += "\t[@energy     <fg force constant> <fg cutoff> <cg force constant> <cg cutoff>]\n";
  usage += "\t@traj       <trajectory files>\n";

  std::ofstream eneout;
  int returnCode = 0;

  try {
    args::Arguments args(argc, argv, knowns, usage);

    // get the @time argument
    utils::Time time(args);

    //  read topology
    args.check("topo", 1);
    gio::InTopology it(args["topo"]);

    gcore::System sys(it.system());

    utils::SoluteAverageDistance sad(sys, args); // :-)

    // define input coordinate
    gio::InG96 ic;

    std::cout << sad.title();
    if (sad.withEnergy()) {
      std::string eneOutFile = "cg_ene.out";
      eneout.open(eneOutFile.c_str());
      eneout << sad.title();
      std::cerr << "# Energy parameter given. Will print energies to "
              << eneOutFile << std::endl;
    } else {
      std::cerr << "# No energy parameters given, but don't panic!\n";
    }

    // loop over all trajectories
    for (args::Arguments::const_iterator iter = args.lower_bound("traj"),
            to = args.upper_bound("traj");
            iter != to; ++iter) {

      // open file
      ic.open((iter->second).c_str());
      ic.select("ALL");

      // loop over single trajectory
      while (!ic.eof()) {
        ic >> sys >> time;
        if (ic.stride_eof())
          break;

        sad.calculate();

        std::cout << time << " \t\t";
        std::cout << sad;
        std::cout << std::endl;
        
        if (sad.withEnergy()){
          eneout << time << " \t\t";
          sad.energies(eneout);
          eneout << std::endl;
        }
      }


      ic.close();
    }

  } catch (const gromos::Exception &e) {
    std::cerr << e.what() << std::endl;
    returnCode = 1;
  }
  if (eneout.is_open()) {
    eneout.close();
  }
  return returnCode;
}