/**
 * @file hbond.cc
 * Monitors the occurrence of hydrogen bonds
 */

/**
 * @page programs Program Documentation
 *
 * @anchor hbond
 * @section hbond monitors the occurrence of hydrogen bonds
 * @author @ref mk
 * @date 9-8-2006
 *
 * Program hbond monitors the occurrence of hydrogen bonds over a molecular 
 * trajectory file. It can monitor conventional hydrogen bonds, as well as 
 * three-centered hydrogen bonds through geometric criteria.
 *
 * A hydrogen bond is considered to be present if the distance between a 
 * hydrogen atom, H, connected to a donor atom D, is within a user specified 
 * distance (typically 0.25 nm) from an acceptor atom A and the D-H-A angle is 
 * larger than another user specified value (typically 135 degrees). Occurrences
 * of three centered hydrogen bonds are defined for a donor atom D, hydrogen
 * atom H and two acceptor atoms A1 and A2 if (i) the distances H-A1 and H-A2 
 * are within a user specified value (typically 0.27 nm); (ii) the angles 
 * D-H-A1 and D-H-A2 are larger than a second user specified value (typically 
 * 90 degrees); (iii) the sum of the angles D-H-A1, D-H-A2 and A1-H-A2 is larger
 * than a third user specified value (typically 340 degrees); and (iv) the 
 * dihedral angle defined by the planes through the atoms D-A1-A2 and H-A1-A2 
 * is smaller than a fourth user specified value (typically 15 degrees).
 *
 * The user can specify two groups of atoms (A and B) between which the 
 * hydrogen bonds are to be monitored. If hydrogen bond donors and acceptors 
 * are not explicitly specified, these can be filtered based on their masses, 
 * as can be specified in a so-called "massfile". If a reference structure is 
 * given, only hydrogen bonds that are observed in the reference structure will
 * be monitored. 
 * 
 * The program calculates average angles, distances and occurrences for all 
 * observed hydrogen bonds over the trajectories and prints out a time series 
 * of the observed hydrogen bonds.
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;boundary type&gt; </td></tr>
 * <tr><td> [\@time</td><td>&lt;@ref utils::Time "time and dt"&gt;] </td></tr>
 * <tr><td> \@DonorAtomsA</td><td>&lt;@ref AtomSpecifier "atoms"&gt; </td></tr>
 * <tr><td> \@AcceptorAtomsA</td><td>&lt;@ref AtomSpecifier "atoms"&gt; </td></tr>
 * <tr><td> \@DonorAtomsB</td><td>&lt;@ref AtomSpecifier "atoms"&gt; </td></tr>
 * <tr><td> \@AcceptorAtomsB</td><td>&lt;@ref AtomSpecifier "atoms"&gt; </td></tr>
 * <tr><td> \@Hbparas</td><td>&lt;distance [nm] and angle; default: 0.25, 135&gt; </td></tr>
 * <tr><td> [\@threecenter</td><td>&lt;distances [nm]&gt; &lt;angles&gt; &lt;sum&gt; &lt;dihedral&gt]; </td></tr>
 * <tr><td> [\@ref</td><td>&lt;reference coordinates for native H-bonds&gt;] </td></tr>
 * <tr><td> [\@massfile</td><td>&lt;massfile&gt;] </td></tr>
 * <tr><td> \@traj</td><td>&lt;trajectory files&gt; </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  hbond
    @topo             ex.top
    @pbc              r
    @time             0 1
    @DonorAtomsA      1:a
    @AcceptorAtomsA   1:a
    @DonorAtomsB      s:a
    @AcceptorAtomsB   s:a
    @Hbparas          0.25 135
    @threecenter      0.27 90 340 15
    @massfile         ../data/hbond.massfile
    @ref              exref.coo
    @traj             ex.tr
 @endverbatim
 *
 * <hr>
 */

#include <cassert>
#include "../src/args/Arguments.h"
#include "../src/gio/InG96.h"
#include "../src/gcore/System.h"
#include "../src/gio/InTopology.h"
#include "../src/utils/Hbond.h"
#include "../src/utils/groTime.h"

using namespace gcore;
using namespace gio;
using namespace args;
using namespace utils;
using namespace std;

int main(int argc, char** argv) {
  Argument_List knowns;
  knowns << "topo" << "pbc" << "ref" << "DonorAtomsA" << "AcceptorAtomsA"
          << "DonorAtomsB" << "AcceptorAtomsB" << "Hbparas" << "threecenter"
          << "time" << "massfile" << "traj";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo         <molecular topology file>\n";
  usage += "\t@pbc            <boundary type> [<gathermethod>]\n";
  usage += "\t[@time          <time and dt>]\n";
  usage += "\t@DonorAtomsA    <atoms>\n";
  usage += "\t@AcceptorAtomsA <atoms>\n";
  usage += "\t@DonorAtomsB    <atoms>\n";
  usage += "\t@AcceptorAtomsB <atoms>\n";
  usage += "\t@Hbparas        <distance [nm] and angle; default: 0.25, 135>\n";
  usage += "\t[@threecenter   <distances [nm]> <angles> <sum> <dihedral>]\n";
  usage += "\t[@massfile      <massfile>]\n";
  usage += "\t[@ref           <reference coordinates for native H-bonds>]\n";
  usage += "\t@traj           <trajectory files>\n";


  try {
    Arguments args(argc, argv, knowns, usage);
    InTopology it(args["topo"]);
    System sys(it.system());
    Time time(args);

    bool hbond3c = false;
    // get the paras
    vector<double> v_hbparas2c = args.getValues<double>("Hbparas", 2, false,
            Arguments::Default<double>() << 0.25 << 135.0);
    vector<double> v_hbparas3c;
    v_hbparas3c.resize(4);
    if (args.count("threecenter") >= 0) {
      hbond3c = true;
      v_hbparas3c = args.getValues<double>("threecenter", 4, false,
              Arguments::Default<double>() << 0.27 << 90.0 << 340.0 << 15.0);
    }
    HBPara2c hbparas2c = HB::mk_hb2c_paras(v_hbparas2c);
    HBPara3c hbparas3c = HB::mk_hb3c_paras(v_hbparas3c);
    HB hb(sys, args, hbparas2c, hbparas3c);
    // initialize the calculation
    hb.init();
    // do native?
    if (args.count("ref") > 0) {
      InG96 ic(args["ref"]);
      ic.select("ALL");
      ic >> sys;
      ic.close();
      // calculate the hb.
      hb.calc();
      // and clear the statistics and reset the time
      hb.clear();
    }
    InG96 ic;
    // loop over all trajectories
    for (Arguments::const_iterator iter = args.lower_bound("traj"),
            to = args.upper_bound("traj"); iter != to; ++iter) {

      // open file
      ic.open((iter->second).c_str());
      ic.select("ALL");

      // loop over single trajectory

      while (!ic.eof()) {
        ic >> sys >> time;
        static int frame = 1;
        // get the number of atoms and break in case these numbers change from
        // one frame to another
        int numSoluAt = 0, numSolvAt = 0;
        static int numSoluAt_old = -1, numSolvAt_old = -1;
        int numSolu = sys.numMolecules();
        int numSolv = sys.numSolvents();
        for (int i = 0; i < numSolu; ++i) {
          numSoluAt += sys.mol(i).numAtoms();
        }
        for (int i = 0; i < numSolv; ++i) {
          numSolvAt += sys.sol(i).numAtoms();
        }
        if (numSoluAt_old != -1 && numSoluAt != numSoluAt_old) {
          stringstream msg;
          msg << "The number of solute atoms changed in " << iter->second.c_str() << ":\n"
                  << "             frame " << frame - 1 << ": " << numSoluAt_old << " solute atoms\n"
                  << "             frame " << frame << ": " << numSoluAt << " solute atoms\n"
                  << "       The calculation of hbond has been stopped therefore.";
          throw gromos::Exception("hbond", msg.str());
        }
        if (numSolvAt_old != -1 && numSolvAt != numSolvAt_old) {
          stringstream msg;
          msg << "The number of solvent atoms changed in " << iter->second.c_str() << ":\n"
                  << "             frame " << frame - 1 << ": " << numSolvAt_old << " solvent atoms\n"
                  << "             frame " << frame << ": " << numSolvAt << " solvent atoms\n"
                  << "       The calculation of hbond has been stopped therefore.";
          throw gromos::Exception("hbond", msg.str());
        }

        numSoluAt_old = numSoluAt;
        numSolvAt_old = numSolvAt;
        hb.settime(time.time());
        hb.calc();
        frame++;
      }
      ic.close();
    }
    hb.printstatistics();
  } catch (const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}