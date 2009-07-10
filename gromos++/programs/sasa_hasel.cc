/**
 * @file sasa_hasel.cc
 * compute sasa using hasel formula
 */

/**
 * @page programs Program Documentation
 *
 * @anchor sasa_hasel
 * @section sasa_hasel compute sasa using hasel formula
 * @author @ref ja
 * @date 23. 4. 2009
 *
 * Program sasa_hasel computes the solvent-accessible surface area (sasa)
 * of all atoms in the solute part of the molecular system according to the
 * method of Hasel et al. [Tetra. Comput. Method., 1, 103-116, (1988)]. This is
 * the same method implemented in the SASA/VOL implicit solvent model. If a
 * single conformation is given, either the atomic sasa values or the total sasa,
 * along with the hydrophilic and hydrophobic contributions (defined by the sign
 * of the sigma values given in the sasaspec file) may be printed. If multiple
 * conformations are given, the averaged totals, the averaged atomic sasa values,
 * or a time-series of the total sasa values may be printed.

 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;boundary type&gt; [&lt;gather method&gt;] </td></tr>
 * <tr><td> [\@time</td><td>&lt;@ref utils::Time "time and dt" (optional and only if time-series)&gt;] </td></tr>
 * <tr><td> [\@timeseries</td><td>&lt;if you want the time-series as well as the average&gt;] </td></tr>
 * <tr><td> [\@timespec</td><td>&lt;timepoints at which to compute the sasa: ALL (default), EVERY or SPEC (if time-series)&gt;] </td></tr>
 * <tr><td> [\@timepts</td><td>&lt;timepoints at which to compute the sasa (if time-series and timespec EVERY or SPEC)&gt;] </td></tr>
 * <tr><td> [\@atomic</td><td>&lt;print atomic sasa (only if not time-series)&gt;] </td></tr>
 * <tr><td> \@sasaspec</td><td>&lt;sasa specification library file&gt; </td></tr>
 * <tr><td> [\@radius</td><td>&lt;radius of water molecule (default: 0.14 nm)&gt;] </td></tr>
 * <tr><td> \@traj</td><td>&lt;trajectory file(s)&gt; </td></tr>
 * </table>
 *
 * Example:
 * @verbatim
   atominfo
     @topo       ex.top
     @pbc        v
     @timeseries
     @timespec   EVERY
     @timepts    100
     @sasaspec   sasaspec45b3.lib
     @traj       ex.trj
   @endverbatim

 * <hr>
 */

#include <cassert>
#include <map>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <set>
#include <fstream>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/bound/Boundary.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gcore/Exclusion.h"
#include "../src/gcore/Bond.h"
#include "../src/gcore/Angle.h"
#include "../src/gcore/Improper.h"
#include "../src/gcore/Dihedral.h"
#include "../src/gio/InG96.h"
#include "../src/gio/InTopology.h"
#include "../src/gio/Ginstream.h"
#include "../src/gmath/Vec.h"
#include "../src/gmath/Physics.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/utils/Time.h"
#include "../src/utils/Neighbours.h"

using namespace std;
using namespace args;
using namespace bound;
using namespace gcore;
using namespace gio;
using namespace gmath;
using namespace utils;

bool compute_sasa(int i, std::string const & timespec, vector<int> const & timepts,
        unsigned int & timesWritten, bool & done);

struct sasa_parameter {
  double radius;
  double probability;
  double sigma;
};

int main(int argc, char **argv) {

  Argument_List knowns;
  knowns << "topo" << "pbc" << "time" << "timeseries" << "timespec"
          << "timepts" << "atomic" << "sasa_spec" << "radius" << "traj";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo        <molecular topology file>\n";
  usage += "\t@pbc         <boundary type> [<gather method>]\n";
  usage += "\t[@time       <time and dt>]\n";
  usage += "\t[@timeseries <if you want the time-series as well as the average>]\n";
  usage += "\t[@timespec   <timepoints at which to compute the SASA: ALL (default), EVERY or SPEC>])\n";
  usage += "\t[@timepts    <timepoints at which to compute the SASA (if timespec EVERY or SPEC)>]\n";
  usage += "\t[@atomic     <print atomic sasa (only if not time-series)>]\n";
  usage += "\t@sasa_spec   <sasa specification file>\n";
  usage += "\t[@radius     <radius of water molecule> (default: 0.14 nm)]\n";
  usage += "\t@traj        <trajectory file(s)>\n";

  try {
    Arguments args(argc, argv, knowns, usage);

    // read topology
    InTopology it(args["topo"]);
    gcore::System sys(it.system());

    // parse boundary conditions
    Boundary *pbc = BoundaryParser::boundary(sys, args);
    // parse gather method
    Boundary::MemPtr gathmethod = args::GatherParser::parse(args);

    // get time
    Time time(args);

    // parse timespec
    string timespec = "ALL";
    vector<int> timepts;
    if (args.count("timespec") > 0) {
      timespec = args["timespec"];
      if (timespec != "ALL" && timespec != "EVERY" && timespec != "SPEC")
        throw gromos::Exception("sasa_hasel",
              "timespec format " + timespec + " unknown.\n");
      if (timespec == "EVERY" || timespec == "SPEC") {
        for (Arguments::const_iterator it = args.lower_bound("timepts");
                it != args.upper_bound("timepts"); ++it) {
          int bla = atoi(it->second.c_str());
          timepts.push_back(bla);
        }
        if (timepts.size() == 0) {
          throw gromos::Exception("sasa_hasel",
                  "if you give EVERY or SPEC you have to use "
                  "@timepts as well");
        }
        if (timepts.size() != 1 && timespec == "EVERY") {
          throw gromos::Exception("sasa_hasel",
                  "if you give EVERY you have to give exactly"
                  " one number with @timepts");
        }
      }
    }

    // check if we want to print the time-series
    bool sasa_ts = false;
    if (args.count("timeseries") != -1)
      sasa_ts = true;

    // check if we want to print the sasa of every atom
    bool sasa_at = false;
    if (args.count("atomic") != -1 ) {
      sasa_at = true;
      if (sasa_ts){
          throw gromos::Exception("sasa_hasel",
                  "printing of atomic SASA currently not implemented with time-series");
        }
    }

    // store sasa specifications according to IAC
    if (args.count("sasa_spec") != 1)
      throw gromos::Exception("sasa_hasel", "No sasa specification file");

    map<int, sasa_parameter> sasa_spec;
    {
      Ginstream spec_file(args["sasa_spec"]);
      vector<string> buffer;
      spec_file.getblock(buffer);
      if (buffer[0] != "SASASPEC")
        throw gromos::Exception("sasa_hasel",
              "sasa specification file does not contain a SASASPEC block!");
      if (buffer[buffer.size() - 1].find("END") != 0)
        throw gromos::Exception("sasa_hasel", "sasa specification file " + spec_file.name() +
              " is corrupted. No END in SASASPEC"
              " block. Got\n"
              + buffer[buffer.size() - 1]);

      vector<string>::const_iterator it = buffer.begin() + 1, to = buffer.end() - 1;

      for (; it != to; ++it) {
        sasa_parameter param;
        istringstream line(*it);
        unsigned int num;
        line >> param.radius >> param.probability >> param.sigma >> num;
        if (line.fail())
          throw gromos::Exception("sasa_spec",
                "bad line in SASASPEC block!");

        for (unsigned int i = 0; i < num; ++i) {
          int iac;
          line >> iac;
          if (line.fail()) {
            ostringstream msg;
            msg << "bad line in SASASPEC block: could not read " << num
                    << " IACs from line.";
            throw gromos::Exception("sasa_spec", msg.str());
          }
          sasa_spec[iac - 1] = param;
        } // for iacs
      } // SASASPEC block
    }

    // get the radius of water (if supplied)
    double R_h2o = 0.14;
    if (args.count("radius") == 1) {
      istringstream is(args["radius"]);
      is >> R_h2o;
    }

    // define input coordinate
    InG96 ic;

    // print titles for time-series (if wanted)
    if (sasa_ts) {
      cout << "# Time                  hydrophobic sasa    hydrophilic sasa   total sasa" << endl;
    }

    // declare some variables for averaging
    double ave_phobic_sasa = 0.0;
    double ave_philic_sasa = 0.0;
    double ave_tot_sasa = 0.0;

    // start at -1 to get times right
    int num_frames = -1;
    // number of time-points for which SASAs have been written
    unsigned int times_written = 0;
    // for SPEC: so that we stop trying when all requested timepoints are written
    bool done = false;

    // loop over all trajectories
    for (Arguments::const_iterator
      iter = args.lower_bound("traj"),
            to = args.upper_bound("traj");
            iter != to; ++iter) {

      // open file
      ic.open((iter->second).c_str());
      ic.select("SOLUTE");

      // loop over single trajectory
      while (!ic.eof()) {

        // declare variables for this frame
        double phobic_sasa = 0.0;
        double philic_sasa = 0.0;
        double tot_sasa = 0.0;

        ic >> sys >> time;
        // gather
        (*pbc.*gathmethod)();

        // check whether to skip or not
        num_frames++;
        if (compute_sasa(num_frames, timespec, timepts, times_written, done)) {

          // write headers for atomic sasa
          if (sasa_at) {
            cout.precision(10);
            cout.setf(ios::right, ios::adjustfield);
            cout << setw(6) << "# Atom" << setw(20) << "SASA" << endl;
          }

          // loop over molecules
          for (int m = 0; m < sys.numMolecules(); ++m) {

            // loop over atoms
            for (int i = 0; i < sys.mol(m).numAtoms(); ++i) {

              // get the bonded neighbours of atom i (pij = 0.8875)
              Neighbours neighbours(sys, m, i);

              // get radius, probability and sigma for atom i
              map<int, sasa_parameter>::const_iterator para_i =
                      sasa_spec.find(sys.mol(m).topology().atom(i).iac());

              const sasa_parameter & s = para_i->second;
              double R_i = s.radius;
              double p_i = s.probability;
              double sigma_i = s.sigma;

              // total surface area of atom i
              double S_i = 4 * pi * (R_i + R_h2o) * (R_i + R_h2o);

              // initialise multiplicative factor
              double factor = 1.0;

              // loop over (all other) atoms
              for (int j = 0; j < sys.mol(m).numAtoms(); ++j) {

                // check we are not looking at the same atom
                if (i != j) {

                  // get radius (and probability and sigma) for atom j
                  map<int, sasa_parameter>::const_iterator para_j =
                          sasa_spec.find(sys.mol(m).topology().atom(j).iac());

                  const sasa_parameter & t = para_j->second;
                  double R_j = t.radius;
                  //double p_j = t.probability;
                  //double sigma_j = t.sigma;

                  // initialise p_ij for non-neighbour case
                  double p_ij = 0.3516;
                  // check whether it is a nearest neighbour
                  Neighbours::const_iterator itn = neighbours.begin(), ton = neighbours.end();
                  for (; itn != ton; ++itn) {
                    if ( *itn == j ) {
                      // if so, pij = 0.8875
                      p_ij = 0.8875;
                      break;
                    }
                  }

                  // compute distance and reduction factor b_ij
                  double b_ij = 0.0;
                  double r_ij = (sys.mol(m).pos(i) - sys.mol(m).pos(j)).abs();
                  double check = R_i + R_j + 2 * R_h2o;
                  if (r_ij < check) {
                    // add "if not zero"?? seems to work OK...
                    b_ij = pi * (R_i + R_h2o)*(check - r_ij) *
                            (1.0 + (R_j - R_i) / r_ij);
                  }
                  // multiply factor by term for atom j
                  factor = factor * (1 - p_i * p_ij * (b_ij / S_i));
                } // i = j
              } //j loop

              // compute SASA for atom i
              double sasa_i = S_i * factor;

              // if atomic, write out SASA for atom i
              if (sasa_at) {
                cout.precision(10);
                cout.setf(ios::right, ios::adjustfield);
                cout << setw(6) << i << setw(20) << sasa_i << endl;
              }

              // add to totals for this timestep (assuming +ve sigma for hydrophobic,
              // -ve sigma for hydrophilic)
              if (sigma_i > 0.0) phobic_sasa += sasa_i;
              else if (sigma_i < 0.0 ) philic_sasa += sasa_i;
              tot_sasa += sasa_i;

            } // atoms i
          } // molecules

          // if time-series, print areas for this frame
          if (sasa_ts) {
            cout.precision(10);
            cout.setf(ios::right, ios::adjustfield);
            cout << setw(6) << time
                    << setw(20) << phobic_sasa
                    << setw(20) << philic_sasa
                    << setw(20) << tot_sasa
                    << endl;
          }

          //store values for averaging averages
          ave_phobic_sasa += phobic_sasa;
          ave_philic_sasa += philic_sasa;
          ave_tot_sasa += tot_sasa;

        }//end if compute_sasa
        if (done)
          break;

      }//end frame
    }// end io

    // print out averages (remember num_frames starts from -1)
    if ( num_frames > 0 ) {
      cout.precision(10);
      cout.setf(ios::right, ios::adjustfield);
      cout << endl << "# Averages:        hydrophobic sasa    hydrophilic sasa     total sasa" << endl;
      cout << "#         " << setw(20) << ave_phobic_sasa / times_written
              << setw(20) << ave_philic_sasa / times_written
              << setw(20) << ave_tot_sasa / times_written
              << endl;
    } else {
      cout.precision(10);
      cout.setf(ios::right, ios::adjustfield);
      cout << endl << "# Totals:          hydrophobic sasa    hydrophilic sasa     total sasa" << endl;
      cout << "#         " << setw(20) << ave_phobic_sasa / times_written
              << setw(20) << ave_philic_sasa / times_written
              << setw(20) << ave_tot_sasa / times_written
              << endl;
    }

  } catch (const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}

bool compute_sasa(int i, std::string const & timespec, vector<int> const & timepts,
        unsigned int & timesWritten, bool & done) {
  if (timespec == "ALL") {
    ++timesWritten;
    return true;
  } else if (timespec == "EVERY" && i % timepts[0] == 0) {
    ++timesWritten;
    return true;
  } else if (timespec == "SPEC") {
    for (unsigned int j = 0; j < timepts.size(); ++j) {
      if (timepts[j] == i) {
        ++timesWritten;
        if (timesWritten == timepts.size())
          done = true;
        return true;
      } // compute
    } // times
  }
  return false;
}

