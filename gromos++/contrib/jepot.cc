/**
 * @file jepot.cc
 * compute the J-value local elevation potential
 */

/**
 * @page programs Program Documentation
 *
 * @anchor jepot
 * @section jepot compute the J-value local elevation potential
 * @author @ref jallison, mc
 * @date 01. 04. 09
 *
 * compute the J-value local elevation potential
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@jval</td><td>&lt;jvalue restraint specifications&gt; </td></tr>
 * <tr><td> \@K</td><td>&lt;force constant&gt; </td></tr>
 * <tr><td> \@ngrid</td><td>&lt;number of grid points&gt; </td></tr>
 * <tr><td> \@fin</td><td>&lt;file containing final coordinates  (if not time-series)&gt; </td></tr>
 * <tr><td> [\@time</td><td>&lt;time and dt (optional and only if you want a time-series of the LE potential)&gt;] </td></tr>
 * <tr><td> [\@traj</td><td>&lt;restraint trajectory file(s) (optional: only if you want a time-series of the LE potential)&gt;] </td></tr
 * </table>
 *
 * Example:
 * @verbatim
  jepot
    @jval   ex.jval
    @K      0.01
    @ngrid  36
    @time   0 0.01
    @traj   ex.trs

    @endverbatim
 *
 * <hr>
 */

#include <cassert>
#include <sstream>

#include "../src/args/Arguments.h"
#include "../src/gio/Ginstream.h"
#include "../src/gmath/Physics.h"
#include "../src/gmath/Vec.h"
#include "../src/utils/RestrTraj.h"
#include "../src/utils/Time.h"
#include <vector>
#include <iomanip>
#include <math.h>
#include <iostream>

using namespace std;
using namespace gio;
using namespace args;
using namespace utils;

class karplus {
public:
  int m_i;
  int m_j;
  int m_k;
  int m_l;
  double weight;
  double j0;
  double delta;
  double A;
  double B;
  double C;

  karplus(int i, int j, int k, int l) {
    m_i = i;
    m_j = j;
    m_k = k;
    m_l = l;
  }

  karplus(const karplus &k) :
  m_i(k.m_i), m_j(k.m_j), m_k(k.m_k), m_l(k.m_l), weight(k.weight),
  j0(k.j0), delta(k.delta), A(k.A), B(k.B), C(k.C) {
  }
};

int main(int argc, char **argv) {

  Argument_List knowns;
  knowns << "jval" << "K" << "ngrid" << "fin" << "time" << "traj";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@jval      <jvalue restraint specifications>\n";
  usage += "\t@K           <force constant>\n";
  usage += "\t@ngrid       <number of grid points>\n";
  usage += "\t@fin         <file containing final coordinates (if not time-series)>\n";
  usage += "\t[@time       <t> <dt>] (optional and only if you want a time-series of the LE potential)\n";
  usage += "\t[@traj       <restraint trajectory files>] (optional: if you want a time-series of the LE potential)\n";

  try {
    Arguments args(argc, argv, knowns, usage);

    // get K and n
    double K = 0.0;
    int ngrid = 0;
    {
      istringstream is(args["K"]);
      is >> K;
    }
    {
      istringstream is(args["ngrid"]);
      is >> ngrid;
    }
    // set bin width
    const double w = 2.0 * M_PI / ngrid;

    // get time
    Time time(args);

    // read in the j-value specifications
    Ginstream jf(args["jval"]);
    vector<string> buffer;
    jf.getblock(buffer);

    if (buffer[0] != "JVALRESSPEC")
      throw gromos::Exception("main", "jval file does not contain an JVALRESSPEC block!");
    if (buffer[buffer.size() - 1].find("END") != 0)
      throw gromos::Exception("jval", "J-value file " + jf.name() +
            " is corrupted. No END in " + buffer[0] +
            " block. Got\n"
            + buffer[buffer.size() - 1]);

    // store atom numbers, jvalues and karplus relation components
    vector<karplus> kps;

    for (unsigned int jj = 1; jj < buffer.size() - 1; jj++) {
      istringstream is(buffer[jj]);
      int i, j, k, l;
      double fdum;
      is >> i >> j >> k >> l >> fdum;
      karplus kp(i - 1, j - 1, k - 1, l - 1);

      is >> kp.weight >> kp.j0 >> kp.delta >> kp.A >> kp.B >> kp.C;

      if (is.fail())
        throw gromos::Exception("jval", "Bad line in jval-file\n" + buffer[jj]);
      kps.push_back(kp);
    }
    jf.close();

    // check if we have (a) trajectory(ies) or just the final configuration
    bool je_ts = false;
    if (args.count("traj") >= 0)
      je_ts = true;
    
    // so we read the resetraint trajectories and compute the LE potential at each point in time

    // get the names of the trajectory file(s)
    if (je_ts == true) {
      RestrTraj je;

      // loop over all trajectories
      for (Arguments::const_iterator
        iter = args.lower_bound("traj"),
              to = args.upper_bound("traj");
              iter != to; ++iter) {

        // open file
        je.open((iter->second).c_str());
        // read in file
        while (!je.eof()) {
          je.read();
          // store jvalue eps data (for all restraints and time-points)
          JValueRestrData epsdata;
          je >> epsdata >> time;
          // loop through the eps data
          for (unsigned int i = 0; i < epsdata.data().size(); ++i ) {
              // write time and atom numbers
              cout << "#time:\t" << time << "\t\tatoms:\t" << epsdata.data()[i].i <<
              "\t" << epsdata.data()[i].j << "\t" << epsdata.data()[i].k <<
              "\t" << epsdata.data()[i].l << "\n\n";

              // loop through possible phi values in 1deg increments
              for (double phi = 0.0; phi < 2 * M_PI; phi += M_PI / 180.0) {
                // compute contribution from each bin to this phi
                double V = 0.0;
                for (int bin = 0; bin < ngrid; ++bin) {
                  // phi0 = midpoint of this bin
                  const double phi0 = (bin + 0.5) * w;
                  // correct for periodicity
                  double phi_bin = phi;
                  while (phi_bin < phi0 - M_PI)
                    phi_bin += 2 * M_PI;
                  while (phi_bin > phi0 + M_PI)
                    phi_bin -= 2 * M_PI;
                  // distance from corrected phi to midpoint of this bin
                  const double delta_phi = phi_bin - phi0;
                  // compute contribution from this phi to this bin
                  V += K * epsdata.data()[i].epsilon[bin] * exp(-delta_phi * delta_phi / (2 * w * w));
                } // end loop over bins
                // print contribution from all bins to the potential for this phi
                cout << setw(18) << 180.0 * phi / M_PI << setw(18) << V << "\n";
              } // end loop over phi
              cout << "\n\n";
            } // end restraint loop
          } // end time loop
        } // end loop over trajectories

    // read the epsilons in the case of only computing the final LE potential
    // we don't use the RestrTraj read function here because the format is different
    } else {
      Ginstream jpot(args["fin"]);

      jpot.getblock(buffer);
      while (buffer[0] != "JVALUERESEPS")
        jpot.getblock(buffer);

      if (buffer[buffer.size() - 1].find("END") != 0) {
        throw gromos::Exception("jepot", "Final coordinate file " + jpot.name() +
                " is corrupted. No END in " + buffer[0] +
                " block. Got\n"
                + buffer[buffer.size() - 1]);
      }

      // loop over jvalues
      for (unsigned int i = 1; i < buffer.size() - 1; i += 1) {
        // get epsilon for each bin
        istringstream is(buffer[i]);
        vector<double> epsilon(ngrid);
        for (int j = 0; j < ngrid; ++j)
          is >> epsilon[j];

        // write restraint number and atom numbers (read from jval file)
        cout << "\n\n# " << (i + 1) / 2 << " atoms: " << kps[i-1].m_i + 1 << " "
        << kps[i-1].m_j + 1 << " " << kps[i-1].m_k + 1 << " " << kps[i-1].m_l + 1 << " " << "\n";

        // loop through possible phi values in 1deg increments
        for (double phi = 0.0; phi < 2 * M_PI; phi += M_PI / 180.0) {
          // compute contribution from each bin to this phi
          double V = 0.0;
          for (int bin = 0; bin < ngrid; ++bin) {
            // phi0 = midpoint of this bin
            const double phi0 = (bin + 0.5) * w;
            // correct for periodicity
            double phi_bin = phi;
            while (phi_bin < phi0 - M_PI)
              phi_bin += 2 * M_PI;
            while (phi_bin > phi0 + M_PI)
              phi_bin -= 2 * M_PI;
            // distance from corrected phi to midpoint of this bin
            const double delta_phi = phi_bin - phi0;
            // compute contribution from this phi to this bin
            V += K * epsilon[bin] * exp(-delta_phi * delta_phi / (2 * w * w));
          } // end loop over bins
          // print contribution from all bins to the potential for this phi
          cout << setw(18) << 180.0 * phi / M_PI << setw(18) << V << "\n";
        } // end loop over phi
      } // end loop over j-value restraints
    } // end of if not timeseries

  }  catch (const gromos::Exception &e) {
    cerr << "EXCEPTION:\t";
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}



