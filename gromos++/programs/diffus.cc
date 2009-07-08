/**
 * @file diffus.cc
 * Calculates the diffusion constant for a set of atoms
 */

/**
 * @page programs Program Documentation
 *
 * @anchor diffus
 * @section diffus Calculates the diffusion constant for a set of atoms
 * @author @ref co
 * @date 21-6-2007
 *
 * Program diffus calculates the diffusion of the centre-of-geometry of a 
 * specified set of atoms, using the Einstein equation:
 * 
 * @f[ D = \lim_{t\to\infty}\frac{<[\vec{r_0} - \vec{r}(t)]^2>}{2 N_d t} @f]
 *
 * where @f$\vec{r_0}@f$ is the centre-of-geometry in the reference
 * configuration (if none is given the program takes the first configuration of
 * the trajectory file). @f$\vec{r}(t)@f$ is the centre-of-geometry at time t.
 * @f$N_d@f$ is the number of dimensions that are being taken into account.
 *
 * The program calculates the diffusion constant by directly applying this 
 * equation as well as by applying a linear regression of the
 * mean-square-displacement. The slope of this linear regression then gives the
 * diffusion constant. 
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;boundary type&gt; [&lt;gather method&gt;] </td></tr>
 * <tr><td> [\@time</td><td>&lt;@ref utils::Time "time and dt"&gt;] </td></tr>
 * <tr><td> \@dim</td><td>&lt;dimensions to consider&gt; </td></tr>
 * <tr><td> \@atoms</td><td>&lt;@ref utils::AtomSpecifier "atom specifier": atoms to follow&gt; </td></tr>
 * <tr><td> \@ref</td><td>&lt;reference frame (r(0))&gt; </td></tr>
 * <tr><td> \@traj</td><td>&lt;trajectory files&gt; </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  diffus
    @topo  ex.top
    @pbc   r
    [@time  0 0.1]
    @dim   x y z
    @atoms s:OW
    @ref   exref.coo
    @traj  ex.tr
 @endverbatim
 *
 * <hr>
 */
//diffus calculates diffusion
#include <vector>
#include <iomanip>
#include <cmath>
#include <iostream>
#include <fstream>
#include <cassert>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/fit/Reference.h"
#include "../src/gio/InG96.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/Solvent.h"
#include "../src/gio/InTopology.h"
#include "../src/bound/Boundary.h"
#include "../src/fit/PositionUtils.h"
#include "../src/gmath/Vec.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/utils/Time.h"
#include "../src/gmath/Stat.h"

using namespace std;
using namespace gcore;
using namespace gio;
using namespace bound;
using namespace args;
using namespace utils;
using namespace gmath;

void calcD(Stat<double> & r, double t, double & diff, double & rmsd, double & err);

int main(int argc, char **argv) {

  Argument_List knowns;
  knowns << "topo" << "pbc" << "time" << "dim" << "atoms" << "ref" << "traj";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo   <molecular topology file>\n";
  usage += "\t@pbc    <boundary type> [<gather method>]\n";
  usage += "\t[@time   <time and dt>]\n";
  usage += "\t@dim    <dimensions to consider>\n";
  usage += "\t@atoms  <AtomSpecifier: atoms to follow>\n";
  usage += "\t@ref    <reference frame (r(0))>\n";
  usage += "\t@traj   <trajectory files>\n";


  try {
    Arguments args(argc, argv, knowns, usage);

    // get the @time argument
    utils::Time time(args);

    // get the relevant dimensions
    int ndim = 3;
    int dim[3] = {0, 1, 2};
    {
      Arguments::const_iterator iter = args.lower_bound("dim");
      if (iter != args.upper_bound("dim")) {
        ndim = 0;

        string dum = iter->second.c_str();
        if (dum == "x") {
          dim[ndim] = 0;
          ndim++;
        } else if (dum == "y") {
          dim[ndim] = 1;
          ndim++;
        } else if (dum == "z") {
          dim[ndim] = 2;
          ndim++;
        } else {
          throw gromos::Exception("diffus",
                  "Wrong argument given for @dim !!\n"
                  "  You can give x, y, and/or z\n"
                  "  For example:\n"
                  "    @dim  x\n"
                  "  or:\n"
                  "    @dim x y z\n");
        }
        iter++;
      }
      if (iter != args.upper_bound("dim")) {
        string dum = iter->second.c_str();
        if (dum == "x") {
          dim[ndim] = 0;
          ndim++;
        } else if (dum == "y") {
          dim[ndim] = 1;
          ndim++;
        } else if (dum == "z") {
          dim[ndim] = 2;
          ndim++;
        } else {
          throw gromos::Exception("diffus",
                  "Wrong argument given for @dim !!\n"
                  "  You can give x, y, and/or z\n"
                  "  For example:\n"
                  "    @dim  x\n"
                  "  or:\n"
                  "    @dim x y z\n");
        }
        iter++;
      }
      if (iter != args.upper_bound("dim")) {
        string dum = iter->second.c_str();
        if (dum == "x") {
          dim[ndim] = 0;
          ndim++;
        } else if (dum == "y") {
          dim[ndim] = 1;
          ndim++;
        } else if (dum == "z") {
          dim[ndim] = 2;
          ndim++;
        } else {
          throw gromos::Exception("diffus",
                  "Wrong argument given for @dim !!\n"
                  "  You can give x, y, and/or z\n"
                  "  For example:\n"
                  "    @dim  x\n"
                  "  or:\n"
                  "    @dim x y z\n");
        }
      }
    }
    //  cout << ndim << endl;
    //  for(int i=0;i<ndim;i++) cout << dim[i] << endl;

    //  read topology
    args.check("topo", 1);
    InTopology it(args["topo"]);

    //  read reference coordinates
    System refsys(it.system());
    InG96 ic;

    try {
      ic.open(args["ref"]);
    }    
    // if it didn't work, take the first frame of the trajectory
    catch (const Arguments::Exception &) {
      Arguments::const_iterator iter = args.lower_bound("traj");
      if (iter != args.upper_bound("traj"))
        ic.open((iter->second).c_str());
    }
    ic.select("ALL");
    ic >> refsys;
    ic.close();

    // we always need the old coordinates to take the nearest image
    System oldsys(refsys);

    // and the current system
    System sys(refsys);

    // parse boundary conditions
    Boundary *pbc = BoundaryParser::boundary(sys, args);

    // set atom number
    AtomSpecifier at(sys);
    {
      Arguments::const_iterator iter = args.lower_bound("atoms");
      Arguments::const_iterator to = args.upper_bound("atoms");
      for (; iter != to; iter++)
        at.addSpecifier(iter->second.c_str());
    }

    // for ease of looping, we make three of these atomspecifiers
    // one for each system
    AtomSpecifier ref_at = at;
    ref_at.setSystem(refsys);
    AtomSpecifier old_at = at;
    old_at.setSystem(oldsys);


    // the reference system already contains coordinates, here we can check
    // if the ref_at actually has a size
    if (!ref_at.size())
      throw gromos::Exception("diffus",
            "No atoms to calculate the diffusion for!");


    // calculate the com of the reference state
    Vec com0(0.0, 0.0, 0.0);
    for (int i = 0; i < ref_at.size(); i++)
      com0 += ref_at.pos(i);
    com0 /= ref_at.size();

    // values to store results
    int frames = 1;
    Vec comx;

    ofstream ts;
    ts.open("diffusts.out");
    ofstream dp;
    dp.open("diffusdp.out");
    ts << "# Time series of the direct diffusion\n# time diffus rmsd error\n";
    dp << "# Time series of the mean square displacement\n";
    vector<double> tdp;
    vector<double> tt;

    double diff, rmsd, ee;
    

    // loop over all trajectories
    for (Arguments::const_iterator
      iter = args.lower_bound("traj"), to = args.upper_bound("traj");
            iter != to; ++iter) {

      // open file
      ic.open((iter->second).c_str());
      ic.select("ALL");

      // loop over single trajectory
      while (!ic.eof()) {
        ic >> sys >> time;
        comx = Vec(0.0, 0.0, 0.0);
        Stat<double> disp_data;

        // loop over all atoms to gather with respect to their previous position
        for (int i = 0; i < at.size(); i++) {
          /*
           *at.coord(i) =
            pbc->nearestImage(*old_at.coord(i),
           *at.coord(i),
                              sys.box());
          comx+=*at.coord(i);
           */

          at.pos(i) = pbc->nearestImage(old_at.pos(i), at.pos(i), sys.box());
          comx += at.pos(i);

        }
        comx /= at.size();


        // loop over the atoms to consider
        for (int i = 0; i < at.size(); i++) {
          // calculate difference to refsys for the relevant dimensions
          // correct for com 
          for (int k = 0; k < ndim; k++) {
            //d=sys.mol(m).pos(a)[dim[k]]-comx[dim[k]]
            //   -refsys.mol(m).pos(a)[dim[k]]+com0[dim[k]];
            const double d = (at.pos(i))[dim[k]]-(ref_at.pos(i))[dim[k]];
            disp_data.addval(d*d);
          }

          // copy the current system to oldsys
          old_at.pos(i) = at.pos(i);
          // *old_at.coord(i)=*at.coord(i);
        }
        if (time.time() != 0) {
          calcD(disp_data, time.time(), diff, rmsd, ee);
          ts << time
                  << setw(20) << diff * (0.01)
                  << setw(20) << rmsd * (0.01)
                  << setw(20) << ee * (0.01) << endl;
        }

        const double disp = disp_data.ave() * ndim;
        dp << time << setw(20) << disp << endl;
        tdp.push_back(disp);
        tt.push_back(time.time());

        frames++;
      }
      ic.close();
    }
    // calculate the diffusion
    // by doing a least square fit to the average displacement
    double sx = 0, sy = 0, sxx = 0, sxy = 0;
    int N = tdp.size();

    for (int i = 0; i < N; i++) {
      sx += tt[i];
      sy += tdp[i];
      sxx += tt[i] * tt[i];
      sxy += tt[i] * tdp[i];
    }
    double a = (sxy - sx * sy / N) / (sxx - sx * sx / N);
    double b = -(a * sx - sy) / N;

    cout << "# Diffusion is calculated from the mean square displacements:\n";
    cout << "#   D = <[r0 - r(t)]^2> / (2*ndim*t)  for t -> inf\n";
    cout << "# " << endl;
    cout << "# Direct application of this relation:\n";
    cout << "#   D = " << diff * (0.01) << " cm^2/s (rmsd: " << rmsd << ", ee: " << ee << ")\n";
    cout << "# " << endl;
    cout << "# Least square fit of the mean square displacement:" << endl;
    cout << "#   <[r-r(t)]^2> = " << b << " + " << a << " * t " << endl;
    cout << "#   D = " << a / 2 / ndim * 0.01 << " cm^2/s\n";
    cout << "# " << endl;

    dp.close();
    ts.close();
  }
 catch (const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}

void calcD(Stat<double> & r, double t, double & diff, double & rmsd, double & err) {
  const double factor = 1.0 / (2 * t);
  diff = r.ave() * factor;
  rmsd = r.rmsd() * factor;
  err = r.ee() * factor;
}
