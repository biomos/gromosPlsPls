/**
 * @file sasa.cc
 * Calculates solvent-accessible surface areas for selected atoms
 */

/**
 * @page programs Program Documentation
 *
 * @anchor sasa
 * @section sasa Calculates solvent-accessible surface areas for selected atoms
 * @author @ref mk
 * @date 21-6-07
 *
 * Program sasa calculates and prints the solvent-accessible surface
 * area (sasa) of all heavy atoms in the solute part of the molecular system.
 * It also calculates the contribution made by a specified set of heavy atoms.
 * The program uses the algorithm of Lee and Richards [J. Mol. Biol., 55, 379-400 (1971)].
 * A spherical probe of given radius is rolled over the surface of the molecule
 * (the size of the probe is typically 0.14~nm for water). The path traced out
 * by its centre gives the accessible surface.
 *
 * In GROMOS, the radii of the heavy atoms are obtained by calculating
 * the minimum energy distance of the interaction between the heavy
 * atom and the first solvent atom. This value is reduced by the specified
 * probe radius to account for the radius of the solvent atom.
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;boundary type&gt; [&lt;gather method&gt;] </td></tr>
 * <tr><td> [\@time</td><td>&lt;@ref utils::Time "time and dt"&gt;]</td></tr>
 * <tr><td> \@atoms</td><td>&lt;@ref AtomSpecifier "atoms" to consider for sasa&gt; </td></tr>
 * <tr><td> [\@zslice</td><td>&lt;distance between the Z-slices through the molecule (default: 0.005~nm)&gt;] </td></tr>
 * <tr><td> \@probe</td><td>&lt;probe IAC and radius&gt; </td></tr>
 * <tr><td> [\@verbose</td><td>(print summaries)] </td></tr>
 * <tr><td> \@traj</td><td>&lt;trajectory file(s)&gt; </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  sasa
    @topo     ex.top
    @pbc      r
    @time     0 1
    @atoms    1:CB
    @zslice   0.005
    @probe    5 0.14
    @verbose
    @traj     ex.tr
 @endverbatim
 *
 * <hr>
 */


#include <cassert>
#include <vector>
#include <iomanip>
#include <iostream>
#include <sstream>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/gio/InG96.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/LJException.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gio/InTopology.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/bound/Boundary.h"
#include "../src/gmath/Vec.h"
#include "../src/gmath/Physics.h"
#include "../src/fit/PositionUtils.h"
#include "../src/utils/groTime.h"
#include "../src/utils/AtomicRadii.h"

using namespace std;
using namespace gcore;
using namespace gmath;
using namespace gio;
using namespace bound;
using namespace args;
using namespace fit;
using namespace utils;

void heapsort(double* values, int n, int* key);

//some constants
double const PI = gmath::physConst.get_pi();
double const twoPI = 2 * PI;

int main(int argc, char **argv) {

  Argument_List knowns;
  knowns << "topo" << "pbc" << "time" << "zslice" << "atoms" << "probe" << "traj"
          << "verbose";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo     <molecular topology file>\n";
  usage += "\t@pbc      <boundary type>\n";
  usage += "\t[@time     <time and dt>]\n";
  usage += "\t@atoms    <atoms to consider for sasa>\n";
  usage += "\t[@zslice  <distance between the Z-slices (default: 0.005)>]\n";
  usage += "\t@probe   <probe IAC and radius>\n";
  usage += "\t[@verbose (print summaries)\n";
  usage += "\t@traj     <trajectory files>\n";

  try {
    Arguments args(argc, argv, knowns, usage);

    //   get simulation time
    Time time(args);
    int numFrames = 0;

    //try for zslice, else set default
    double zslice = args.getValue<double>("zslice", false, 0.005);
    //get probe IAC and radius
    int probe_iac;
    double probe;
    if (args.count("probe") > 1) {
      vector<double> probearg = args.getValues<double>("probe", 2);
      probe_iac = int(probearg[0]);
      probe = probearg[1];
    } else {
      throw gromos::Exception("sasa", "You need to specify the probe IAC and radius!");
    }

    //  read topology
    InTopology it(args["topo"]);
    System sys(it.system());

    utils::compute_atomic_radii_vdw(probe_iac, probe, sys, it.forceField());

    //get sasa atoms
    AtomSpecifier sasaatoms(sys);
    {
      Arguments::const_iterator iter = args.lower_bound("atoms");
      Arguments::const_iterator to = args.upper_bound("atoms");
      for (; iter != to; iter++) {
        string spec = iter->second.c_str();
        sasaatoms.addSpecifier(spec);
      }
    }

    System refSys(it.system());

    // parse boundary conditions
    Boundary *pbc = BoundaryParser::boundary(sys, args);
    // parse gather method
    Boundary::MemPtr gathmethod = args::GatherParser::parse(sys, refSys, args);

    //get radii and other things
    AtomSpecifier heavyatoms(sys);
    vector<double> radheavy;
    vector<double> radheavysq;
    int count = 0;
    double rmax = 0;
    int natoms = 0;
    //get all heavy atoms...
    for (int i = 0; i < sys.numMolecules(); ++i) {
      natoms += sys.mol(i).numAtoms();
      for (int j = 0; j < sys.mol(i).numAtoms(); ++j) {
        if (!sys.mol(i).topology().atom(j).isH()) {
          ++count;
          heavyatoms.addAtom(i, j);
          double rad = sys.mol(i).topology().atom(j).radius();
          rad += probe;
          rmax = ((rad > rmax) ? (rad) : (rmax));
          radheavy.push_back(rad);
          radheavysq.push_back(rad * rad);
        }
      }
    }
    // define input coordinate
    InG96 ic;
    
    // declare some variables
    vector<int> itab; // stores the number of atoms per cube index
    vector<int> inov; // stores neighbors
    vector<int> empty; // dummy vector for natm
    vector<double> dx; // stores x-coordinate of distances
    vector<double> dy; // stores y-coordinate of distances
    vector<double> dsq; // stores square distances
    vector<double> d; // stores distances
    vector<int> cube(natoms); // stores the cube index for each atom
    vector<vector<int> > natm(natoms);
    vector<double> accs(natoms, 0.0); // stores the accessibility

    // print title
    cout << "#     "
            << setw(15) << "selected"
            << setw(10) << "heavy" << endl
            << "# time"
            << setw(15) << "atoms"
            << setw(10) << "atoms" << endl;



    // loop over all trajectories
    for (Arguments::const_iterator
      iter = args.lower_bound("traj"),
            to = args.upper_bound("traj");
            iter != to; ++iter) {

      // open file
      ic.open((iter->second).c_str());

      // loop over single trajectory
      while (!ic.eof()) {
        ic >> sys >> time;
        (*pbc.*gathmethod)();

        double totSASA = 0;
        double totSASA_all = 0;

        //determine max and min positions using the PositionsUtils
        //using a one-dimensional array for this is not so nice...
        //well, just remember that, MAX position is stored from 0 to 2 (starting at index 0)
        //                          MIN position is stored from 3 to 5

        //this routine is fucked...
        Vec min = fit::PositionUtils::getmincoordinates(&sys, true);
        Vec max = fit::PositionUtils::getmaxcoordinates(&sys, true);

        //set up cubicals containing the atoms
        //this is analogous to the hanging spanish moss algorithm from the 70's saturday night live

        double idim = rint((max[0] - min[0]) / rmax + 0.1);
        if (idim < 3) idim = 3;
        double jidim = rint((max[1] - min[1]) / rmax + 0.1);
        if (jidim < 3) jidim = 3;
        jidim = idim * jidim;
        double kjidim = rint((max[2] - min[2]) / rmax + 0.1);
        if (kjidim < 3) kjidim = 3;
        kjidim = jidim * kjidim;

        // set the vectors to zero and resize
        itab.clear();
        inov.clear();
        dx.clear();
        dy.clear();
        dsq.clear();
        d.clear();
        empty.resize((int)kjidim, 0);
        for (int v = 0; v < natoms; v++) {
          natm[v] = empty;
          cube[v] = 0;
        }
        itab.resize((int)kjidim, 0);
        inov.resize((int)kjidim, 0);
        dx.resize((int)kjidim, 0.0);
        dy.resize((int)kjidim, 0.0);
        dsq.resize((int)kjidim, 0.0);
        d.resize((int)kjidim, 0.0);

        for (int l = 0; l < (int) heavyatoms.size(); ++l) {
          Vec tmp = *heavyatoms.coord(l);
          int i = (int) ((tmp[0] - min[0]) / rmax + 1.0);
          int j = (int) ((tmp[1] - min[1]) / rmax);
          int k = (int) ((tmp[2] - min[2]) / rmax);

          // this is the cube index
          int kji = ((int) k) * ((int) jidim) + ((int) j) * ((int) idim) + ((int) i);
          itab[kji]++;
          natm[itab[kji]][kji] = l;
          cube[l] = kji;
        }

        double nzp = rint(0.1 / (zslice) + 0.05);
        
        // loop over atoms
        for (int ir = 0; ir < count; ++ir) {

          int io = 0;
          //first loop over cubes to get the neighbors
          for (int k = -1; k <= 1; ++k) {
            for (int j = -1; j <= 1; ++j) {
              for (int i = -1; i <= 1; ++i) {
                int mkji = cube[ir] + k * ((int) jidim) + j * ((int) idim) + i;
                if (mkji >= 1) {
                  if (mkji > kjidim) {
                    i = k = j = 2; // we break the whole loop
                  } else {
                    int nm = itab[mkji];
                    if (nm >= 1 && nm < kjidim) {
                      //     -- record the atoms in inov that neighbor atom ir
                      for (int m = 1; m <= nm; ++m) {
                        int in = natm[m][mkji];
                        if (in != ir) {
                          io++;
                          if (io > kjidim) {
                            ostringstream os;
                            os << "Problem: io > kjidim";
                            throw (gromos::Exception("sasa", os.str()));
                          }
                          Vec tmp = *heavyatoms.coord(in);
                          dx[io] = heavyatoms.pos(ir)[0] - tmp[0];
                          dy[io] = heavyatoms.pos(ir)[1] - tmp[1];
                          dsq[io] = dx[io] * dx[io] + dy[io] * dy[io];
                          d[io] = sqrt(dsq[io]);
                          inov[io] = in;
                        }
                      } // nm loop
                    } // if nm >= 1
                  } // if mkji > kjidim
                } // if mkji >= 1
              } // i loop
            } // j loop
          } // k loop 

          // set some variables
          double area = 0.0; // sums the area
          
          double rr = radheavy[ir];
          double rrx2 = rr * 2;
          double rrsq = radheavysq[ir];
          
          if (io >= 1) { // we have some neighbors
            // z resolution determined
            double zres = rrx2 / nzp;
            Vec atmvec = *heavyatoms.coord(ir);
            double zgrid = atmvec[2] - rr - zres / 2;
            
            // section atom spheres perpendicular to the z axis
            // main inner loop
            for (int i = 0; i < nzp; ++i) {
              bool breakmain = false;
              double arcsum = 0; // sums the length of the arc
              zgrid += zres;

              //     find the radius of the circle of intersection of 
              //     the ir sphere on the current z-plane
              double rsec2r = rrsq - (zgrid 
                     - heavyatoms.pos(ir)[2]) * (zgrid - heavyatoms.pos(ir)[2]);
              double rsecr = sqrt(rsec2r);

              // vectors to store the start and end points of the arcs
              vector<double> arcf((int)kjidim, 0.0);
              vector<double> arci((int)kjidim, 0.0);

              int karc = -1;

              // inner loop over neighbors
              for (int j = 1; j <= io; ++j) {
                
                //find radius of circle locus
                Vec tmp = *heavyatoms.coord(inov[j]);
                double rsec2n = radheavysq[inov[j]] - ((zgrid - tmp[2]) * (zgrid - tmp[2]));
                double rsecn = sqrt(rsec2n);
                double diff_rsec = rsecr - rsecn;

                // find intersections of n.circles with ir circles in section
                // do the circles intersect, or is one circle completely inside the other?
                if (rsec2n > 0.0 && d[j] < (rsecr + rsecn)) {

                  if (d[j] < abs(diff_rsec) && diff_rsec <= 0.0) { // we break the inner loop
                    j = io;
                    breakmain = true;
                  } else if (d[j] >= abs(diff_rsec)) {
                    //if the circles intersect, find the points of intersection
                    karc++;
                    if (karc >= kjidim) {
                      ostringstream os;
                      os << "Problem: karc >= kjidim";
                      throw (gromos::Exception("sasa", os.str()));
                    }

                    //     Initial and final arc endpoints are found for the ir circle intersected
                    //     by a neighboring circle contained in the same plane. The initial endpoint
                    //     of the enclosed arc is stored in arci, and the final arc in arcf
                    //     law of cosines
                    double trig_test = (dsq[j] + rsec2r - rsec2n) / (2 * d[j] * rsecr);
                    if (trig_test >= 1.0) trig_test = 0.99999;
                    if (trig_test <= -1.0) trig_test = -0.99999;
                    double alpha = acos(trig_test);

                    //     alpha is the angle between a line containing a point of intersection and
                    //     the reference circle center and the line containing both circle centers
                    double beta = atan2(dy[j], dx[j]) + PI;
                    //     beta is the angle between the line containing both circle centers and the x-axis
                    double ti = beta - alpha;
                    double tf = beta + alpha;
                    if (ti < 0.0) ti += twoPI;
                    if (tf > twoPI) tf -= twoPI;
                    arci[karc] = ti;

                    if (tf < ti) {
                      //if the arc crosses zero, then it is broken into two segments.
                      //the first ends at twoPI and the second begins at zero
                      arcf[karc] = twoPI;
                      karc++;
                    }
                    arcf[karc] = tf;

                  } // d[j] >= abs(b)
                } // rsec2n > 0.0 && d[j] < (rsecr + rsecn)
              } // j loop

              //find the accessible surface area for the sphere ir on this section
              // only do this if (d[j] >= abs(b) && b > 0.0)
              if (!breakmain) {

                //  sum contributions
                if (karc != -1) {
                  vector<int> tag((int)kjidim, 0);
                  
                  //The arc endpoints are sorted on the value of the initial arc endpoint
                  heapsort(&arci[0], karc + 1, &tag[0]);
                  
                  arcsum = arci[0];
                  double t = arcf[tag[0]];

                  if (karc != 0) {
                    for (int k = 1; k <= karc; ++k) {
                      if (t < arci[k]) arcsum += arci[k] - t;
                      double tt = arcf[tag[k]];
                      if (tt > t) t = tt;
                    }
                  }
                  //calculate the accessible area
                  //The area/radius is equal to the accessible arc length x the section thickness.
                  arcsum += twoPI - t;
                  
                } else { // no overlap with other circles
                  arcsum = twoPI;
                }
                double parea = arcsum*zres;

                //Add the accessible area for this atom in this section to the area for this
                //atom for all the section encountered thus far
                area += parea;
                
              } // breakmain = false
            } // i loop
            
          } else { // we don't have neighbors, so calculate the exact area of the sphere
            area = twoPI * rrx2;
          }

          //scale area to vdw shell
          double scaled_area = area*rr;

          // add it for averaging
          accs[ir] += scaled_area;
          if (sasaatoms.findAtom(heavyatoms.mol(ir), heavyatoms.atom(ir)) >= 0) totSASA += scaled_area;
          totSASA_all += scaled_area;
          
        } // loop over atoms: ir < count
        
        cout.precision(8);
        cout << time;
        cout.precision(5);
        cout << setw(15) << totSASA
                << setw(20) << totSASA_all << endl;

        numFrames++;
        
      } // loop over frames in a single trajectory
      ic.close();
      
    } // loop over trajectories
    
    // calculate and print averages
    double totSASA = 0.0;
    double totSASA_all = 0.0;

    // loop over all heavy atoms
    for (int i = 0; i < count; ++i) {
      accs[i] /= numFrames;
      totSASA_all += accs[i];
      if (sasaatoms.findAtom(heavyatoms.mol(i), heavyatoms.atom(i)) >= 0) totSASA += accs[i];
    }

    cout.precision(5);
    cout << "#\n# ave."
            << setw(15) << totSASA
            << setw(10) << totSASA_all << endl;
    if (args.count("verbose") >= 0) {

      cout << "#\n# average contribution per selected heavy atom\n";
      cout << "#\n# "
              << setw(6) << "atom"
              << setw(10) << "residue"
              << setw(5) << "name" << ' '
              << setw(10) << "SASA" << ' '
              << setw(10) << "\% selected" << ' '
              << setw(10) << "\% heavy" << endl;

      double sumaccs = 0.0, sumper = 0.0, sumpera = 0.0;

      for (int i = 0; i < sasaatoms.size(); ++i) {
        int index = heavyatoms.findAtom(sasaatoms.mol(i), sasaatoms.atom(i));
        if (index >= 0) {
          cout << "# "
                  << setw(6) << sasaatoms.toString(i)
                  << setw(5) << sasaatoms.resnum(i) + 1
                  << setw(5) << sasaatoms.resname(i)
                  << setw(5) << sasaatoms.name(i) << ' '
                  << setw(10) << accs[index] << ' '
                  << setw(10) << accs[index] / totSASA * 100.0 << ' '
                  << setw(20) << accs[index] / totSASA_all * 100.0 << ' '
                  << endl;
          sumaccs += accs[index];
          sumper += accs[index] / totSASA * 100.0;
          sumpera += accs[index] / totSASA_all * 100.0;
        }
      }
      cout << "#\n# total                 "
              << setw(10) << sumaccs << ' '
              << setw(10) << sumper << ' '
              << setw(10) << sumpera << ' ' << endl;
    }
  } catch (const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}

void heapsort(double* values, int n, int* key) {

  //this implements a heapsort, which also returns the keys...
  //adapted to arrays that start indexing from 0...
  //     initialize index into the original ordering

  for (int i = 0; i < n; ++i) key[i] = i;
  //     perform the heapsort of the input list
  //		left = values.length/2;
  //		right = values.length-1;

  int k = n / 2 + 1;
  int index = n;
  double lists;
  int keys;

  do {

    if (k > 1) {
      k = k - 1;
      lists = values[k - 1];
      keys = key[k - 1];
    } else {
      lists = values[index - 1];
      keys = key[index - 1];
      values[index - 1] = values[0];
      key[index - 1] = key[0];
      index = index - 1;
    }
    if (index <= 1) {
      values[0] = lists;
      key[0] = keys;
      return;
    }

    int i = k;
    int j = k + k;
    do {

      if (j < index) {
        if (values[j - 1] < values[j]) ++j;
      }
      if (lists < values[j - 1]) {
        values[i - 1] = values[j - 1];
        key[i - 1] = key[j - 1];
        i = j;
        j = j + j;
      } else {
        j = index + 1;
      }
    } while (j <= index);

    values[i - 1] = lists;
    key[i - 1] = keys;

  } while (n > 1);


} //end void heapsort
