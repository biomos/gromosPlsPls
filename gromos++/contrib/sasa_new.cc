#include <cassert>
#include <iostream>
#include <vector>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/bound/Boundary.h"
#include "../src/gcore/System.h"
#include "../src/gio/InG96.h"
#include "../src/gio/InTopology.h"
#include "../src/gromos/Exception.h"
#include "../src/gmath/Physics.h"
#include "../src/utils/AtomicRadii.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/utils/groTime.h"

/* Current status:
The protein is sliced, the arcs are calculated for each slice, and the cross points of the arcs
are calculated. What is missing is the following:
- calculation of the angles between the cross points and the origin of the arc
- storing of the angles in map<double, double> without
- at the end, find the piece of the arc which is remaining (loop over all cross points/angles)
- sum up all remaining pieces of the arcs in a slice
- sum up all slices
- in a second phase: think what to do with inclusions (discuss at GROMOS meeting?)

*/

// here all the stuff used for this program only is defined
namespace sasa {

  // a sphere to represent an atom with radius r

  class vec2 {
  protected:
    double pos[2];
  public:
    vec2() {};
    vec2(double x, double y);
    void setPos(double x, double y);
    double get_x();
    double get_y();
    double norm2();
    vec2 operator+(vec2 v);
    vec2 operator-(vec2 v);
  };
  
  class sphere {
  protected:
    double r;
    Vec pos;
  public:
    sphere(double r_, Vec pos_);
    Vec get_pos();
    double get_radius();
  };

  class arc {
  protected:
    vec2 pos;
    double r;
    double overlapAngle;
    map<double, double> without;
    double len;
  public:
    arc(double x_, double y_, double r_, double len_, double overlapAngle_);
    double get_r();
    double get_x();
    double get_y();
    vec2 get_pos();
    double get_len();
    void set_len(double l);
    void overlap(vector<arc>::iterator a1, vector<arc>::iterator a2);
  };

}

using namespace args;
using namespace bound;
using namespace gcore;
using namespace gio;
using namespace sasa;
using namespace std;
using namespace utils;

const double pi = gmath::physConst.get_pi();

int main(int argc, char ** argv) {

  try {

    // the (known) arguments
    Argument_List knowns;
    knowns << "topo" << "atoms" << "time" << "probe" << "zslice" << "pbc" << "traj";

    string usage = "# " + string(argv[0]);
    usage += "\n\t@topo      <molecular topology file>\n";
    usage += "\t@atoms    <solut atoms to be considered to calculate the SASA>\n";
    usage += "\t[@probe    <the IAC and LJ radius of the solvent/probe; first solvent atom is taken if not specified>]\n";
    usage += "\t[@time     <time and dt>]\n";
    usage += "\t@pbc      <periodic boundary and gathering>\n";
    usage += "\t[@zslice   <distance between the Z-slices (default: 0.005)>]\n";
    usage += "\t@traj     <trajectory files>\n";

    Arguments args(argc, argv, knowns, usage);

    // read in the topology to build the system
    args.check("topo", 1);
    InTopology it(args["topo"]);
    System sys(it.system());

    // get the atom/LJ radius of the solvent
    double solviac = sys.sol(0).topology().atom(0).iac();
    double solvrad = sys.sol(0).topology().atom(0).radius();
    {
      if (args.count("probe") == 2) {
        Arguments::const_iterator iter = args.lower_bound("probe");
        stringstream ss;
        ss << iter->second;
        ++iter;
        ss << iter->second;
        ss >> solviac >> solvrad;
        if (ss.fail() || ss.bad()) {
          stringstream msg;
          msg << args["probe"] << ": wrong arguments for @probe";
          throw gromos::Exception("sasa", msg.str());
        }
      } else if (args.count("probe") > -1) {
        throw gromos::Exception("sasa", "wrong arguments for @probe");
      }
    }

    // calculate the van der Waals radii
    compute_atomic_radii_vdw(solviac, solvrad, sys, it.forceField());

    // get the list of atoms to be considered for the SASA calculation
    AtomSpecifier atoms(sys);
    {
      Arguments::const_iterator iter = args.lower_bound("atoms");
      Arguments::const_iterator to = args.upper_bound("atoms");
      for (; iter != to; iter++) {
        atoms.addSpecifier(iter->second.c_str());
      }
    }
    if (atoms.size() < 1) {
      stringstream msg;
      throw gromos::Exception("sasa", "no atoms found for given atom specifier (@atoms)");
    }

    // is there a trajectory file?
    if (args.count("traj") < 1) {
      throw gromos::Exception("sasa", "no trajectory file specified (@traj)");
    }

    // get the boundary and gather method
    Boundary *pbc = BoundaryParser::boundary(sys, args);
    // parse gather method
    Boundary::MemPtr gathmethod = args::GatherParser::parse(sys, pbc->refSys(), args);

    // get the simulation time
    Time time(args);

    // get the distance between the x/y-planes
    double dz = 0.0005;
    if (args.count("zslice") > 0) {
      stringstream ss;
      ss << args.count("zslice");
      ss >> dz;
      if (ss.fail() || ss.bad()) {
        throw gromos::Exception("sasa", "bad value for @zslice");
      }
    }

    // loop over the trajectory file
    Arguments::const_iterator iter = args.lower_bound("traj");
    Arguments::const_iterator to = args.upper_bound("traj");
    for (; iter != to; ++iter) {
      InG96 ic;
      ic.open((iter->second).c_str());
      // loop over the different configurations/frames of the trajectory file
      while (!ic.eof()) {

        // read the current frame/configuration
        ic >> sys >> time;

        // gather the current configurations
        (*pbc.*gathmethod)();

        // build a sphere for each (solute) atom to be considered
        vector<sphere> spheres;
        for (int a = 0; a < atoms.size(); ++a) {
          double r = atoms.radius(a);
          Vec pos = atoms.pos(a);
          //cerr << "r = " << r << endl;
          spheres.push_back(sphere(r, pos));
        }

        // get z_min and z_max (to span the x/y-planes)
        double z_min = spheres[0].get_pos()[2] -
                spheres[0].get_radius() - solvrad;
        double z_max = spheres[0].get_pos()[2] +
                spheres[0].get_radius() + solvrad;
        for (unsigned int s = 0; s < spheres.size(); ++s) {
          double z = spheres[s].get_pos()[2];
          if ((z - spheres[s].get_radius() - solvrad) < z_min) {
            z_min = z - spheres[s].get_radius() - solvrad;
          } else if ((z + spheres[s].get_radius() + solvrad) > z_max) {
            z_max = z + spheres[s].get_radius() + solvrad;
          }
        }

        // overlap the spheres with the N x/y-planes
        // first get the number of planes needed, N
        int N = floor((z_max - z_min) / dz);
        if (N == 0) {
          N = 1;
        }
        //cerr << "N = " << N << endl;
        // get a vector containing the circles/arcs per plane
        vector<vector<arc> > arcs(N);
        // loop over planes
        for (int n = 0; n < N; ++n) {
          // the z-position of the current plane
          double z0 = z_min + (n + 1) * dz; // we don't start at z_min since there
          // is by definition no sphere-plane overlap
          // loop over all spheres and add a circle/arc if the sphere cuts the current plane
          for (unsigned int s = 0; s < spheres.size(); ++s) {
            double R = spheres[s].get_radius() + solvrad; // sphere radius (including solvent/probe radius)
            double z_sphere = spheres[s].get_pos()[2]; // z-position of centre of sphere
            double dist = abs(z0 - z_sphere); // distance of centre of sphere to the current plane
            //cerr << "R = " << R << ", z_sphere = " << z_sphere << endl;
            //cerr << "dist = " << dist << endl;
            if (dist < R) { // the sphere cuts the current plane
              double x = spheres[s].get_pos()[0]; // x-position of circle/arc
              double y = spheres[s].get_pos()[1]; // y-position of circle/arc
              double r = R * sin(acos(dist / R)); // radius of circle/arc
              double len = 2 * pi * r; // the length of the arc, actually the whole circle
              arc a(x, y, r, len, 0.0);
              arcs[n].push_back(a);
            }
          } // end of loop over spheres

          // now overlap the circles within one plane and remove intersections
          // loop over circles/arcs within a plane
          vector<arc>::iterator it1 = arcs[n].begin();
          //vector<arc>::iterator it2 = arcs[n].begin();
          for (int a1 = 0; a1 < ((int) arcs[n].size() - 1); ++a1, ++it1) {
            vector<arc>::iterator it2 = it1;
            ++it2;
            for (unsigned int a2 = a1 + 1; a2 < arcs[n].size(); ++a2, ++it2) {
              it1->overlap(it1, it2);
            } // end of loop over (second) circles
          } // end of loop over (first) circles
          //cerr << "plane " << n << ": " << arcs[n].size() << " circles\n";
        } // end of loop over planes

      } // end of loop over configurations/frames
    } // end of loop over trajectory files
    
    //vec2 v1(1,1);
    //vec2 v2(2,3);
    
    //cerr << "v1 + v2 = " << (v1 + v2).get_x() << "/" << (v1 + v2).get_y() << endl;

  } catch (const gromos::Exception &e) {

    // quit with an error message
    cerr << e.what() << endl;
    exit(1);

  }

  return 0;
}

// -----------------------------------------------------------------------------
// ======================= Function, Class, ... Definitions ====================
//
// first all the stuff in the sasa namespace
namespace sasa {

  vec2::vec2(double x, double y) {
    pos[0] = x;
    pos[2] = y;
  }
  
  void vec2::setPos(double x, double y) {
    pos[0] = x;
    pos[1] = y;
  }
  
  double vec2::get_x() {
    return pos[0];
  }
  
  double vec2::get_y() {
    return pos[1];
  }
  
  double vec2::norm2() {
    return pos[0] * pos[0] + pos[1] * pos[1];
  }
  
  vec2 vec2::operator+(vec2 v) {
    vec2 r(v.get_x() + get_x(), v.get_y() + get_y());
    return r;
  }
  
  vec2 vec2::operator-(vec2 v) {
    v.setPos(v.get_x() - pos[0], v.get_y() - pos[1]);
    return v;
  }
  
  sphere::sphere(double r_, Vec pos_) {
    r = r_;
    pos = pos_;
  }

  Vec sphere::get_pos() {
    return pos;
  }

  double sphere::get_radius() {
    return r;
  }

  arc::arc(double x_, double y_, double r_, double len_, double overlapAngle_) {
    pos.setPos(x_, y_);
    r = r_;
    len = len_;
    overlapAngle = overlapAngle_;
  }

  double arc::get_r() {
    return r;
  }
  
  double arc::get_x() {
    return pos.get_x();
  }
  
  double arc::get_y() {
    return pos.get_y();
  }

  vec2 arc::get_pos() {
    return pos;
  }
  
  double arc::get_len() {
    return len;
  }
  
  void arc::set_len(double l) {
    len = l;
  }

  void arc::overlap(vector<arc>::iterator a1, vector<arc>::iterator a2) {
    // distance of the two centres
    double dx = (a1->get_x() - a2->get_x());
    double dy = (a1->get_y() - a2->get_y());
    double d = sqrt(dx * dx + dy * dy);
    // now treat the different cases of overlap and no overlap
    if (d >= a1->get_r() + a2->get_r()) { // no overlap, circles too far away
      return;
    } else if (abs(a1->get_r() - a2->get_r()) >= d) { // there is no overlap, one circle in the other
      if (a1->get_r() < a2->get_r()) {
        a1->set_len(0.0);
      } else {
        a2->set_len(0.0);
      }
    } else { // there is an overlap
        // overlap the two circles:
      // r1^2 = (x - x1)^2 + (y - y1)^2     (I)
      // r2^2 = (x - x2)^2 + (y - y2)^2    (II)
      // (II) - (I) => y = mx + n with    (III)
      double R1 = (a1->get_r() * a1->get_r()) - (a1->get_x() * a1->get_x()) - (a1->get_y() * a1->get_y());
      double R2 = (a2->get_r() * a2->get_r()) - (a2->get_x() * a2->get_x()) - (a2->get_y() * a2->get_y());
      double m = -dx / dy;
      double n = (R2 - R1) / (2 * dy);
      // insertion of y = mx + n in (I) leads to x_{1,2}
      // calculation of x_{1,2} using the p-q-formula:
      double m2i = 1.0 / (1 + m * m);
      double p = 2 * (m * n - a1->get_x() - a1->get_y() * m) * m2i;
      double q = (n * n - 2 * n * a1->get_y() - R1) * m2i;
      double factor2 = sqrt((p * p / 4.0) - q);
      double factor1 = -p / 2.0;
      double x1 = factor1 + factor2;
      double x2 = factor1 - factor2;
      // insert x_{1,2} into (III)
      double y1 = m * x1 + n;
      double y2 = m * x2 + n;
      //vec2 p1(x1 - a1->get_x(), y1 - a1->get_y());
      //vec2 p2(x2 - a1->get_x(), y2 - a1->get_y());
      //vec2 pp1(x1, y1);
      //vec2 pp2(x2, y2);
      // compare the norms to r^2
      cerr.precision(9);
      cerr << "point 1: " << x1 << ", " << y1 << endl;
      cerr << "point 2: " << x2 << ", " << y2 << endl;
      cerr << "center 1: " << a1->get_x() << ", " << a1->get_y() << ", radius = " << a1->get_r() << endl;
      cerr << "center 2: " << a2->get_x() << ", " << a2->get_y() << ", radius = " << a2->get_r() << endl;
      //cerr << "r * r = " << a1->get_r() * a1->get_r() << endl;
      //cerr << "p1.norm2() = " << p1.norm2() << endl;
      //cerr << "p2.norm2() = " << p2.norm2() << endl;
      //cerr << "(c - pp1).norm2() = " << (a1->get_pos() - pp1).norm2() << endl;
      //cerr << "(c - pp2).norm2() = " << (a1->get_pos() - pp2).norm2() << endl << endl;

    }
  }

}
