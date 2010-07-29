/**
 * @file rdf.cc
 * calculates a radial distribution function
 */

/**
 * @page programs Program Documentation
 *
 * @anchor rdf
 * @section rdf calculates a radial distribution function
 * @author @ref co
 * @date 28.7.2006
 *
 * Program rdf calculates radial distribution functions over structure files or
 * trajectories. The radial distribution function, g(r), is defined here as the
 * probability of finding a particle of type J at distance r from a central
 * particle I relative to the same probability for a homogeneous distribution
 * of particles J around I. Program rdf calculates g(r) for a number of
 * discreet distances r(k), separated by distance dr as
 *
 * @f[ g(r) = \frac{N_J(k)}{4\pi r^2 dr \rho_J} @f]
 *
 * where @f$N_J(k)@f$ is the number of particles of type J found at a distance 
 * between r(k) - 1/2 dr and r(k) + 1/2 dr and @f$\rho_J@f$ is the number 
 * density of particles J. If particles I and J are of the same type, 
 * @f$\rho_J@f$ is corrected for that. At long distances, g(r) will generally 
 * tend to 1. 
 *
 * Both atoms of type I and J can be solute atoms, solvent atoms as well as 
 * @ref VirtualAtom "virtual atoms". If more than one particle of type I is
 * specified, rdf calculates the average radial distribution function for all
 * specified atoms.
 *
 * This program is parallelised.
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;boundary type&gt; </td></tr>
 * <tr><td> \@centre</td><td>&lt;@ref AtomSpecifier "atoms" to take as centre&gt; </td></tr>
 * <tr><td> \@with</td><td>&lt;@ref AtomSpecifier "atoms" to calculate distances for&gt; </td></tr>
 * <tr><td> \@cut</td><td>&lt;maximum distance&gt; </td></tr>
 * <tr><td> \@grid</td><td>&lt;number of points&gt; </td></tr>
 * <tr><td> [\@nointra</td><td>&lt;skip all intramolecular contributions&gt;] </td></tr>
 * <tr><td> \@traj</td><td>&lt;trajectory files&gt; </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  rdf
    @topo   ex.top
    @pbc    r
    @centre 1:45
    @with   s:OW
    @cut    3.0
    @grid   100
    @traj   ex.tr
 @endverbatim
 *
 * <hr>
 */
#include <vector>
#include <string>
#include <iomanip>
#include <cmath>
#include <iostream>
#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/bound/TruncOct.h"
#include "../src/bound/Vacuum.h"
#include "../src/bound/RectBox.h"
#include "../src/bound/Boundary.h"
#include "../src/fit/Reference.h"
#include "../src/fit/PositionUtils.h"
#include "../src/gio/InG96.h"
#include "../src/gio/InTopology.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/AtomPair.h"
#include "../src/gcore/LJExcType.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gcore/Solvent.h"
#include "../src/gcore/SolventTopology.h"
#include "../src/gcore/Box.h"
#include "../src/gmath/Vec.h"
#include "../src/gmath/Distribution.h"
#include "../src/utils/AtomSpecifier.h"

using namespace std;
using namespace fit;
using namespace gcore;
using namespace gmath;
using namespace gio;
using namespace bound;
using namespace args;
using namespace utils;


int main(int argc, char **argv){

  Argument_List knowns;
  knowns << "topo" << "pbc" << "centre" << "with"
         << "cut" << "grid" << "nointra" << "traj";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo   <molecular topology file>\n";
  usage += "\t@pbc    <boundary type>\n";
  usage += "\t@centre <atoms to take as centre>\n";
  usage += "\t@with   <atoms to calculate distances for>\n";
  usage += "\t@cut    <maximum distance>\n";
  usage += "\t@grid   <number of points>\n";
  usage += "\t[@nointra   <skip intramolecular atoms>]\n";
  usage += "\t@traj   <trajectory files>\n";
  
 
try{
  Arguments args(argc, argv, knowns, usage);

  // read topology
  args.check("topo",1);
  InTopology it(args["topo"]);
  System sys(it.system());

  // set centre atoms
  AtomSpecifier centre(sys);
  {
    Arguments::const_iterator iter=args.lower_bound("centre");
    Arguments::const_iterator to=args.upper_bound("centre");
    for(;iter!=to;iter++)
      centre.addSpecifier(iter->second.c_str());
  }
  
  // set atom to consider
  AtomSpecifier with(sys);
  {
    Arguments::const_iterator iter=args.lower_bound("with");
    Arguments::const_iterator to=args.upper_bound("with");
    for(;iter!=to;iter++)
      with.addSpecifier(iter->second.c_str());
  }
  
  // read in cut-off distance
  double cut=1.0;
  if(args.count("cut")>0) cut=atof(args["cut"].c_str());
  
  // read in grid number
  int grid=100;
  if(args.count("grid")>0) grid=atoi(args["grid"].c_str());

  // Check if intramolecular rdf is included
  bool nointra = false;
  if(args.count("nointra") >=0) nointra = true;

  // parse boundary conditions
  double vol_corr=1;
  Boundary *pbc = BoundaryParser::boundary(sys, args);

  if(pbc->type()=='t') vol_corr=0.5;

  // define input coordinate
  InG96 ic;

  // set up distribution arrays
  std::vector<double> rdf(grid);
  double correct=4*acos(-1.0)*cut/double(grid);
  
  for(int i=0;i<grid;i++) rdf[i]=0;
  
  // loop over all trajectories
  int count_frame=0;
  
  for(Arguments::const_iterator 
        iter=args.lower_bound("traj"),
        to=args.upper_bound("traj");
      iter!=to; ++iter){
    
    // open file
    ic.open((iter->second).c_str());
    ic.select("ALL");
    
   
    // loop over single trajectory
    while(!ic.eof()){
      ic >> sys;

      // here we have to check whether we really got the atoms we want
      // maybe the solvent is missing.
      if(centre.size()==0 || with.size()==0) {
        string argument = centre.size() == 0 ? "centre" : "width";
        throw gromos::Exception("rdf", "No atoms specified for "+argument+" atoms!");
      }
      
      // just to make absolutely sure that we treat the centre and with
      // always in the same way, sort them
      centre.sort();
      with.sort();

      // calculate the volume
      const double vol=sys.box().K().abs()*sys.box().L().abs()*sys.box().M().abs()*vol_corr;
      if(nointra == false) {
          // loop over the centre atoms
#ifdef OMP
#pragma omp parallel for
#endif
          for (int i = 0; i < centre.size(); i++) {
            gmath::Distribution dist(0, cut, grid);

            // to know if this atom is also in the with set.
            int inwith = 0;
            if(with.findAtom(centre.mol(i), centre.atom(i))>-1) inwith = 1;

            // loop over the atoms to consider
            const Vec & centre_coord = *centre.coord(i);
            for(int j = 0; j < with.size(); j++) {
              if(!(with.mol(j) == centre.mol(i) && with.atom(j) == centre.atom(i))) {
                const Vec & tmp = pbc->nearestImage(centre_coord,
                        *with.coord(j),
                        sys.box());
                dist.add((tmp - centre_coord).abs());
              }
            }
            // now calculate the g(r) for this atom
            const double dens = (with.size() - inwith) / vol;
            for(int k = 0; k < grid; k++) {
              const double r = dist.value(k);
              const double rdf_val = double(dist[k]) / (dens * correct * r * r);
#ifdef OMP
#pragma omp critical
#endif
              {
                rdf[k] += rdf_val;
              }
            }
          }
        }
        else {
#ifdef OMP
#pragma omp parallel for
#endif
          for (int i = 0; i < centre.size(); i++) {
            gmath::Distribution dist(0, cut, grid);

            // to know if this atom is also in the with set.
            int inwith = 0;
            if(with.findAtom(centre.mol(i), centre.atom(i))>-1) inwith = 1;

            // loop over the atoms to consider
            const Vec & centre_coord = *centre.coord(i);
            for(int j = 0; j < with.size(); j++) {
              if(!(with.mol(j) == centre.mol(i))) { //Here if the molecule is
                                                    //the same, the rdf is not
                                                    //calculated
                const Vec & tmp = pbc->nearestImage(centre_coord,
                        *with.coord(j),
                        sys.box());
                dist.add((tmp - centre_coord).abs());
              }
            }
            // now calculate the g(r) for this atom
            const double dens = (with.size() - inwith) / vol;
            for(int k = 0; k < grid; k++) {
              const double r = dist.value(k);
              const double rdf_val = double(dist[k]) / (dens * correct * r * r);
#ifdef OMP
#pragma omp critical
#endif
              {
                rdf[k] += rdf_val;
              }
            }
          }
        }

      count_frame++;
    }
    ic.close();
  }

  //now correct the distribution for the number of frames and the number 
  //of centre atoms
  cout << "# number of frames considered: " << count_frame << endl;
  int divide=count_frame*centre.size();
  
  for(int i=0;i<grid;i++){
    double r=(double(i)+0.5)*cut/grid;
    cout << r << "\t" << rdf[i]/double(divide) << endl;
  }
  
  }
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}

