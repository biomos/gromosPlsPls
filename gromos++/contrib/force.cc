/**
 * @file force.cc
 * Calculates the total electrostatic force on atoms of group A from atoms of group B. 
 * 
 * @page contrib Contrib program documentation
 *
 * @anchor force
 * @section calculates the electrostatic force between two groups of atoms
 * @author @ref ae
 * @date November 16, 2012
 *
 * <hr>
 */

#include <iostream>
#include <cassert>
#include <iomanip>
#include <list>
#include "../src/args/Arguments.h"
#include "../src/gio/InTopology.h"
#include "../src/gcore/System.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/gcore/GromosForceField.h"
#include "../src/bound/Boundary.h"
#include "../src/args/GatherParser.h"
#include "../src/args/BoundaryParser.h"
#include "../src/gio/InG96.h"
#include "../src/utils/groTime.h"
#include "../src/utils/SimplePairlist.h"
#include "../src/gmath/Physics.h"
#include "../src/utils/Value.h"
#include "../src/utils/VectorSpecifier.h"
#include "gcore/LJType.h"
#include "gcore/Exclusion.h"

using namespace std;
using namespace args;
using namespace gcore;
using namespace gio;
using namespace utils;
using namespace bound;
using namespace gmath;

int main(int argc, char **argv) {

  Argument_List knowns;
  knowns << "topo" << "pbc" << "pairlist" << "cut" << "cutmin" << "epskap"
          << "atomsA" << "atomsB" << "trc" << "verbose" << "projvec";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo   <molecular topology file>\n";
  usage += "\t[@pbc     <periodic boundary conditions>]\n"
           "\t          (only useful in case you want to overrule the entries of GENBOX)\n";
  usage += "\t@pairlist <type (CHARGEGROUP or ATOMIC)>\n";
  usage += "\t@cut      <cut-off radius for the force calculations>\n";
  usage += "\t[@cutmin  <only consider pairs i,j with r_ij >= cutmin>]\n";
  usage += "\t@epskap   <epsilon and kappa to be used Coulomb/RF calculation>\n";
  usage += "\t@atomsA   <atoms of group A (atom specifier)>\n";
  usage += "\t@atomsB   <atoms of group B (atom specifier)>\n";
  usage += "\t[@projvec  <a vector specifier to project the force vector to (vector specifier)>]\n";
  usage += "\t           (e.g. atom(1:1,2) the vector pointing from atom 2 to 1,\n";
  usage += "\t                 cart(x,y,z) the vector with Cartesian coordinates x, y, and z)\n";
  usage += "\t@verbose  (prints some information about the variables used)\n";
  usage += "\t@trc      <positional simulation trajectory files>\n";

  try {
    Arguments args(argc, argv, knowns, usage);

    // read topology
    args.check("topo", 1);
    InTopology it(args["topo"]);
    System sys(it.system());
    System refSys(it.system()); // a reference topology, in our case the same as the actual topology

    // The GROMOS force field as read from the topology
    GromosForceField gff(it.forceField());

    // parse boundary conditions (from GENBOX, if no @pbc is given, from @pbc else)
    Boundary *pbc = BoundaryParser::boundary(sys, args);
    
    // parameters of the non-bonded interactions
    args.check("cut", 1);
    args.check("epskap", 2);
    double cut, eps, kappa;
    {
      stringstream ss;
      ss << args["cut"];
      ss >> cut;
      ss.clear();
      ss.str("");
      Arguments::const_iterator start = args.lower_bound("epskap");
      Arguments::const_iterator end = args.upper_bound("epskap");
      for (; start != end; start++) {
        ss << start->second << endl;
      }
      ss >> eps >> kappa;
    }
    
    double cutmin = 0.0;
    if(args.count("cutmin") >= 1) {
      stringstream ss;
      ss << args["cutmin"];
      ss >> cut;
    }

    // create Time object to read the time from trajectory
    Time time(args);

    // add the atoms of groups A and B to the atom specifiers
    args.check("atomsA", 1);
    args.check("atomsB", 1);
    AtomSpecifier atomsA(sys);
    AtomSpecifier atomsB(sys);
    {
      Arguments::const_iterator iter = args.lower_bound("atomsA");
      Arguments::const_iterator to = args.upper_bound("atomsA");
      for (; iter != to; iter++) {
        atomsA.addSpecifier(iter->second.c_str());
      }
      iter = args.lower_bound("atomsB");
      to = args.upper_bound("atomsB");
      for (; iter != to; iter++) {
        atomsB.addSpecifier(iter->second);
      }
    }
    // make sure the atoms are nicely sorted
    atomsA.sort();
    atomsB.sort();
    
    // print some information in case of @verbose
    if(args.count("verbose") >= 0) {
      cerr << "# cut = " << cut << endl;
      cerr << "# eps = " << eps << endl;
      cerr << "# kap = " << kappa << endl << "#" << endl;
      cerr << "# number of atoms in group A: " << atomsA.size() << endl;
      cerr << "# number of atoms in group B: " << atomsB.size() << endl << "#" << endl;
    }
    
    // loop over the trajectory files
    for (Arguments::const_iterator iter = args.lower_bound("trc");
            iter != args.upper_bound("trc"); ++iter) {

      // define input coordinates
      InG96 ic;
      ic.open(iter->second);
      ic.select("ALL");
      

      // initiate the pair list to be used later
      map<int, SimplePairlist> pl; // key (int) = gromos atom number
                                   //   for each atom of group A there is one key and one pair list
                                   // value (SimplePairlist): pair list from center
                                   //   atom "key" to all atoms within group B
      args.check("pairlist", 1);
      string type = args["pairlist"];
      for (int a = 0; a < atomsA.size(); a++) {
        int gromosNumA = atomsA.gromosAtom(a);
        int molNumA = atomsA.mol(a);
        int atomNumA = atomsA.atom(a);
        SimplePairlist _pl(sys, *pbc, cut);
        pl[gromosNumA] = _pl;
        pl.find(gromosNumA)->second.setType(type);
        pl.find(gromosNumA)->second.setAtom(molNumA, atomNumA); // set the center atom
      }
      
      // a vector to keep the forces
      vector<Vec> force_CRF(atomsA.size());
      vector<Vec> force_LJ(atomsA.size());
      
      // the reaction field constant
      const double crf = ((2 - 2 * eps) * (1 + kappa * cut) - eps * (kappa * kappa * cut * cut)) /
              ((1 + 2 * eps)*(1 + kappa * cut) + eps * (kappa * kappa * cut * cut));
      
      // print the header of the table
      if(args.count("projvec")>=1) {
        cout << "#" << setw(14) << "time" << setw(20) << "f_CRF_x" << setw(20) << "f_CRF_y" << setw(20) << "f_CRF_z" << setw(20)
                << "f_LJ_x" << setw(20) << "f_LJ_y" << setw(20) << "f_LJ_z" << setw(20)
                << "f_tot_x" << setw(20) << "f_tot_y" << setw(20) << "f_tot_z" << setw(20)
                << "projection_CRF" << setw(20) << "projection_LJ"  << setw(20) << "projection_tot" << endl;
      } else {
        cout << "#" << setw(14) << "time" << setw(20) << "f_CRF_x" << setw(20) << "f_CRF_y" << setw(20) << "f_CRF_z" << setw(20)
                << "f_LJ_x" << setw(20) << "f_LJ_y" << setw(20) << "f_LJ_z" << setw(20)
                << "f_tot_x" << setw(20) << "f_tot_y" << setw(20) << "f_tot_z" << endl;
      }
      
      // loop over all frames
      while (!ic.eof()) {
        
        // read the configuration and the time of the current frame
        ic >> sys >> time;
        
        // the reference vector for the projection, if specified
        Vec e(0.0, 0.0, 0.0);
        if(args.count("projvec") >= 1) {
          VectorSpecifier vs(sys, pbc, args["projvec"]);
          e = vs().normalize();
        }
        
        // calculate the pair list for all atoms of group A
#ifdef OMP
#pragma omp parallel for
#endif
        for(int a = 0; a < atomsA.size(); a++) {
          
          // set the forces on atom a (within group A) to zero at the beginning
          force_CRF[a] = Vec(0.0, 0.0, 0.0);
          force_LJ[a] = Vec(0.0 ,0.0 ,0.0);
          
          int gromosNumA = atomsA.gromosAtom(a);
          Vec posA = atomsA.pos(a);
          double chargeA = atomsA.charge(a);
          map<int, SimplePairlist>::iterator it_pl = pl.find(gromosNumA); // the pair list to the current atom a within the group A
          it_pl->second.clear();
          it_pl->second.calc(atomsB, cutmin);
          it_pl->second.removeExclusions();
          
          // and the force acting on the atoms of group A by looping over
          // the according pair list
          
          for(int b = 0; b < it_pl->second.size(); b++) {
            Vec posB = pbc->nearestImage(posA, it_pl->second.pos(b), sys.box());
            Vec r_vec = posA - posB;
            double r = r_vec.abs();
            
            // the CRF force (without RF terms from excluded atoms)
            double chargeB = it_pl->second.charge(b);
            //cerr << "chargeA = " << chargeA << endl << "chargeB = " << chargeB << endl;
            double qq = chargeA * chargeB;
            Vec f = (qq * (physConst.get_four_pi_eps_i())) * (1 / (r * r * r) + crf / (cut * cut * cut)) * r_vec;
            force_CRF[a] += f;
            //cerr << "qq = " << qq << endl;
            //cerr << "f = (" << f[0] << "," << f[1] << "," << f[2] << ")\n";
            //cerr << "force_CRF[" << a << "] = (" << force_CRF[a][0] << "," << force_CRF[a][1] << "," << force_CRF[a][2] << ")\n";
            
            // the LJ force
            LJType lj(gff.ljType(AtomPair(atomsA.iac(a), it_pl->second.iac(b))));
            // exclusions were already kicked out after building the pair list,
            // but we have to check for the special 1,4-interactions
            bool special14 = false;
            {
              // (14) exclusions are from the atom with the smaller gromos number to that one with the higher, so
              // make sure we agree with that convention ("smaller" atom excludes "bigger" atom, in terms of sequential numbers)
              int atomNum1 = atomsA.atom(a) <= it_pl->second.atom(b) ? atomsA.atom(a) : it_pl->second.atom(b);
              int atomNum2 = atomsA.atom(a) <= it_pl->second.atom(b) ? it_pl->second.atom(b) : atomsA.atom(a);
              // loop over the 14 exclusions of the "smaller" atom to see if the other one belong to the set of special 14 exclusion
              for(int e = 0; e < sys.mol(atomsA.mol(a)).topology().atom(atomNum1).exclusion14().size(); e++) {
                // if the two atoms are in a different molecule, there is nothing to do at all
                if (atomsA.mol(a) != it_pl->second.mol(b)) {
                  continue;
                }
                if(atomNum2 == sys.mol(atomsA.mol(a)).topology().atom(atomNum1).exclusion14().atom(e)) {
                  special14 = true;
                  break;
                }
              }
            }
            // get the LJ parameters, according to if it is a normal or special 14 LJ atom pair
            double c12, c6;
            if(special14) {
              c12 = lj.cs12();
              c6 = lj.cs6();
            } else {
              c12 = lj.c12();
              c6 = lj.c6();
            }
            //cerr << "c12 = " << c12 << endl << "c6 = " << c6 << endl;
            // and finally calculate the force of the LJ interaction
            double r2 = r*r;
            f = ((2*c12/(r2 * r2 *r2))-c6) * (6 * r_vec/(r2 * r2 * r2 * r2 ));
            force_LJ[a] += f;
          }
          
          // add the RF terms from the excluded atoms to the force_CRF
          // (one may improve the code's speed, but at the moment we just want a working solution)
          for(int b = 0; b < atomsB.size(); b++) {
            // there is nothing to do if A and B are the same atom
            if(atomsA.gromosAtom(a) == atomsB.gromosAtom(b)) {
              continue;
            }
            // we have to loop over the exclusions of the atom with the lower
            // GROMOS number to see if the other one is in the set of excluded atoms
            bool excluded = false;
            if(atomsA.gromosAtom(a) < atomsB.gromosAtom(b)) {
              // loop over the exclusions of atom A to see if B is in this set
              int gromosNum1 = atomsA.gromosAtom(a);
              int gromosNum2 = atomsB.gromosAtom(b);
              for(int e = 0; e < atomsA.sys()->mol(atomsA.mol(a)).topology().atom(atomsA.atom(a)).exclusion().size(); e++) {
                int excludedGromosAtomNum = gromosNum1 + atomsA.sys()->mol(atomsA.mol(a)).topology().atom(atomsA.atom(a)).exclusion().atom(e) - atomsA.atom(a);
                if(excludedGromosAtomNum == gromosNum2) {
                  excluded = true;
                  break;
                }
              }
            } else {
              // loop over the exclusions of B to see if A is in this set
              int gromosNum1 = atomsB.gromosAtom(b);
              int gromosNum2 = atomsA.gromosAtom(a);
              for(int e = 0; e < atomsB.sys()->mol(atomsB.mol(b)).topology().atom(atomsB.atom(b)).exclusion().size(); e++) {
                int excludedGromosAtomNum = gromosNum1 + atomsB.sys()->mol(atomsB.mol(b)).topology().atom(atomsB.atom(b)).exclusion().atom(e) - atomsB.atom(b);
                if(excludedGromosAtomNum == gromosNum2) {
                  excluded = true;
                  break;
                }
              }
            }
            if (excluded) {
              double chargeB = atomsB.charge(b);
              double qq = chargeA * chargeB;
              Vec posB = pbc->nearestImage(posA, atomsB.pos(b), sys.box());
              Vec r_vec = posA - posB;
              Vec f = qq * (physConst.get_four_pi_eps_i()) * crf / (cut * cut * cut) * r_vec;
              force_CRF[a] += f;
            }
          }

        }
        
        // calculate the average CRF and LJ force over all atoms of group A
        Vec f_CRF(0.0, 0.0, 0.0);
        Vec f_LJ(0.0, 0.0, 0.0);
        for(int a = 0; a < atomsA.size(); a++) {
          f_CRF += force_CRF[a];
          f_LJ += force_LJ[a];
        }
        f_CRF /= (double) atomsA.size();
        f_LJ /= (double) atomsA.size();
        
        // in case there is a projection vector given , the output prints the force vector and its projection to e, otherwise
        cout.precision(9);
        if(e.abs2() > 0) {
          cout << setw(15) << time << scientific << setw(20)
                  << f_CRF[0] << setw(20) << f_CRF[1] << setw(20) << f_CRF[2] << setw(20)
                  << f_LJ[0] << setw(20) << f_LJ[1] << setw(20) << f_LJ[2] << setw(20)
                  << f_CRF[0] + f_LJ[0] << setw(20) << f_CRF[1] + f_LJ[1] << setw(20) << f_CRF[2] + f_LJ[2] << setw(20)
                  << f_CRF.dot(e) << setw(20) << f_LJ.dot(e) << setw(20) << f_CRF.dot(e) + f_LJ.dot(e) << endl;
        } else {
          cout << setw(15) << time << scientific << setw(20)
                  << f_CRF[0] << setw(20) << f_CRF[1] << setw(20) << f_CRF[2] << setw(20)
                  << f_LJ[0] << setw(20) << f_LJ[1] << setw(20) << f_LJ[2] << setw(20)
                  << f_CRF[0] + f_LJ[0] << setw(20) << f_CRF[1] + f_LJ[1] << setw(20) << f_CRF[2] + f_LJ[2] << endl;
        }
        
      }
    }
    
  } catch (const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}

