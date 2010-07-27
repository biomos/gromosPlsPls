
/**
 * @file oparam.cc
 * balbal oparam.cc
 */

/**
 * @page contrib Contrib Program Documentation
 *
 * @anchor oparam
 * @section oparam calculates order parameter for lipids
 * @author @ref bh
 * @date 28-04-2009
 *
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;boundary type&gt; </td></tr>
 * <tr><td> [\@time</td><td>&lt;@ref utils::Time "time and dt"&gt;] </td></tr>
 * <tr><td> \@atoms</td><td>&lt;@ref AtomSpecifier&gt; </td></tr>
 * <tr><td> \@refvec</td><td>&lt;@ref x y z&gt; </td></tr>
 * <tr><td> \@traj</td><td>&lt;trajectory files&gt; </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  stacking
    @topo             ex.top
    @pbc              r
    [@time            0 1]
    @atoms            1-72:11-24
    @refvec           0 0 1
    @traj             ex.tr
 @endverbatim
 *
 * <hr>
 */

#include <cassert>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/gio/InG96.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gio/InTopology.h"
#include "../src/bound/TruncOct.h"
#include "../src/bound/Vacuum.h"
#include "../src/bound/RectBox.h"
#include "../src/gcore/AtomPair.h"
#include "../src/gcore/LJExcType.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/utils/AtomSpecifier.h"


#include <vector>
#include <iomanip>
#include <sstream>
#include <iostream>

using namespace std;
using namespace gcore;
using namespace gio;
using namespace bound;
using namespace args;
using namespace utils;
using namespace gmath;

int main(int argc, char** argv) {
  Argument_List knowns;
  knowns << "topo" << "pbc" << "atoms" 
          << "refvec" << "time"
          << "traj";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo           <molecular topology file>\n";
  usage += "\t@pbc            <boundary type> [<gathermethod>]\n";
  usage += "\t[@time           <time and dt>]\n";
  usage += "\t@atoms          <atomspecifier>\n";
  usage += "\t@refvec        <x y z>]\n";
  usage += "\t@traj           <trajectory files>\n";
  

  try {
    Arguments args(argc, argv, knowns, usage);

    //get reference vector, normalize it
    Vec refvec(0.0, 0.0, 0.0);

    {
      Arguments::const_iterator iter = args.lower_bound("refvec");
      if(iter != args.upper_bound("refvec")) {
        refvec[0] = atof(iter->second.c_str());
        ++iter;
      }
      if(iter != args.upper_bound("refvec")) {
        refvec[1] = atof(iter->second.c_str());
        ++iter;
      }
      if(iter != args.upper_bound("refvec")) {
        refvec[2] = atof(iter->second.c_str());
      }
    }
    //normalize
    refvec = refvec.normalize();
    // read topology
    InTopology it(args["topo"]);
    //make system out of topology
    System sys(it.system());

    AtomSpecifier bilayer_atoms(sys);
    {
      Arguments::const_iterator to = args.upper_bound("atoms");
      for(Arguments::const_iterator iter = args.lower_bound("atoms"); iter != to; iter++)
        bilayer_atoms.addSpecifier(iter->second);
    }


    int moln = bilayer_atoms.mol(bilayer_atoms.size() - 1) + 1;
    int num_atperlip = (bilayer_atoms.size() + 1) / moln;

    if(moln <= 0) {
      throw gromos::Exception("oparam", "molecule number cannot be <= 0!\n");
    }
    

    vector<int> atoms;
    for(int i = 0; i < num_atperlip; i++) {
      atoms.push_back(bilayer_atoms.atom(i));
    }

    if(int (atoms.size()) == 0) {
      throw gromos::Exception("oparam", "at least one atom needs to be defined!\n");
    }
    // define vector
    vector<int> at;

    for(int i = 0; i < int (atoms.size()); ++i) {
      at.push_back(atoms[i] - 1);
      if(atoms[i] == 0) {
        throw gromos::Exception("oparam", "cannot calculate the oparam for the 1st atom!\n");
      }
      at.push_back(atoms[i]);
      at.push_back(atoms[i] + 1);
    }


    // Parse boundary conditions
    Boundary *pbc = BoundaryParser::boundary(sys, args);
    // parse gather method
    Boundary::MemPtr gathmethod = args::GatherParser::parse(args);
    // loop over all trajectories
    int numFrames = 0;
    Vec z(0.0, 0.0, 0.0);
    Vec y(0.0, 0.0, 0.0);
    Vec x(0.0, 0.0, 0.0);
    Vec cos(0.0, 0.0, 0.0);
    Vec scos2(0.0, 0.0, 0.0);
    vector<Vec> S;
    for(int i = 0; i< int (at.size() / 3); ++i) {
      S.push_back(z);
    }
    vector<Vec> avcos2;
    for(int i = 0; i< int (at.size() / 3); ++i) {
      avcos2.push_back(z);
    }
    // define input coordinate
    InG96 ic;
    for(Arguments::const_iterator iter = args.lower_bound("traj");
            iter != args.upper_bound("traj"); ++iter) {
      ic.open(iter->second);

      // loop over all frames
      while(!ic.eof()) {
        numFrames++;
        ic >> sys;
        (*pbc.*gathmethod)();

        // calculate the z-vector between atom i-1 and i+1, normalize
        int cc = -2;
        for(int i = 0; i < moln; ++i) {
          cc = -2;
          for(int j = 0; j< int (at.size() / 3); j++) {
            cc += 2;
            // cout << " i: " << i << " j: " << j << " " << "cc: " << cc << endl;
            // cout << sys.mol(i).pos(at[j+1+cc])[0] << endl;
            z = sys.mol(i).pos(at[j + 2 + cc]) - sys.mol(i).pos(at[j + cc]);
            z = z.normalize();
            //calculate y, normalize
            y = (sys.mol(i).pos(at[j + 1 + cc]) - sys.mol(i).pos(at[j + cc]));
            y = y - ((sys.mol(i).pos(at[j + 1 + cc])-(sys.mol(i).pos(at[j + cc]))).dot(z)) * z;
            y = y.normalize();
            //calculate x
            x = z.cross(y);
            x = x.normalize(); //is this nescessary?
            // determine the angle
            cos[0] = refvec.dot(x);
            cos[1] = refvec.dot(y);
            cos[2] = refvec.dot(z);
            scos2 = cos*cos;
            avcos2[j] += scos2;
          }
        }

      }
    }
    ic.close();

    // average the avcos2 finally
    for(int i = 0; i < int (at.size() / 3); ++i) {
      avcos2[i] /= (numFrames * moln);
    }

    //get the S components
    for(int i = 0; i< int (at.size() / 3); ++i) {
      Vec tmp = avcos2[i];
      tmp[0] = ((3 * tmp[0]) - 1) / 2;
      tmp[1] = ((3 * tmp[1]) - 1) / 2;
      tmp[2] = ((3 * tmp[2]) - 1) / 2;
      S[i] = tmp;
    }

    //determine the  SCDOP and SCDAP, SCDEP
    vector<double> SCDOP, SCDAP, SCDEB;
    for(int i = 0; i< int (at.size() / 3); ++i) {
      Vec tmp = avcos2[i];
      Vec tmp1 = S[i];
      SCDOP.push_back(-(tmp[0]+(tmp[1] - 1) / 2));
      SCDAP.push_back(tmp1[2] / 2);
      SCDEB.push_back((2 * tmp1[0] / 3) + tmp1[1] / 3);
    }

    //print out results
    cout << setw(4) << "Atom"
            << setw(8) << "SX" << setw(8) << "SY" << setw(8) << "SZ"
            << setw(8) << "SCDOP"
            << setw(8) << "SCDAP"
            << setw(8) << "SCDEB" << endl;
    cout << endl;

    int c = -1;
    for(int i = 0; i < int (at.size() / 3); ++i) {
      c += 2;
      Vec tmp = S[i];
      cout.setf(ios::floatfield, ios::fixed);
      cout.setf(ios::right, ios::adjustfield);
      cout.precision(4);
      cout << setw(4) << at[i + c] + 1
              << setw(8) << tmp[0] << setw(8) << tmp[1] << setw(8) << tmp[2]
              << setw(8) << SCDOP[i] << setw(8) << SCDAP[i] << setw(8) << SCDEB[i] << endl;
    }
  }  catch(const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}


