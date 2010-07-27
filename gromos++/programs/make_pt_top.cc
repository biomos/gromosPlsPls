/**
 * @file make_pt_top.cc
 * Creates perturbation topologies from molcular topologies
 */

/**
 * @page programs Program Documentation
 *
 * @anchor make_pt_top
 * @section make_pt_top Create a perturbation topology from two molecular topologies
 * @author @ref ns
 * @date 22.01.08
 * Program make_pt_top takes two or more molecular topologies and writes the 
 * differences in the perturbation topology format. Both topologies must contain
 * the same number of solute atoms. The softness parameters @f$\alpha_{LJ}@f$ 
 * @f$\alpha_{CRF}@f$ can be specified. 
 * 
 * The resulting perturbation topology file is written out to the standard output.
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology files for states&gt; </td></tr>
 * <tr><td> \@softpar</td><td>&lt;softness parameters (@f$\alpha_{LJ}@f$ @f$\alpha_{CRF}@f$)&gt; </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  make_pt_top
    @topo     exA.top exB.top
 @endverbatim
 *
 * <hr>
 */

#include <cassert>
#include <map>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <algorithm>
#include <iterator>
#include <map>
#include <set>
#include "../src/args/Arguments.h"
#include "../src/gio/InTopology.h"
#include "../src/gcore/System.h"
#include "../src/gcore/AtomPair.h"
#include "../src/gcore/Bond.h"
#include "../src/gcore/Angle.h"
#include "../src/gcore/Improper.h"
#include "../src/gcore/Dihedral.h"
#include "../src/gcore/CrossDihedral.h"
#include "../src/gcore/PtTopology.h"
#include "../src/gio/OutPtTopology.h"

using namespace std;
using namespace gcore;
using namespace gio;
using namespace args;


int main(int argc, char *argv[]){

  Argument_List knowns; 
  knowns << "topo" << "softpar";
  
  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo      <topology files for states>\n";
  usage += "\t[@softpar   <alpha_lj alpha_crf (default 1.51 0.5)>]\n"; 
  
  try{
    Arguments args(argc, argv, knowns, usage);
    
    ostringstream title;
    title << argv[0] << " generated perturbation topology." << endl;

    // get the topologies
    vector<System*> systems;
    {
      unsigned int i = 1;
      for (Arguments::const_iterator it = args.lower_bound("topo"),
           to = args.upper_bound("topo"); it != to; ++it, ++i) {
        title << "State " << i << ": " << it->second << endl;
        InTopology intopo(it->second);
        systems.push_back(new System(intopo.system()));
      }
    }
    
    if (systems.size() < 2) {
      throw gromos::Exception("make_pt_top", "Give at least two states.");
    }
    
    double alpha_lj = 1.51, alpha_crf = 0.5;
    // get the softness parameters
    if (args.count("softpar") > 0) {
      if (args.count("softpar") != 2) {
        throw gromos::Exception("make_pt_top", "Both softness parameters "
                "(alpha_LJ and alpha_CRF) have to be given.");
      }
      Arguments::const_iterator iter=args.lower_bound("softpar");
      istringstream arg1(iter->second), arg2((++iter)->second);
      if (!(arg1 >> alpha_lj) || !(arg2 >> alpha_crf)) {
        throw gromos::Exception("make_pt_top", "Make sure the softness "
                "parameters are floating point numbers.");
      }
    }
    
    OutPtTopology out_pt(cout);
    out_pt.setTitle(title.str());
    
    // create the perturbation from the two systems
    PtTopology *pttop;
    if (systems.size() == 2) {
      pttop = new PtTopology(*systems[0], *systems[1]);
    } else { // or from more systems
      cerr << "Creating multiple perturbation topology - ignoring bonded and exclusions." << endl;
      pttop = new PtTopology(systems);
    }
    
    // set the softness
    for(int i = 0; i < pttop->numAtoms(); ++i) {
      pttop->setAlphaLJ(i, alpha_lj);
      pttop->setAlphaCRF(i, alpha_crf);
    }
    
    // write it out
    if (systems.size() == 2)
      out_pt.write(*pttop, systems[0]);
    else
      out_pt.write_multiple(*pttop);
    
    // clean up
    delete pttop;
    for(vector<System*>::const_iterator it = systems.begin(),
        to = systems.end(); it != to; ++it)
      delete *it;
    
    return 0;
  }
  catch(gromos::Exception e){
    cerr << e.what() << endl;
    return 1;
  }
}
