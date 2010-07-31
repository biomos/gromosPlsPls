/**
 * @file pt_top.cc
 * Combine topologies and perturbation topologies
 */

/**
 * @page programs Program Documentation
 *
 * @anchor pt_top
 * @section pt_top Combine topologies and perturbation topologies
 * @author @ref co
 * @date 7-6-07
 *
 * Combines topologies with perturbation topologies to produce new topologies 
 * or perturbation topologies. Reads a topology and a perturbation topology to
 * produce a new (perturbation) topology. The perturbation topology can contain
 * a PERTATOMPARAM or MPERTATOM block (see volume IV). The atom numbers 
 * in the perturbation topology do not need to match the numbers in the topology
 * exactly. If the topology and perturbation topology do not match in their 
 * atom numbering, a shift can be applied using the \@firstatom option.
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@pttopo</td><td>&lt;perturbation topology&gt; </td></tr>
 * <tr><td> \@type</td><td>&lt;output format: TOPO, PERTTOPO&gt; </td></tr>
 * <tr><td> \@npt</td><td>&lt;sequence number of the perturbation to apply, default 1 (state B)&gt; </td></tr>
 * <tr><td> \@firstatom</td><td>&lt;first @ref AtomSpecifier "atom" to which the perturbation will be applied&gt; </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  pt_top
    @topo      ex.top
    @pttopo    ex.pttop
    @type      P
    @npt       1
    @firstatom 1:1
 @endverbatim
 *
 * <hr>
 */
#include <cassert>
#include <iostream>
#include <sstream>
#include <iomanip>
#include "../src/args/Arguments.h"
#include "../src/gio/InTopology.h"
#include "../src/gio/InPtTopology.h"
#include "../src/gio/OutTopology.h"
#include "../src/gio/OutPtTopology.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Angle.h"
#include "../src/gcore/Bond.h"
#include "../src/gcore/Improper.h"
#include "../src/gcore/Dihedral.h"
#include "../src/gcore/CrossDihedral.h"
#include "../src/gcore/AtomPair.h"
#include "../src/gcore/PtTopology.h"
#include "../src/utils/AtomSpecifier.h"

using namespace std;
using namespace gcore;
using namespace gio;
using namespace args;

int main(int argc, char *argv[]){

  Argument_List knowns;
  knowns << "topo" << "pttopo" << "firstatom" << "npt" << "type";
  
  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo      <molecular topology file>\n";
  usage += "\t@pttopo    <perturbation topology>\n";  
  usage += "\t@type      <output format: TOPO, PERTTOPO>\n";
  usage += "\t@npt       <sequence number of the perturbation to apply, default 1 (state B)>\n";
  usage += "\t@firstatom <first atom to which the perturbation will be applied>\n";

  try{
    Arguments args(argc, argv, knowns, usage);

    // read in topology
    InTopology it(args["topo"]);
    System sys(it.system());

    //with which atom do we start?
    int start=0;
    {
      Arguments::const_iterator iter=args.lower_bound("firstatom");
      if(iter!=args.upper_bound("firstatom")){
        utils::AtomSpecifier at(sys, iter->second.c_str());
        start = at.gromosAtom(0);
      }
    }
    
    // read in multiple-perturbation-topology-file
    InPtTopology ipt(args["pttopo"]);
    PtTopology pt(ipt.ptTopo());
    
    // which perturbation do we use 
    // default is 1 which is state B
    int iipt=1;
    {
      Arguments::const_iterator iter=args.lower_bound("npt");
      if(iter!=args.upper_bound("npt"))
        iipt=atoi(iter->second.c_str())-1;
      if(iipt < 0 || iipt>=pt.numPt()) throw gromos::Exception("pt_top", 
	  "Higher perturbation specified than available");
    }

    enum output_format_enum {of_topo, of_pttopo} output_format;
    char first_letter = args["type"][0];
    if(first_letter == 'P' || first_letter == 'p')
      output_format = of_pttopo;
    else if(first_letter == 'T' || first_letter == 't')
      output_format = of_topo;
    else {
      throw gromos::Exception("pt_top", "Unkown output format.");
    }
    
    // create a title    
    ostringstream title;
    if(output_format == of_pttopo) title << "Perturbation t";
    else title << "T";
    // note the t/T on the previous lines. Beauty is in the details -
    // Comment: LMAO. Beauty in this program!?
    title << "opology based on" << endl;
    title << args["topo"] << " and " << endl;
    title << args["pttopo"];
    title << " (perturbation " << iipt + 1;
    if(start > 0) title << "; ";
    if(start > 0) {
      title << "shifted to start at atomnumber " << start+1;
    }
    title << ")";
    
    
    if(output_format == of_pttopo) {
      // create a new perturbation topology: The strategy is a bit compilcated
      // but appears to work:
      // first we create a copy of the system and apply the perturbation to it.
      // Then we use the new and the original topology to contruct a perturbation
      // topology.
      
      // copy the system
      System new_sys = sys;
      // apply the perturbation
      pt.apply(new_sys, iipt, start);
      // get the new perturbation topology from the two systems
      PtTopology new_pt(sys, new_sys);
      // write out the new perturbation topology
      OutPtTopology out(cout);
      out.setTitle(title.str());
      out.write(new_pt, &sys);
    } else if (output_format == of_topo) {
      // The strategy is simple: just apply and write out.
      // apply the perturbation
      pt.apply(sys, iipt, start);
      // write out the new topology
      OutTopology out(cout);
      out.setTitle(title.str());
      out.write(sys, it.forceField());
    } else {
      throw gromos::Exception("pt_top", "Output format not implemented.");
    }    
    return 0;
  }
  catch(gromos::Exception e){
    cerr << e.what() << endl;
    return 1;
  }
}
