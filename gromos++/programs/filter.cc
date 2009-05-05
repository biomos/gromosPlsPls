/**
 * @file filter.cc
 * Filters a coordinate trajectory for a specific set of atoms
 */

/**
 * @page programs Program Documentation
 *
 * @anchor filter
 * @section filter Filters a coordinate trajectory for a specific set of atoms
 * @author @ref co
 * @date 11-6-07
 *
 * Program filter reduces coordinate trajectory file(s) and writes out a
 * trajectory (in gromos or pdb format) or position restraints file in which
 * for every frame, the coordinates are only kept for atoms that are within a 
 * specific distance of a specified part of the system. To determine if 
 * interatomic distances are within the specified cut-off, either an atomic or a
 * charge-group based cut-off scheme can be employed. Additionally, parts of the
 * system can be specified for which in all cases, the atomic coordinates 
 * should either be kept or rejected.
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;boundary type&gt; [&lt;gather method&gt;] </td></tr>
 * <tr><td> \@atoms</td><td>&lt;@ref AtomSpecifier "atomspecifier": atoms to consider as reference part of the system&gt; </td></tr>
 * <tr><td> \@cutoff</td><td>&lt;cut-off distance (nm, default: 0.0)&gt; </td></tr>
 * <tr><td> \@select</td><td>&lt;@ref AtomSpecifier "atomspecifier": atoms to keep&gt; </td></tr>
 * <tr><td> \@reject</td><td>&lt;@ref AtomSpecifier "atomspecifier": atoms not to keep&gt; </td></tr>
 * <tr><td> \@pairlist</td><td>&lt;cut-off scheme (ATOMIC (default) or CHARGEGROUP)&gt; </td></tr>
 * <tr><td> \@outformat</td><td>&lt;g96, position, posres, posresspec or pdb&gt; </td></tr>
 * <tr><td> \@traj</td><td>&lt;input trajectory file(s)&gt; </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  filter
    @topo      ex.top
    @pbc       r
    @atoms     1:1
    @cutoff    0.5
    @select    1:1-20
    @reject    1:21-51
    @pairlist  ATOMIC
    @outformat pdb
    @traj      ex.tr
 @endverbatim
 *
 * <hr>
 */

#include <vector>
#include <iomanip>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cassert>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/gio/InG96.h"
#include "../src/gio/OutG96.h"
#include "../src/bound/Boundary.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/Bond.h"
#include "../src/gcore/Constraint.h"
#include "../src/gcore/Solvent.h"
#include "../src/gcore/SolventTopology.h"
#include "../src/gcore/Box.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/utils/SimplePairlist.h"
#include "../src/gio/InTopology.h"
#include "../src/gmath/Vec.h"

using namespace gcore;
using namespace gio;
using namespace args;
using namespace bound;
using namespace std;

enum outformatType { ofPdb, ofG96, ofG96red, ofPosres, ofPosresspec };

int main(int argc, char **argv){

  Argument_List knowns;
  knowns << "topo" << "pbc" << "traj" << "cutoff" << "atoms" << "select" 
         << "reject" << "pairlist" << "outformat";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo       <molecular topology file>\n";
  usage += "\t@pbc        <boundary type> [<gather method>]\n";
  usage += "\t[@atoms     <Atomspecifier atoms to consider as reference point>]\n";
  usage += "\t[@cutoff    <cutoff distance>]\n";
  usage += "\t[@select    <Atomspecifier atoms to keep>]\n";
  usage += "\t[@reject    <Atomspecifier atoms not to keep>]\n";
  usage += "\t[@pairlist  <ATOMIC or CHARGEGROUP>]\n";
  usage += "\t[@outformat <g96, position, posres, posresspec or pdb>]\n";
  usage += "\t@traj       <input trajectory files>\n";
   

  try{
    Arguments args(argc, argv, knowns, usage);

    // read topology
    InTopology it(args["topo"]);
    System sys(it.system());
    System osys(it.system());
    
    // parse boundary conditions
    Boundary *pbc = BoundaryParser::boundary(sys, args);
    //parse gather method
    Boundary::MemPtr gathmethod = args::GatherParser::parse(args);

    // define in and output coordinates
    InG96 ic;
    OutG96 oc;

    // read in the atom list the has to be kept definitely
    utils::AtomSpecifier ls(sys);
    {
      Arguments::const_iterator iter=args.lower_bound("select"),
	to=args.upper_bound("select");
      for(;iter!=to; ++iter)
	ls.addSpecifier(iter->second);
    }    
    // read in the atom list for atoms to throw away for sure
    utils::AtomSpecifier rej(sys);
    {
      Arguments::const_iterator iter=args.lower_bound("reject"),
	to=args.upper_bound("reject");
      for(;iter!=to; ++iter)
	rej.addSpecifier(iter->second);
    }
    // check that there are no doubles in ls and rej
    for(int i=0; i<rej.size(); i++)
      if( ls.findAtom(rej.mol(i), rej.atom(i))!=-1 )
	throw gromos::Exception("filter", "select and reject show overlap");

    // read in the reference atoms for the additional distance criterion
    utils::AtomSpecifier ref(sys);
    double cut=0;
    {
      Arguments::const_iterator iter=args.lower_bound("atoms"),
	to=args.upper_bound("atoms");
      for(;iter!=to; ++iter)
	ref.addSpecifier(iter->second);
    }
    if(ref.size()){
      // then we need a cutoff
      if(args.count("cutoff")<=0)
	throw gromos::Exception("filter", "If you specify reference atoms, then you need a cutoff");
      cut=atof(args["cutoff"].c_str());
    }

    // read in the type
    std::string t="ATOMIC";
    if(args.count("pairlist")>0){
      if(args["pairlist"]=="CHARGEGROUP") t=args["pairlist"];
      else if(args["pairlist"]!="ATOMIC") throw gromos::Exception("filter",
                "only ATOMIC and CHARGEGROUP are accepted for pairlist");
    }

    outformatType outformat=ofG96red;
    
    if(args.count("outformat")>0){
      if(args["outformat"]=="pdb")
        outformat = ofPdb;
      else if(args["outformat"]=="position")
        outformat = ofG96;
      else if(args["outformat"]=="g96")
        outformat = ofG96red;
      else if(args["outformat"]=="posres")
        outformat = ofPosres;
      else if(args["outformat"]=="posresspec")
        outformat = ofPosresspec;
      else throw gromos::Exception("filter", string("Unknown outformat: ") + args["outformat"]);
    }

    
    // loop over all trajectories
    ostream &os(cout);

    if(outformat != ofPdb){
      os <<"TITLE" << endl;
      
      os << "Filtered trajectory. Keeping" << endl;
      if(ls.size()){
	vector<string> s=ls.toString();
	os << "* atoms ";
	for(unsigned int i=0; i<s.size(); i++) os << s[i] << " ";
	os << endl;
      }
      if(ref.size()){
	vector<string> s=ref.toString();
	os << "* atoms within " << cut << " of atoms ";
	for(unsigned int i=0; i<s.size(); i++) os << s[i] << " ";
	os << endl;
	os << "  using a";
	if(t=="ATOMIC") os << "n atomic"; else os << " chargegroup based";
	os << " cutoff" << endl;
      }
      if(rej.size()){
	vector<string> s=rej.toString();
	os << "* as long as the above atoms do not belong to ";
	for(unsigned int i=0; i<s.size(); i++) os << s[i] << " ";
	os << endl;
      }
      os << "END" << endl;
    }
    
    for(Arguments::const_iterator iter=args.lower_bound("traj");
      iter!=args.upper_bound("traj"); ++iter){

      ic.open(iter->second);
      ic.select("ALL");
      oc.open(os);

      // loop over all frames
      while(!ic.eof()){
        ic >> sys;
	(*pbc.*gathmethod)();

	utils::AtomSpecifier rls=ls;
	
	Vec center(0.0,0.0,0.0);
	// add all atoms that need to be added according to the 
	// distances to reference atoms
        #pragma omp parallel for 
	for(int i=0; i<ref.size(); i++){
	  utils::SimplePairlist spl(sys, *pbc, cut);
	  spl.setAtom(*ref.atom()[i]);
	  spl.setType(t);
	  spl.calc();
	  
	  if((*ref.atom()[i]).type() != utils::spec_virtual)
	    spl.addAtom(ref.mol(i), ref.atom(i));
          
          #pragma omp critical 
          {
	    rls = rls + spl;
	    center += *ref.coord(i);
          }
	  
	}
	if(ref.size()) center/=ref.size();
	
	// remove atoms that are to be rejected
	for(int i=0; i<rej.size(); i++){
	  rls.removeAtom(rej.mol(i), rej.atom(i));
	}
	
	rls.sort();
        if(outformat == ofPdb) {
	  os.setf(ios::fixed, ios::floatfield);
	  os.setf(ios::unitbuf);
	  os.precision(3);
	  int res=0;
	  int count=0;
	  int resoff=0;
	  
	  for (int i=0; i<rls.size(); ++i){
	    
	    int maxmol=rls.mol(i);
	    if(maxmol<0) maxmol=sys.numMolecules();
	    count=rls.atom(i);
	    resoff=0;
	    for(int j=0; j< maxmol; j++) {
	      count+=sys.mol(j).numAtoms();
	      resoff+=sys.mol(j).topology().numRes();
	    }
	    
	    if(rls.mol(i)<0) res=rls.atom(i)/sys.sol(0).topology().numAtoms();
	    else res=sys.mol(rls.mol(i)).topology().resNum(rls.atom(i));
	    os << "ATOM";
	    os.setf(ios::right, ios::adjustfield);
	    os << setw(7) << count+1;
	    os.setf(ios::left, ios::adjustfield);
	    os << "  " <<setw(4) << rls.name(i).substr(0,3).c_str();
	    if(rls.mol(i)<0) os << setw(4) << "SOLV";
	    else os << setw(4) << sys.mol(rls.mol(i)).topology().resName(res).substr(0,4).c_str();
	    os.setf(ios::right, ios::adjustfield);
	    Vec pos=pbc->nearestImage(center, *rls.coord(i), sys.box()) 
	      - center;
	    
	    int resn = res+resoff+1;
	    if (resn > 9999) resn = 9999;
	    
	    os << " " << setw(4) << resn << "    "
	       << setw(8) << pos[0]*10
	       << setw(8) << pos[1]*10
	       << setw(8) << pos[2]*10
	       << "  1.00  0.00" << endl;
	  }
	  // now get the bonds
	  int molcount=0;
	  
	  for(int m=0; m<sys.numMolecules(); m++){
	    BondIterator bi(sys.mol(m).topology());
	    for(;bi;++bi){
	      int index_i=rls.findAtom(m,bi()[0]);
	      int index_j=rls.findAtom(m,bi()[1]);
	      int count_i=0, count_j=0;
	      
	      if(index_i!=-1 && index_j!=-1){
		  count_i=rls.atom(index_i)+molcount+1;
		  count_j=rls.atom(index_j)+molcount+1;

		  os << "CONECT " << count_i << " " << count_j << endl;
		  
	      }
	    }
	    molcount+=sys.mol(m).numAtoms();
	    
	  }
	  //also for solvent
	  int oldoffset=-1;
	  
	  for(int j=0; j<rls.size(); j++){
	    // check if it is a solvent
	    if(rls.mol(j)<0){
	      // get the atom number in this solvent
	      int si = rls.atom(j) % sys.sol(0).topology().numAtoms();
	      int offset=rls.atom(j) - si;
	      if(offset!=oldoffset){
		//new molecule
		ConstraintIterator ci(sys.sol(0).topology());
		for(; ci;++ci){
		  int index_i=rls.findAtom(rls.mol(j),offset+ci()[0]);
		  int index_j=rls.findAtom(rls.mol(j),offset+ci()[1]);
		  int count_i=0, count_j=0;
	      
		  if(index_i!=-1 && index_j!=-1){
		    count_i=rls.atom(index_i)+molcount+1;
		    count_j=rls.atom(index_j)+molcount+1;

		    os << "CONECT " << count_i << " " << count_j << endl;
		  }
		}
	      }
	      oldoffset=offset;
	    }
	  }		
	  os << "TER\n";
        } else if(outformat == ofG96 || outformat == ofPosres){
          if (outformat == ofG96)
            os  << "POSITION" << endl;
          else
            os << "POSRES" << endl;
          
	  os.setf(ios::fixed, ios::floatfield);
	  os.precision(9);
	  os << "# filter selected " << rls.size() << " atoms" << endl;
          os.setf(ios::unitbuf);
          
	  for (int i=0;i<rls.size();++i){
            os.setf(ios::right, ios::adjustfield);
            int offset = 1;
            for(int j=0;j < ((rls.mol(i) >= 0) ? rls.mol(i) : sys.numMolecules()); ++j)
              offset += sys.mol(j).topology().numRes();
            os << setw(5) << rls.resnum(i) + offset;
            os.setf(ios::left, ios::adjustfield);
            string res = rls.resname(i);
            
            if (rls.mol(i) < 0) res = "SOLV";
            os << ' ' <<setw(6) <<  res 
	       << setw(6) << rls.name(i);
            os.setf(ios::right, ios::adjustfield);
            os << setw(6) << rls.gromosAtom(i) + 1
	       << setw(15) << (*rls.coord(i))[0]
	       << setw(15) << (*rls.coord(i))[1]
	       << setw(15) << (*rls.coord(i))[2]<< endl;
	  }
	  os << "END" << endl;
        } else if (outformat == ofPosresspec) {
          os << "POSRESSPEC" << endl;
	  os.setf(ios::fixed, ios::floatfield);
	  os.precision(9);
	  os << "# filter selected " << rls.size() << " atoms" << endl;
          os.setf(ios::unitbuf);
          
	  for (int i=0;i<rls.size();++i){
            os.setf(ios::right, ios::adjustfield);
            int offset = 1;
            for(int j=0;j < ((rls.mol(i) >= 0) ? rls.mol(i) : sys.numMolecules()); ++j)
              offset += sys.mol(j).topology().numRes();
            os << setw(5) << rls.resnum(i) + offset;
            os.setf(ios::left, ios::adjustfield);
            string res = rls.resname(i);
            
            if (rls.mol(i) < 0) res = "SOLV";
            os << ' ' <<setw(6) <<  res 
	       << setw(6) << rls.name(i);
            os.setf(ios::right, ios::adjustfield);
            os << setw(6) << rls.gromosAtom(i) + 1 << endl;
	  }
	  os << "END" << endl;          
        } else {
	  os  << "POSITIONRED" << endl;
	  os.setf(ios::fixed, ios::floatfield);
	  os.precision(9);
	  os << "# filter selected " << rls.size() << " atoms" << endl;
	  for (int i=0;i<rls.size();++i){
	    os << setw(15) << (*rls.coord(i))[0]
	       << setw(15) << (*rls.coord(i))[1]
	       << setw(15) << (*rls.coord(i))[2] 
	       << "  # ";
	    if(rls.mol(i)<0) os << "s"; else os << rls.mol(i)+1;
	    os << ":" << rls.atom(i)+1 << endl;
	    if(!((i+1)%10))
	      os << "#" << setw(10) << i+1 << endl;
	  }
	  os << "END" << endl;
        }
        if (outformat != ofPdb) {
	  // and write the box block
          os << "GENBOX" << endl;
          oc.writeGenBox(sys.box());
          os << "END" << endl;
        }
      }    
      
      ic.close();
      oc.close();
    }
  }
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}
