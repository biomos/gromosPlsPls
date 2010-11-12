/**
 * @file prep_bb.cc
 * prepares a building block
 */

#include <map>


/**
 * @page contrib Contrib Program Documentation
 *
 * @anchor prep_bb
 * @section prep_bb prepares a building block
 * @author @ref co
 * @date 16. 3. 2005
 *
 * Program prep_bb takes information from a specification or PDB file and creates
 * a building block. If required, it makes suggestions for the choice of parameters.
 * 
 * The specification file has to be given in a certain format:
 * @verbatim
TITLE
ETOH
END
ATOMS
# specification of the atoms
# name    element
  C1      C
  C2      C
  O2      O
  H2      H
 END
 BONDS
 # specification of bonds
 # i    j    type (1=single, 2=double, 3=triple, 4=delocalized/aromatic)
   1    2    1
   2    3    1
   3    4    1
 END
 @endverbatim
 * The file consists of three blocks: a TITLE block which holds the name of the 
 * building block, an ATOMS block which specifies the names of the atoms and a 
 * BONDS block which specifies the bonds. Atomic elements and bond types are 
 * only required for graph based parameter suggestions (see below) and can be omitted. 
 *
 * If a PDB file is given bonds are automatically detected: all atoms pairs
 * with a interatomic distance lower than a given limit (\@bound) are considered
 * to be bonded.
 *
 * The order of the atoms in the resulting building block can be controled by
 * specifying the first atom using the \@reorder argument.
 *
 * If force-field files (and a graph library) are given, the 
 * program interacts (\@interact) with the user and provides advice for the choice of 
 * appropriate parameters. Else, the program chooses invalid default parameters
 * and writes a building block skeletton for further manual processing. 
 *
 * prep_bb is able the predict parameters using a graph based algorithm in a 
 * accurate way. In order to use this features you have to provide a specification
 * file (the information in a PDB is not sufficient) and give element and bond 
 * type information (see above). In addition, a graph library is needed to translate
 * the existing building blocks in the MTB files to a graph. The graph library should
 * look like this:
  * @verbatim
TITLE
Graph library for force field
END
FORCEFIELD
53A6
END
ELEMENTMAPPING
# mass  element (capitals only. i.e. Cl should be CL)
1       H
3       C
4       C
5       C
...
END
BONDMAPPING
# bond type   bond order 
# bond orders: 1.0 : single
#              1.5 : aromatic bond/delocalized
#              2.0 : double bond, 3.0: triple bond.
1   1.0
2   1.0
...
 END
 @endverbatim
 *
 * The resulting building block is written to BUILDING.out 
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@spc</td><td>&lt;specification file&gt; OR</td></tr>
 * <tr><td> \@pdb</td><td>&lt;PDB file&gt; </td></tr>
 * <tr><td> \@bound</td><td>&lt;upper bound for bond length (for PDB)&gt; </td></tr>
 * <tr><td> \@reorder</td><td>&lt;first atom&gt;</td></tr>
 * <tr><td>[\@build</td><td>&lt;building block file&gt;]</td></tr>
 * <tr><td>[\@param</td><td>&lt;parameter file&gt;]</td></tr>
 * <tr><td>[\@interact</td><td></td></tr>
 * <tr><td>[\@graph_library</td><td>&lt;library to translate BB into graphs.&gt;]</td></tr>
 * </table>
 *
 * Example:
 * @verbatim
   prep_bb 
     @pdb           ligand.pdb
     @build         mtb53a6.dat
     @param         ifp53a6.dat
     @interact
     @graph_library graphlib.53a6
   @endverbatim

 * <hr>
 */

#include <iostream>
#include <vector>
#include <set>
#include <string>
#include <iomanip>
#include <cassert>

#include "../src/args/Arguments.h"
#include "../src/gio/Ginstream.h"
#include "../src/gcore/BuildingBlock.h"
#include "../src/gcore/AtomPair.h"
#include "../src/gcore/GromosForceField.h"
#include "../src/gcore/Bond.h"
#include "../src/gcore/BondType.h"
#include "../src/gcore/Angle.h"
#include "../src/gcore/AngleType.h"
#include "../src/gcore/Improper.h"
#include "../src/gcore/ImproperType.h"
#include "../src/gcore/Dihedral.h"
#include "../src/gcore/DihedralType.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gcore/Exclusion.h"
#include "../src/gcore/LJException.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/BbSolute.h"
#include "../src/gio/InBuildingBlock.h"
#include "../src/gio/InParameter.h"
#include "../src/gio/InBuildingBlock.h"
#include "../src/gio/InParameter.h"
#include "../src/utils/FfExpert.h"
#include "../src/utils/FfExpertGraph.h"
#include "../src/gmath/Vec.h"
#include "../src/utils/CommandLine.h"
#include "../src/gio/OutBuildingBlock.h"

using namespace std;
using namespace args;
using namespace gcore;
using namespace gio;
using namespace utils::CommandLine;

static const string program_name = "prep_bb";

template<typename key_type, typename value_type>
map<value_type, key_type> swap_key_value(const map<key_type, value_type> & m) {
  map<value_type, key_type> res;
  for(typename map<key_type, value_type>::const_iterator it = m.begin(), to = m.end(); it != to; ++it) {
    res[it->second] = it->first;
  }
  return res;
}

string read_input(string file, vector<string> & atom, vector<string> & element, vector<Bond> & bond);
string read_pdb(string file, vector<string> &atom, vector<string> & element, vector<Bond> & bond, 
		double bondbound);
/*
string read_mol(string file, vector<string> &atom, vector<string> & element, vector<Bond> & bond);
*/
set<int> neighbours(int a, vector<Bond> & bond);
void forward_neighbours(int a, vector<Bond> & bond, set<int> &cum, int prev);
void add_neighbours(int i, vector<Bond> & bond, set<int> &at, vector<int> & order);
int count_bound(int i, vector<Bond> & bond, set<int> &at, set<int>& j);
set<int> ring_atoms(vector<Bond> & bond, int numatoms);

string & operator++(string & s)
{
  return s = s + "  ";
}

string & operator--(string & s)
{
  return s = s.substr(0, s.length() - 2);
}

int bignumber=1000;

int main(int argc, char **argv){
  Argument_List knowns;
  knowns << "spec" << "pdb" << "bound" << "reorder" << "build" << "param"
         << "interact" << "graph_library";

  string usage = "# " + program_name;
  usage += "\n\t@spec         <specifications> OR\n";
  usage += "\t@pdb            <pdb-file>\n";
  usage += "\t@bound          <upper bound for bondlength (for pdb)>\n";
  usage += "\t@reorder        <first atom>\n";
  usage += "\t[@build         <building block file>]\n";
  usage += "\t[@param         <parameter file>]\n";
  usage += "\t[@interact]\n";
  usage += "\t[@graph_library <file>]\n";
 
  try{
    Arguments args(argc, argv, knowns, usage);

    // check whether the user wants interaction
    bool interact=false;
    if(args.count("interact") >= 0) interact=true;
    
    // set up the expert system that knows everything
    utils::FfExpert exp;
    bool suggest = false;
    bool has_graph = false;
    if(args.count("build")>0){
      suggest = true;
      BuildingBlock mtb;
      Arguments::const_iterator iter=args.lower_bound("build"), to=args.upper_bound("build");
      for( ; iter!=to; ++iter){
	InBuildingBlock ibb(iter->second);
	mtb.addBuildingBlock(ibb.building());
      }
      
      if (args.count("graph_library") && !args.count("pdb")) {
        utils::FfExpertGraphMapper mapper(args["graph_library"]);
        exp.learn(mtb, &mapper);
        has_graph = true;
      } else {
        exp.learn(mtb);
      }
      if(!interact)
	throw gromos::Exception(program_name, "specification of building block "
				"and or parameter file for suggested "
				"parameters only takes effect if interact "
				"is specified");
      
    }
    GromosForceField *gff = NULL;
    if(args.count("param")>0){
      InParameter ipp(args["param"]);
      gff= new GromosForceField(ipp.forceField());
    }
    else if(suggest){
      throw gromos::Exception(program_name, "If you specify a building block "
			      "for suggestions, I also need a parameter file");
    }
    

    // read input    
    vector<string> atom_name, elements;
    vector<Bond> bond;
    string name;
    if(args.count("spec")>0)
      name=read_input(args["spec"], atom_name, elements, bond);
    else if(args.count("pdb")>0){
      // let's disable graph suggestions for now. PDB doesn't know about
      // bond order
      has_graph = false;
      double bondbound=0.0;
      if(args.count("bound")>0)	bondbound=atof(args["bound"].c_str());
      name=read_pdb(args["pdb"], atom_name, elements, bond, bondbound);
    }
    if(atom_name.size()>1 && bond.size()==0)
      throw gromos::Exception(program_name, "No bonds defining determined");
    
    set<int> at;
    vector<int> order;
    for(unsigned int i=0; i< atom_name.size(); i++)
      at.insert(i);

    // If necessary, get first atom and reorder the atoms
    if(args.count("reorder")>0){
      int first_atom=atoi(args["reorder"].c_str())-1;
      
      order.push_back(first_atom);
      at.erase(first_atom);
      add_neighbours(first_atom, bond, at, order);
    }
    else{
      for(unsigned int i=0; i< atom_name.size(); i++)
	order.push_back(i);
    }
    
    if(order.size()!=atom_name.size())
      throw gromos::Exception(program_name, "ordering of atoms failed");

    // do some bookkeeping
    map<int, int> newnum;
    for(unsigned int i=0; i<order.size(); i++){
      newnum[order[i]]=i;
    }
    vector<Bond> bondnew;
    for(unsigned int i=0; i< bond.size(); i++){
      Bond b(newnum[bond[i][0]], newnum[bond[i][1]]);
      bondnew.push_back(b);
    }
    
    if(interact && args.count("reorder") >0){
      cerr << "\n--------------------------------------------------------"
	   << "----------------------\n";
      cerr << "Atoms have been renumbered starting with atom "
	   << order[0]+1 << " (" << atom_name[order[0]] << ")\n";
      cerr << "New atom list:\n";
      for(unsigned int i=0; i<order.size(); i++){
	cerr << "\t" << setw(8) << i+1
	     << setw(8) << atom_name[order[i]]
	     << " (was " << order[i]+1 << ")\n";
      }
      cerr << "\n";
      cerr << "Bonds have been adopted accordingly\n";
      cerr << "\n========================================================"
	   << "======================\n";
    }
    
    // determine atoms that are in rings
    set<int> ra=ring_atoms(bondnew, order.size());
    set<int> aro;
    set<int> bt_aro;

    if(interact && ra.size()){
      cerr << "\n--------------------------------------------------------"
	     << "----------------------\n";
      cerr << "Bonds indicate ring system involving atoms ";
      for(set<int>::iterator it=ra.begin(), to=ra.end(); it!=to; ++it)
	cerr << *it + 1 << " ";
      cerr << "\n";
      if(getYesNo("Is this an aromatic ringsystem? (y/n) ")){
	for(set<int>::iterator it=ra.begin(), to=ra.end(); it!=to; ++it){
	  set<int> nb=neighbours(*it, bondnew);
	  if(nb.size()<4){
	    aro.insert(*it);
	    for(set<int>::iterator it2=nb.begin(), to2=nb.end(); it2!=to2; 
		++it2){
	      if(!ra.count(*it2)) bt_aro.insert(*it2);
	    }
	  }
	}
      } else {
	if(getYesNo("Do you want to specify different atoms as being part of an aromatic ring? (y/n) ")){
	  int number=0;
	  getValue<int>(number, "Give number of aromatic atoms: ");
	  cerr << "Give " << number << " atom numbers involving an "
	       << "aromatic ring: ";
	  for(int j=0; j< number; j++){
	    int n;
            ostringstream prompt;
            prompt << " - " << (j+1) << " atom: ";
	    getValue<int>(n, prompt.str());
	    aro.insert(n-1);
	  }
	  for(set<int>::iterator it=aro.begin(), to=aro.end(); it!=to; ++it){
	    set<int> nb=neighbours(*it, bondnew);
	    for(set<int>::iterator it2=nb.begin(), to2=nb.end(); it2!=to2; 
		++it2){
	      if(!aro.count(*it2)) bt_aro.insert(*it2);
	    }
	  }
	}
      }
      cerr << "\n\033[1;34m"
	   << "========================================================"
	   << "======================\033[22;0m\n";
    }

    // Now create a building block 
    ofstream fout("BUILDING.out");
    
    BbSolute bb;
    bb.setResName(name.substr(0,name.size()-1));

    std::vector<utils::FfExpert::counter> ocList;
    
    // now let's create the graph
    utils::FfExpertGraph * graph = NULL;
    if (has_graph)
      graph = new utils::FfExpertGraph(atom_name, elements, bond, order);

    double totcharge=0;
    int numchargegroup=0;
    for(unsigned int i=0; i< order.size(); i++){

      string aname=atom_name[order[i]].substr(0,1);
      int mass=0;
      int iac=0;
      double charge=0.0;
      int chargegroup=0;
      if(interact){
	cerr << "\n\033[1;34m"
	     << "--------------------------------------------------------"
	     << "----------------------\n";
	cerr << "Considering atom " << i+1 << " ("
	     << atom_name[order[i]] << ")\n\033[22;0m";

        if (has_graph) {
          vector<vector<utils::Vertex> > hits;

          exp.substructure2iac(order[i], *graph, hits);
          vector<map<int, double> > stat(hits.size());

          bool found = false;
          for (unsigned int radius = hits.size() - 1; radius != 0; --radius) {
            const unsigned int r = radius - 1;
            if (hits[r].size()) {
              found = true;
            }
            for (unsigned int hit = 0; hit < hits[r].size(); ++hit) {
              if (stat[r].find(hits[r][hit].iac) != stat[r].end())
                stat[r][hits[r][hit].iac] += 100.0 / hits[r].size();
              else
                stat[r][hits[r][hit].iac] += 100.0 / hits[r].size();
            }
          }

          if (found) {
            cerr << "\n"
                    << "Suggested Integer Atom Codes for atom (" << aname << ") "
                    << "based on substructure matching: \n";
            bool first = true;
            for (unsigned int radius = hits.size() - 1; radius != 0; --radius) {
              const unsigned int r = radius - 1;
              if (!hits[r].size())
                continue;

              cerr << setw(8) << radius << ": ";
              cerr.precision(2);
              cerr.setf(ios::fixed, ios::floatfield);

              if (first)
                cerr << "\033[1;31m";

              map<double, int> sorted_stat = swap_key_value(stat[r]);

              for (map<double, int>::reverse_iterator it = sorted_stat.rbegin(), to = sorted_stat.rend();
                      it != to; ++it) {
                cerr << (it->second + 1) << " : " << setw(5) << it->first << " %  ";
              }
              if (first) {
                first = false;
                cerr << "\033[22;0m";
              }

              for (unsigned int hit = 0; hit < hits[r].size(); ++hit) {
                if (hit % 6 == 0)
                  cerr << endl << setw(10) << " ";

                cerr << hits[r][hit].iac + 1 << " " << hits[r][hit].residue
                        << "[" << hits[r][hit].name << "]";
                if (hit != hits[r].size() - 1)
                  cerr << ", ";
              }
              cerr << endl;
            }
          } else {
            cerr << "\n"
                    << "No Integer Atom Code found based on topological similarity!\n";
          }
        }
        cerr.precision(4);
        
	// get IAC
	exp.name2iac(aname, ocList);
	if(ocList.size()){
	  unsigned int ap = utils::sort(ocList, true);
	  
	  cerr << "\n"
	       << "Suggested Integer Atom Code based on the first letter "
	       << "of the name (" << aname << ") :\n";
	  
	  for(unsigned int i=0; i< ocList.size(); ++i){
	    if(i==ap) cerr << "\033[1;31m";
	    cerr << setw(10) << ocList[i].type+1 << " :"
		 << setw(5) << gff->atomTypeName(ocList[i].type) 
		 << " occurs " << setw(3)  << ocList[i].occurence
		 << " times\n";
	    if(i==ap) cerr << "\033[22;0m";
	  }
	}
	ostringstream prompt;
        prompt << "Give IAC ( 1 - " << gff->numAtomTypeNames() << " ): ";
	getValue<int>(iac, prompt.str());
	iac--;
	if(iac<0 || iac >=gff->numAtomTypeNames()){
	  cerr << "\033[1;31mThis IAC ("<< iac+1 << ") is not "
	       << "defined in the given parameter file\033[22;0m\n";
	}
        
        if (has_graph) {
          // save the iac in the vertex
          graph->vertices()[order[i]].iac = iac;
        }

	// get mass
	exp.iac2mass(iac, ocList);
	if(ocList.size()){
	  unsigned int ap = utils::sort(ocList, true);

	  cerr << "\n" << "Suggested Mass Code based on the IAC "
	       << "of the atom (" << iac+1 << ") :\n";

	  for(unsigned int i=0; i< ocList.size(); i++){
	    if(i==ap) cerr << "\033[1;31m";
	    cerr << setw(10) << ocList[i].type+1 << " :"
		 << setw(8) << gff->findMass(ocList[i].type) << " occurs " 
		 << setw(3) << ocList[i].occurence << " times\n";
	    if(i==ap) cerr << "\033[22;0m";
	  }
	}
	getValue<int>(mass, "Give Masstype: ");
	mass--;
	if(gff->findMass(mass)==0.0)
	  cerr << "\033[1;31mThis Masstype (" << mass+1 << ") is "
	       << "not defined in the parameter file\n";
        
        // charge
        if (has_graph) {
          vector<vector<utils::Vertex> > hits;
          exp.substructure2iac(order[i], *graph, hits);
          vector<map<double, double> > stat(hits.size());

          bool found = false;
          for (unsigned int radius = hits.size() - 1; radius != 0; --radius) {
            const unsigned int r = radius - 1;
            if (hits[r].size()) {
              found = true;
            }
            for (unsigned int hit = 0; hit < hits[r].size(); ++hit) {
              if (stat[r].find(hits[r][hit].charge) != stat[r].end())
                stat[r][hits[r][hit].charge] += 100.0 / hits[r].size();
              else
                stat[r][hits[r][hit].charge] += 100.0 / hits[r].size();
            }
          }

          if (found) {
            cerr << "\n"
                    << "Suggested charge for atom (" << aname << ") "
                    << "based on substructure matching: \n";
            bool first = true;
            for (unsigned int radius = hits.size() - 1; radius != 0; --radius) {
              const unsigned int r = radius - 1;
              if (!hits[r].size())
                continue;

              cerr << setw(8) << radius << ": ";
              cerr.setf(ios::fixed, ios::floatfield);

              if (first)
                cerr << "\033[1;31m";

              map<double, double> sorted_stat = swap_key_value(stat[r]);
              for (map<double, double>::reverse_iterator it = sorted_stat.rbegin(), to = sorted_stat.rend();
                      it != to; ++it) {
                cerr.precision(4);
                cerr << it->second;
                cerr.precision(2);
                cerr << " : " << setw(5) << it->first << " %  ";
              }
              if (first) {
                first = false;
                cerr << "\033[22;0m";
              }

              cerr.precision(3);
              for (unsigned int hit = 0; hit < hits[r].size(); ++hit) {
                if (hit % 6 == 0)
                  cerr << endl << setw(10) << " ";

                cerr << setw(5) << hits[r][hit].charge << " " << hits[r][hit].residue
                        << "[" << hits[r][hit].name << "]";
                if (hit != hits[r].size() - 1)
                  cerr << ", ";
              }
              cerr << endl;
            }
          } else {
            cerr << "\n" << "No charge found based on topological similarity!\n";
          }
        }

        cerr.precision(4);
        
	// get charge
	exp.iac2charge(iac, ocList);
	if(ocList.size()){
	  unsigned int ap = utils::sort(ocList, false);

	  cerr << "\n" << "Suggested Charge based on the IAC "
	       << "of the atom (" << iac+1 << ") :\n";

	  for(unsigned int i=0; i< ocList.size(); ++i){
	    if(i==ap) cerr <<"\033[1;31m";
	    cerr << setw(10) << exp.charge(ocList[i].type)
		 << " occurs " << setw(3) << ocList[i].occurence 
		 << " times\n";
	    if(i==ap) cerr << "\033[22;0m";
	  }
	  if(numchargegroup){
	    cerr << "Charge group so far (" << numchargegroup
		 << " atom";
	    if(numchargegroup > 1) cerr << "s";
	    cerr << ") has net charge of " << totcharge << "\n\n";
	  }
	}
	getValue<double>(charge, "Give Charge: ");
	totcharge+=charge;
	numchargegroup++;
	
	cerr << "\n"
	     << "Suggested Chargegroup code based on the total charge "
	     << "of this\n" << "charge group (" << totcharge 
	     << "; " << numchargegroup << " atom";
	if(numchargegroup>1) cerr << "s";
	cerr << "): ";
	if((int(rint(totcharge*1000))%1000) == 0)
	  cerr << "1\n\n";
	else
	  cerr << "0\n\n";

        std::set<int> zeroone;
        zeroone.insert(0); zeroone.insert(1);
        getValueFromSet<int>(chargegroup, zeroone, "Give Chargegroup code (0/1): ");

	if(chargegroup==1){
	  numchargegroup=0;
	  totcharge=0.0;
	}
	
      }
      // done gathering data. Prepare atom
      AtomTopology a_top;
      a_top.setName(atom_name[order[i]]);
      a_top.setMass(mass);
      a_top.setIac(iac);
      a_top.setCharge(charge);
      a_top.setChargeGroup(chargegroup);
      bb.setResNum(i,0);
      // It is a building block -> 14 exclusions are not strictly necessary
      // but in the case of aromaticity we might need it
      set<int> exclusions;
      set<int> exclusions14;
      set<int> nb1=neighbours(order[i], bond);
      set<int> nb2;
      set<int> nb3;
      set<int> tmp;
      
      for(set<int>::iterator it=nb1.begin(), to=nb1.end(); it!=to; ++it){
	exclusions.insert(newnum[*it]);
	nb2=neighbours(*it, bond);
	for(set<int>::iterator it2=nb2.begin(), to2=nb2.end(); it2!=to2; ++it2){
	  exclusions.insert(newnum[*it2]);
	  tmp=neighbours(*it2, bond);
	  for(set<int>::iterator it3=tmp.begin(), to3=tmp.end();
	      it3!=to3; ++it3){
	    nb3.insert(*it3);
	  }
	}
      }
      for(set<int>::iterator it=nb3.begin(), to=nb3.end(); it!=to; ++it){
	if(!exclusions.count(newnum[*it])){
	  exclusions14.insert(newnum[*it]);
	}
      }
      Exclusion ex, ex14;
      for(set<int>::iterator it=exclusions.begin(), to=exclusions.end();
	  it!=to; ++it)
	if(*it > int(i)) ex.insert(*it);
	
      for(set<int>::iterator it=exclusions14.begin(), to=exclusions14.end();
	  it!=to; ++it){
	if((aro.count(i)   || bt_aro.count(i)) && 
	   (aro.count(*it) || bt_aro.count(*it))){
	  
	  if(*it > int(i)) ex.insert(*it);
	}
	
	else
	  if(*it > int(i)) ex14.insert(*it);
      }
      a_top.setExclusion(ex);
      a_top.setExclusion14(ex14);
      if(interact)
	cerr << "\n\tDetermined exclusions and 1,4-exclusions.\n";
      
      bb.addAtom(a_top);
      
    }
    
    if(interact){
      
      cerr << "\n\033[1;34m"
	   << "======= ATOMIC INFORMATION GATHERED ===================="
	   << "======================\033[22;0m\n";
      
      cerr << "\033[1;34m======= IS WRITTEN TO FILE ============================="
	   << "======================\n\n\n";
      
      cerr << "======= BOND PARAMETERS ================================"
	   << "======================\033[22;0m\n\n";
    }

    for(unsigned int i=0; i< bondnew.size(); i++){
      if(interact){
	
	Bond iacbond(bb.atom(bondnew[i][0]).iac(), 
		     bb.atom(bondnew[i][1]).iac());
	cerr << "\n\033[1;34m"
	     << "--------------------------------------------------------"
	     << "----------------------\n";
	cerr << "Considering bond " << bondnew[i][0]+1 << " ("
	     << bb.atom(bondnew[i][0]).name() << ")  -  "
	     << bondnew[i][1]+1 << " (" 
	     << bb.atom(bondnew[i][1]).name()
	     << ")\033[22;0m\n";

	exp.iac2bond(iacbond, ocList);
	if(ocList.size()){
	  unsigned int ap=utils::sort(ocList);
	  cerr << "\n" << "Suggested Bondtype based on the IAC "
	       << "of the atoms (" << iacbond[0]+1 << ", " 
	       << iacbond[1]+1 << ") :\n";
	  for(unsigned int i=0; i < ocList.size(); ++i){
	    if(ap==i) cerr << "\033[1;31m";
	    cerr << setw(10) << ocList[i].type +1
		 << " : Kb = " << setw(10) << setprecision(1) 
		 << gff->bondType(ocList[i].type).fc()
		 << "; b0 = " << setw(7) << setprecision(3)
		 << gff->bondType(ocList[i].type).b0()
		 << " occurs " << setw(3) << ocList[i].occurence 
		 << " times\n";
	    if(ap==i) cerr << "\033[22;0m";
	  }
	}
	ostringstream prompt;
        prompt << "Give Bondtype ( 1 - " << gff->numBondTypes()
	     << " ) : ";
	int bt;
	getValue<int>(bt, prompt.str());
	if(bt<1 || bt >gff->numBondTypes()){
	  cerr << "\033[1;31mThis Bondtype ("<< bt << ") is not "
	       << "defined in the given parameter file\033[22;0m\n";
	}
	bondnew[i].setType(bt-1);
      }
      bb.addBond(bondnew[i]);
    }
    //for the angles and dihedrals we also first put them in a vector
    vector<Angle> newangles;
    vector<Improper> newimpropers;
    vector<Dihedral> newdihedrals;
    
    for(unsigned int i=0; i< order.size(); i++){
      set<int> nb=neighbours(i, bondnew);
      vector<int> vnb;
      for(set<int>::iterator it=nb.begin(), to=nb.end(); it!=to; ++it){
	vnb.push_back(*it);
      }
      switch(vnb.size()){
	case 1: break;
	case 2: 
	  newangles.push_back(Angle(vnb[0],i,vnb[1]));
	  break;
	case 3:
	  newangles.push_back(Angle(vnb[0],i,vnb[1]));
	  newangles.push_back(Angle(vnb[0],i,vnb[2]));
	  newangles.push_back(Angle(vnb[1],i,vnb[2]));
	  newimpropers.push_back(Improper(i, vnb[0],vnb[1],vnb[2]));
	  break;
	case 4:
	  newangles.push_back(Angle(vnb[0],i,vnb[1]));
	  newangles.push_back(Angle(vnb[0],i,vnb[2]));
	  newangles.push_back(Angle(vnb[0],i,vnb[3]));
	  newangles.push_back(Angle(vnb[1],i,vnb[2]));
	  newangles.push_back(Angle(vnb[1],i,vnb[3]));
	  newangles.push_back(Angle(vnb[2],i,vnb[3]));
	  break;
	case 5:
	  newangles.push_back(Angle(vnb[0],i,vnb[1]));
	  newangles.push_back(Angle(vnb[0],i,vnb[2]));
	  newangles.push_back(Angle(vnb[0],i,vnb[3]));
	  newangles.push_back(Angle(vnb[0],i,vnb[4]));
	  newangles.push_back(Angle(vnb[1],i,vnb[2]));
	  newangles.push_back(Angle(vnb[1],i,vnb[3]));
	  newangles.push_back(Angle(vnb[1],i,vnb[4]));
	  newangles.push_back(Angle(vnb[2],i,vnb[3]));
	  newangles.push_back(Angle(vnb[2],i,vnb[4]));
	  newangles.push_back(Angle(vnb[3],i,vnb[4]));
	  break;
	default:
	  ostringstream os;
	  os << "Don't know how to create angles for 0 or 6 bonds to atom "
	     << i+1;
	  
	  throw gromos::Exception(program_name, os.str());
	  
      }
    }
    
    for(unsigned int i=0; i< bondnew.size(); i++){
      set<int> nb1=neighbours(bondnew[i][0], bondnew);
      nb1.erase(bondnew[i][1]);
      set<int> nb2=neighbours(bondnew[i][1], bondnew);
      nb2.erase(bondnew[i][0]);
      
      if(nb1.size() && nb2.size()){
	if(aro.count(bondnew[i][0]) && aro.count(bondnew[i][1])){
	  

	  // we are in an aromatic ring
	  int a=-1, b=-1;
	  for(set<int>::iterator it=nb1.begin(), to=nb1.end();
	      it!=to; ++it){
	    if(aro.count(*it)) a=*it;
	  }
	  for(set<int>::iterator it=nb2.begin(), to=nb2.end();
	      it!=to; ++it){
	    if(aro.count(*it)) b=*it;
	  }
	  if(a==-1 || b==-1)
	  
	    throw gromos::Exception(program_name, "Error trying to determine "
				    "improper dihedrals for aromatic ring");
	  newimpropers.push_back(Improper(a, bondnew[i][0], bondnew[i][1], b));
	  
	}
	else
	  
	  newdihedrals.push_back(Dihedral(*nb1.begin(), bondnew[i][0], 
					 bondnew[i][1], *nb2.begin()));
      }
      
    }
    // now get the types
    if(interact && newangles.size())
      cerr << "\n\n\033[1;34m"
	   << "======= ANGLE PARAMETERS ==============================="
	   << "======================\033[22;0m\n\n";
    for(unsigned int i=0; i< newangles.size(); i++){
      if(interact){
	Angle iacangle(bb.atom(newangles[i][0]).iac(), 
		       bb.atom(newangles[i][1]).iac(),
		       bb.atom(newangles[i][2]).iac());
	cerr << "\n\033[1;34m"
	     << "--------------------------------------------------------"
	     << "----------------------\n";
	cerr << "Considering angle " << newangles[i][0]+1 << " ("
	     << bb.atom(newangles[i][0]).name() << ")  -  "
	     << newangles[i][1]+1 << " (" 
	     << bb.atom(newangles[i][1]).name()
	     << ")  -  "  << newangles[i][2]+1 << " ("
	     << bb.atom(newangles[i][2]).name() << ")\033[22;0m\n";
	
	exp.iac2angle(iacangle, ocList);
	if(ocList.size()){
	  unsigned int ap= utils::sort(ocList);
	  
	  cerr << "\n" << "Suggested Angletype based on the IAC "
	       << "of the atoms (" << iacangle[0]+1 << ", " 
	       << iacangle[1]+1 << ", " << iacangle[2]+1 << ") :\n";

	  for(unsigned int i=0; i< ocList.size(); ++i){
	    if(ap==i) cerr << "\033[1;31m";
	    cerr << setw(10) << ocList[i].type +1
		 << " : Kb = " << setw(10) << setprecision(1) 
		 << gff->angleType(ocList[i].type).fc()
		 << "; t0 = " << setw(8) << setprecision(3)
		 << gff->angleType(ocList[i].type).t0()
		 << " occurs " 
		 << ocList[i].occurence << " times\n";
	    if(ap==i) cerr << "\033[22;0m";
	  }
	  
	}
	ostringstream prompt;
        prompt << "Give Angletype ( 1 - " << gff->numAngleTypes()
	     << " ) : ";
	int bt;
	getValue<int>(bt, prompt.str());
	if(bt<1 || bt >gff->numAngleTypes()){
	  cerr << "\033[1;31mThis Angletype ("<< bt << ") is not "
	       << "defined in the given parameter file\033[22;0m\n";
	}
	newangles[i].setType(bt-1);
      }
      bb.addAngle(newangles[i]);
    }    
    if(interact && newimpropers.size())
      cerr << "\n\n\033[1;34m"
	   << "======= IMPROPER PARAMETERS ============================"
	   << "======================\033[22;0m\n\n";
    for(unsigned int i=0; i< newimpropers.size(); i++){
      if(interact){
	Improper iacimproper(bb.atom(newimpropers[i][0]).iac(), 
			     bb.atom(newimpropers[i][1]).iac(),
			     bb.atom(newimpropers[i][2]).iac(),
			     bb.atom(newimpropers[i][3]).iac());
	
	cerr << "\n\033[1;34m"
	     << "--------------------------------------------------------"
	     << "----------------------\n";
	cerr << "Considering improper " << newimpropers[i][0]+1 << " ("
	     << bb.atom(newimpropers[i][0]).name() << ")  -  "
	     << newimpropers[i][1]+1 << " (" 
	     << bb.atom(newimpropers[i][1]).name()
	     << ")  -  " << newimpropers[i][2]+1 << " (" 
	     << bb.atom(newimpropers[i][2]).name() 
	     << ")  -  " << newimpropers[i][3]+1 << " (" 
	     << bb.atom(newimpropers[i][3]).name() 
	     << ")\033[22;0m\n";
	
	exp.iac2improper(iacimproper, ocList);
	if(ocList.size()){
	  unsigned int ap = utils::sort(ocList);
	  
	  cerr << "\n" << "Suggested Impropertype based on the IAC "
	       << "of the atoms (" << iacimproper[0]+1 << ", " 
	       << iacimproper[1]+1 << ", " << iacimproper[2]+1 << ", " 
	       << iacimproper[3]+1 << ") :\n";
	  for(unsigned int i=0; i< ocList.size(); ++i){
	    if(i==ap) cerr << "\033[1;31m";
	    cerr << "\t\t" << setw(10) << ocList[i].type +1
		 << " : Kb = " << setw(10) << setprecision(1) 
		 << gff->improperType(ocList[i].type).fc()
		 << "; q0 = " << setw(10) << setprecision(3)
		 << gff->improperType(ocList[i].type).q0()
		 << " occurs " 
		 << ocList[i].occurence << " times\n";
	    if(i==ap) cerr << "\033[22;0m";
	  }
	}
	ostringstream prompt;
        prompt << "Give Impropertype ( 1 - "
	     << gff->numImproperTypes() 
	     << " ) : ";
	int bt;
	getValue<int>(bt, prompt.str());
	if(bt<1 || bt >gff->numImproperTypes()){
	  cerr << "\033[1;31mThis Impropertype ("<< bt
	       << ") is not defined in the given "
	       << "parameter file\033[22;0m\n";
	}
	if(gff->improperType(bt-1).q0()!=0){
	  cerr << "\033[1;31mYou might want to modify the order of "
	       << " the atoms for this improper dihdral\n" 
	       << "to make sure you get the correct stereochemistry\n"
	       << "\033[22;0m";
	}
	newimpropers[i].setType(bt-1);
      }
      bb.addImproper(newimpropers[i]);
    }   
    if(interact && newdihedrals.size())
      cerr << "\n\n\033[1;34m"
	   << "======= DIHEDRAL PARAMETERS ============================"
	   << "======================\033[22;0m\n\n";
    for(unsigned int i=0; i< newdihedrals.size(); i++){
      if(interact){
	Dihedral iacdihedral(bb.atom(newdihedrals[i][0]).iac(), 
			     bb.atom(newdihedrals[i][1]).iac(),
			     bb.atom(newdihedrals[i][2]).iac(),
			     bb.atom(newdihedrals[i][3]).iac());
	cerr << "\n\033[1;34m"
	     << "--------------------------------------------------------"
	     << "----------------------\n";
	cerr << "Considering dihedral "
	     << newdihedrals[i][0]+1 
	     << " (" << bb.atom(newdihedrals[i][0]).name() << ")  -  "
	     << newdihedrals[i][1]+1 
	     << " (" << bb.atom(newdihedrals[i][1]).name() << ")  -  "
	     << newdihedrals[i][2]+1
	     << " (" << bb.atom(newdihedrals[i][2]).name() << ")  -  "
	     << newdihedrals[i][3]+1
	     << " (" << bb.atom(newdihedrals[i][3]).name() << ")\033[22;0m\n";
	
	exp.iac2dihedral(iacdihedral, ocList);
	if(ocList.size()){
	  unsigned int ap=utils::sort(ocList);
	  
	  cerr << "\n" << "Suggested Dihedraltype based on the IAC "
	       << "of the atoms (" 
	       << iacdihedral[0]+1 << ", " 
	       << iacdihedral[1]+1 << ", "
	       << iacdihedral[2]+1 << ", "
	       << iacdihedral[3]+1 << ") :\n";
	  for(unsigned int i=0; i< ocList.size(); ++i){
	    if(ap==i) cerr << "\033[1;31m";
	    cerr << setw(10) << ocList[i].type +1
		 << " : Kd = " << setw(6) << setprecision(1) 
		 << gff->dihedralType(ocList[i].type).fc()
		 << "; pd = " << setw(4) << setprecision(1)
		 << gff->dihedralType(ocList[i].type).pd()
		 << "; m = " << setw(2) << setprecision(1)
		 << gff->dihedralType(ocList[i].type).np()
		 << " occurs " 
		 << ocList[i].occurence << " times\n";	
	    if(ap==i) cerr << "\033[22;0m" ;
	  }
	}
        ostringstream prompt;
	prompt << "Give Dihedraltype ( 1 - "
	     << gff->numDihedralTypes()  << " ) : ";
	int bt;
	getValue<int>(bt, prompt.str());
	if(bt<1 || bt >gff->numDihedralTypes()){
	  cerr << "\033[1;31mThis Dihedraltype ("<< bt << ") is not defined "
	       << "in the given parameter file\033[22;0m\n";
	}
	newdihedrals[i].setType(bt-1);
      }
      bb.addDihedral(newdihedrals[i]);
      
    }
    if(interact)
      cerr << "\n\n\033[1;34m"
	   << "======= BONDED PARAMETERS COLLECTED ======================="
	   << "===================\033[22;0m\n\n";

    OutBuildingBlock obb(fout);
    obb.writeSingle(bb, OutBuildingBlock::BBTypeSolute);
    fout.close();
    cerr << "Building block was written to BUILDING.out" << endl;
    if(interact)
      cerr << "\n\n\033[1;34m"
	   << "======= PREPBBB HAS GATHERED ALL INFORMATION AND WRITTEN"
	   << " A BUILDING BLOCK=====\033[22;0m\n\n";
    
    
  }
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}

string read_input(string file, vector<string> & atom, vector<string> & element, vector<Bond> & bond)
{
  vector<string> buffer;
  int a,b;
  
  Ginstream gin(file);
  gin.getblock(buffer);
  if(buffer[buffer.size()-1].find("END")!=0)
      throw gromos::Exception(program_name,"Input file " + gin.name() +
			      " is corrupted. No END in "+buffer[0]+
			      " block. Got\n"
			      + buffer[buffer.size()-1]);
  if(buffer[0]!="ATOMS")
    throw gromos::Exception(program_name, "ATOMS block expected in file "+file);
  for(unsigned int i=1; i< buffer.size()-1; i++){
    std::string name, elem;
    istringstream is(buffer[i]);
    is >> name;
    atom.push_back(name);
    if (is >> elem)
      element.push_back(elem);
  }
  gin.getblock(buffer);
  if(buffer[0]!="BONDS")
    throw gromos::Exception(program_name, "BONDS block expected in file "+file);
  for(unsigned int i=1; i< buffer.size()-1; i++){
    istringstream is(buffer[i]);
    is >> a >> b;
    if (is.fail())
      throw gromos::Exception(program_name, "bad bond in BOND block");
    bond.push_back(Bond(a-1,b-1));
    int type;
    if (is >> type)
      bond.back().setType(type);
  }
  
  return gin.title();
}
string read_pdb(string file, vector<string> &atom, vector<string> & element, vector<Bond> & bond, 
		double bondbound)
{
  ifstream fin(file.c_str());
  string inPdbLine;
  string resName;
  double bondbound2=bondbound*bondbound;
  
  vector<gmath::Vec> pos;
  
  while(!fin.eof()){
    getline(fin, inPdbLine);
    if(inPdbLine.substr(0,4) == "ATOM" ||
       inPdbLine.substr(0,6) == "HETATM"){
      istringstream is(inPdbLine.substr(12, 5));
      string name;
      is >> name;
      atom.push_back(name);
      resName = inPdbLine.substr(17, 4);
      gmath::Vec v;
      v[0]=atof(inPdbLine.substr(30, 8).c_str());
      v[1]=atof(inPdbLine.substr(38, 8).c_str());
      v[2]=atof(inPdbLine.substr(46, 8).c_str());
      pos.push_back(v);
      
      // get the element and remove white space
      string elementStr = inPdbLine.substr(77, 2).c_str();
      string::size_type i;
      while((i = elementStr.find(" ")) != string::npos )
        elementStr.erase(i);
      
      element.push_back(elementStr);
    }
    if(inPdbLine.substr(0,6) == "CONECT"){
      istringstream is(inPdbLine.substr(6,80));
      int i,j;
      is >> i >> j;
      bond.push_back(Bond(i-1,j-1));
    }
  }
  if(bond.size()==0){
    for(unsigned int i=0; i< atom.size(); i++){
      for(unsigned int j=i+1; j< atom.size(); j++){
	if((pos[i]-pos[j]).abs2()< bondbound2){
	  bond.push_back(Bond(i,j));
	}
      }
    }
  }

  return resName;
}

/*
string read_mol(string file, vector<string> &atom, vector<string> & element, vector<Bond> & bond)
{
  ifstream fin(file.c_str());
  vector<string> buffer;
  for(unsigned int i = 0; !fin.eof(); ++i) {
    string line;
    getline(fin, line);
    if (i > 3) // skip the header
      buffer.push_back(line);
  }
  fin.close();
  
  vector<string>::const_iterator it = buffer.begin(),
          to = buffer.end();
  
  istringstream _lineStream(*it);
  unsigned int num_atoms, num_bonds;
  _lineStream >> num_atoms >> num_bonds;
  if (_lineStream.fail())
    throw gromos::Exception(program_name, "Could not read number of atoms/bonds from molfile.");
  
  for(unsigned int i = 0; i < num_atoms; ++i) {
    _lineStream.str(*(++it));
    double pos;
    std::string e;
    _lineStream >> pos >> pos >> pos >> e; // discard the position information
    if (_lineStream.fail())
      throw gromos::Exception(program_name, "bad atom in molfile.");
    
    element.push_back(e);
    atom.push_back(e + (i + 1));
  }
  
  for(unsigned int i = 0; i < num_bonds; ++i) {
    _lineStream.str(*(++it));
    unsigned int i, j, t;
    _lineStream >> i >> j >> t;
    if (_lineStream.fail())
      throw gromos::Exception(program_name, "bad bond in molfile.");
    i--; j--;
    if (i < 0 || i >= atom.size() || j < 0 || j >= atom.size())
      throw gromos::Exception(program_name, "bad atoms in bond in molfile.");

    Bond bond(i,j);
    bond.setType(t);
    element.push_back(bond);
  }
  return "CMPD";
}
*/

set<int> neighbours(int a, vector<Bond> & bond)
{
  set<int> tmp;
  for(unsigned int i=0; i< bond.size(); i++){
    if(bond[i][0]==a) tmp.insert(bond[i][1]);
    if(bond[i][1]==a) tmp.insert(bond[i][0]);
  }
  return tmp;
}

void forward_neighbours(int a, vector<Bond> & bond, set<int> & cum, int prev)
{
  set<int> tmp=neighbours(a, bond);
  tmp.erase(prev);
  for(set<int>::iterator it=tmp.begin(), to=tmp.end(); it!=to; ++it){
    //cout << "\tadding " << *it << " (from " << a << ")" << endl;
    if(cum.count(*it)){
      
      //cout << " is this a ring! atom? " << *it << endl;
      return;
    }
    
    else{
      cum.insert(*it);
      forward_neighbours(*it, bond, cum, a);
    }
  }
}
  

set<int> ring_atoms(vector<Bond> & bond, int numatoms)
{
  set<int> ra;
  for(int i=0; i< numatoms; i++){
    //cout << "analysing ring-possibility for " << i << endl;
    
    set<int> tmp;
    forward_neighbours(i,bond, tmp,-1);
    if(tmp.count(i)) ra.insert(i);
  }
  return ra;
  
    
}

  
void add_neighbours(int i, vector<Bond> & bond, set<int> &at, vector<int> & order)
{
  //cout << "called add_neighbours with i: " << i << endl;
  
  set<int> nb = neighbours(i, bond);
  //cout << "\t" << i << " has " << nb.size() << " neighbours" << endl;
  
  map<int,int> cb;
  set<int>::const_iterator it=nb.begin(), to=nb.end();
  for( ; it!=to; ++it){
    //cout << "\t\tcounting the bonds for every neighbour" << endl;
    set<int> counted;
    counted.insert(i);
    
    if(at.count(*it))
      cb[*it]=count_bound(*it, bond, at, counted);
    else
      cb[*it]=bignumber;
    
    //cout << "\t\t" << *it << "\t" << cb[*it] << endl;
  }
  for(map<int,int>::iterator mit=cb.begin(), mto=cb.end(); mit!=mto; ++mit){
    if(mit->second == bignumber )
      nb.erase(mit->first);
  }
  //cout << "\tcorrected neighbour size " << nb.size() << endl;
  
  while(nb.size()){
    int min=bignumber;
    int nextatom=-1;
    
    for(it=nb.begin(), to=nb.end(); it!=to; ++it){
      if(cb[*it] < min) {
	min=cb[*it];
	nextatom=*it;
      }
    }
    if(nextatom!=-1){
      if(at.count(nextatom)){
	 nb.erase(nextatom);
	 at.erase(nextatom);
	 //cout << "added atom " << nextatom << endl;
	 order.push_back(nextatom);
	 add_neighbours(nextatom, bond, at, order);
      }
      else{
	nb.erase(nextatom);
	//cout << "ring! atom was gone already" << i << " - " << nextatom << endl;
      }
      
    }
    
  }
}

int count_bound(int i, vector<Bond> & bond, set<int> &at, set<int> & j)
{
  set<int> nb = neighbours(i, bond);
  int counter=0;
  set<int>::const_iterator it=nb.begin(), to=nb.end();
  for(; it!=to; ++it){
    if(at.count(*it) && !j.count(*it)){
      counter++;
      j.insert(i);
      counter+=count_bound(*it, bond, at, j);
    }
  }
  return counter;
}
