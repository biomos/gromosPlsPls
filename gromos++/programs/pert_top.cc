/**
 * @file pert_top.cc
 * Creates a perturbation topology to remove interactions for specified atoms
 */

/**
 * @page programs Program Documentation
 *
 * @anchor pert_top
 * @section pert_top Creates a perturbation topology to alter interactions for 
 * specified atoms
 * @author @ref co @ref ns
 * @date 7-6-07
 *
 * Creates a perturbation topology to perturb specified atoms. A perturbation 
 * topology is written that defines a perturbation to alter the specified atoms
 * to the specified atom types, charges and masses. Each of the arguments
 * \@types, \@masses and \@charges can be omitted. In this case the values from
 * the topology are taken. If not sufficient values are given, the last given 
 * value is taken for all the remaining atoms.

 * Use program @ref pt_top to convert the resulting perturbation topology to a
 * different format or to a regular molecular topology.
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@atoms</td><td>&lt;@ref AtomSpecifier: atoms to be modified </td></tr>
 * <tr><td> \@types</td><td>&lt;IACS of the perturbed atoms&gt; </td></tr>
 * <tr><td> \@charges</td><td>&lt;charges of the perturbed atoms&gt; </td></tr>
 * <tr><td> \@masses</td><td>&lt;masses of the perturbed atoms&gt; </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  pert_top
    @topo    ex.top
    @atoms   1:34-51
    @types   13 19
    @charges 0.1 0.0
    @masses  1.008
 @endverbatim
 *
 * <hr>
 */
#include <cassert>
#include <map>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>
#include "../src/args/Arguments.h"
#include "../src/gio/InTopology.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/System.h"
#include "../src/utils/AtomSpecifier.h"

using namespace std;
using namespace gcore;
using namespace gio;
using namespace args;


int main(int argc, char *argv[]){

  Argument_List knowns;
  knowns << "topo" << "atoms" << "types" << "charges" << "masses";
  
  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo  <molecular topology file>\n";
  usage += "\t@atoms   <AtomSpecifier: atoms to be modified\n";
  usage += "\t@types   <IACS of the perturbed atoms>\n";
  usage += "\t@charges <charges of the perturbed atoms>\n";
  usage += "\t@masses  <masses of the perturbed atoms>\n";
   
  try{
    Arguments args(argc, argv, knowns, usage);

    InTopology it(args["topo"]);
    System sys(it.system());
    
    // get the atoms to perturb
    utils::AtomSpecifier at(sys);
    {
      Arguments::const_iterator iter=args.lower_bound("atoms");
      for(; iter!=args.upper_bound("atoms"); ++iter){
	at.addSpecifier(iter->second);
      }
    }
    
    if (at.empty())
      throw gromos::Exception(argv[0], "please give the atoms you want to perturb.");

    // get the new IACs
    vector<int> iacs;
    {
      Arguments::const_iterator it = args.lower_bound("types"),
          to = args.upper_bound("types");
      for (; it != to; ++it) {
        istringstream is(it->second);
        int type;
        if (!(is >> type))
          throw gromos::Exception(argv[0], "value in @types is not numeric.");
        iacs.push_back(type);
      }
    }

    // get the new charges
    vector<double> charges;
    {
      Arguments::const_iterator it = args.lower_bound("charges"),
          to = args.upper_bound("charges");
      for (; it != to; ++it) {
        istringstream is(it->second);
        double charge;
        if (!(is >> charge))
          throw gromos::Exception(argv[0], "value in @charges is not numeric.");
        charges.push_back(charge);
      }
    }
    
    // get the new masses
    vector<double> masses;
    {
      Arguments::const_iterator it = args.lower_bound("masses"),
          to = args.upper_bound("masses");
      for (; it != to; ++it) {
        istringstream is(it->second);
        double mass;
        if (!(is >> mass))
          throw gromos::Exception(argv[0], "value in @masses is not numeric.");
        masses.push_back(mass);
      }
    }
    
    // this just makes no sense
    if (iacs.empty() && charges.empty() && masses.empty())
      throw gromos::Exception(argv[0], "please give at least one IAC, charge or mass to perturb.");


    cout << "TITLE\n";
    cout << "Perturbation of atoms: ";
    for (int i = 0; i < at.size(); i++) cout << at.toString(i) << " ";
    cout << "From topology: " << args["topo"] << endl;
    cout << "END\n";
    
    cout << "PERTATOMPARAM\n";
    cout << at.size() << endl;
    cout << "#  NR  RES  NAME IAC(A)  MASS(A) CHARGE(A) IAC(B)  MASS(B) CHARGE(B)       ALJ       ACRF\n";
    for (int i = 0; i < at.size(); i++) {
      cout.setf(ios::fixed);
      cout.precision(5);
      cout << setw(5) << at.gromosAtom(i) + 1
          << setw(5) << at.resnum(i) + 1
          << " " << setw(5) << at.name(i)
          << setw(4) << at.iac(i) + 1 << setw(11) << at.mass(i) << setw(11)
          << at.charge(i);
      
      // perturbation in iac?
      int iac;
      if (iacs.empty()) { // no perturbation
        iac = at.iac(i);
      } else {
        // take the right value or the last one
        iac = i < int(iacs.size()) ? iacs[i] : iacs.back();
      }
      cout << setw(4) << iac;
      
      // perturbation in mass?
      double mass;
      if (masses.empty()) { // no perturbation
        mass = at.mass(i);
      } else {
        // take the right value or the last one
        mass = i < int(masses.size()) ? masses[i] : masses.back();
      }
      cout << setw(11) << mass;
      
      // perturbation in charge?
      double charge;
      if (charges.empty()) { // no perturbation
        charge = at.charge(i);
      } else {
        // take the right value or the last one
        charge = i < int(charges.size()) ? charges[i] : charges.back();
      }
      cout << setw(11) << charge;
      
      // alphas
      cout << setw(11) << 0.0 << setw(11) << 0.0 << endl;
    }
    cout << "END\n";

    cout << "PERTPOLPARAM\n   0\nEND\n";
    cout << "PERTATOMPAIR\n   0\nEND\n";
    cout << "PERTBONDSTRETCHH\n   0\nEND\n";
    cout << "PERTBONDSTRETCH\n   0\nEND\n";
    cout << "PERTBONDANGLEH\n   0\nEND\n";
    cout << "PERTBONDANGLE\n   0\nEND\n";
    cout << "PERTIMPROPERDIHH\n   0\nEND\n";
    cout << "PERTIMPROPERDIH\n   0\nEND\n";
    cout << "PERTPROPERDIHH\n   0\nEND\n";
    cout << "PERTPROPERDIH\n   0\nEND\n";
    cout << "PERTCROSSDIH\n   0\nEND\n";
    cout << "PERTCROSSDIHH\n   0\nEND\n";

    return 0;
  } catch (gromos::Exception e) {
    cerr << e.what() << endl;
    return 1;
  }
}
