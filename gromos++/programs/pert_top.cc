/**
 * @file pert_top.cc
 * Creates a perturbation topology to remove interactions for specified atoms
 */

/**
 * @page programs Program Documentation
 *
 * @anchor pert_top
 * @section pert_top Creates a perturbation topology to remove interactions for specified atoms
 * @author @ref co
 * @date 7-6-07
 *
 * Creates a perturbation topology to perturb specified atoms to neutral dummy
 * atoms. A perturbation topology is written that defines a perturbation to
 * change the specified atoms into a specified atom type. The charges of these
 * atoms are set to 0.0. For the first atom, a different value of IACB can be
 * given. This allows the user to change the last atom attached to the atoms
 * that will be disappearing into the appropriate real atoms (e.g. CH2 to CH3).
 * In these cases, the mass of the first perturbed atom will be adapted as well.
 * Use program @ref pt_top to convert the resulting perturbation topology to a
 * different format or to a regular molecular topology.
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;topology&gt; </td></tr>
 * <tr><td> \@atoms</td><td>&lt;@ref AtomSpecifier: atoms to be modified </td></tr>
 * <tr><td> \@types</td><td>&lt;IACB1, IACB of the first and following perturbed atoms&gt; </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  pert_top
    @topo   ex.top
    @atoms  1:34-51
    @types  13 19
 @endverbatim
 *
 * <hr>
 */
#include <cassert>
#include <iostream>
#include <iomanip>
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

  char *knowns[] = {"topo", "atoms", "types"};
  int nknowns = 3;
  
  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo  <topology>\n";
  usage += "\t@atoms <AtomSpecifier: atoms to be modified\n";
  usage += "\t@types <IACB1, IACB of the first and following perturbed atoms>\n";
   
  try{
    Arguments args(argc, argv, nknowns, knowns, usage);

    InTopology it(args["topo"]);
    System sys(it.system());
    utils::AtomSpecifier at(sys);
    {
      Arguments::const_iterator iter=args.lower_bound("atoms");
      
      for(; iter!=args.upper_bound("atoms"); ++iter){
	at.addSpecifier(iter->second);
      }
    }
    int iacb1=0;
    int iacb=0;
    {
      args.check("types",2);
      Arguments::const_iterator iter=args.lower_bound("types");
      iacb1=atoi(iter->second.c_str());
      iter++;
      iacb =atoi(iter->second.c_str());
    }
    
    double cgb=0.0;
    if(at.size()){
      
      cout << "TITLE\n";
      cout << "Perturbation of atoms: ";
      for(int i=0; i< at.size(); i++) cout << at.toString(i) << " ";
      
      if(iacb!=iacb1){
	cout << "\nfirst atom "<< at.toString(0) << " perturbed into IACB = " << iacb1 
	     << " (mass adapted!);\n";
	if(at.size()>1) cout << "subsequent atoms ";
      }
      else
	cout << "\nall atoms ";
      
      
      if(at.size()>1 || iacb==iacb1) {
	cout << "perturbed into IACB = " << iacb << "; all neutral\n";
      }
      
      cout << "From topology: " << args["topo"] << endl;
      cout << "END\n";
      cout << "PERTATOM\n";
      cout << at.size() << endl;
      cout << "# JLA RESNR ATNAM     IACB      WMB      CGB ISCLJ  ISCC\n";
      for(int i=0;i<at.size();i++){
	cout.setf(ios::fixed);
	
	cout << setw(5) << at.gromosAtom(i)+1
	     << setw(6) << at.resnum(i)+1
	     << " " << setw(4) << at.name(i) <<"\t";
	if(i==0&&iacb1!=iacb) {
	  double mass=at.mass(i)+1.008*(iacb1-at.iac(i)-1);
	  cout << setw(2) << iacb1
	       << setw(9) << setprecision(4) << mass;
	}
	else 
	  cout << setw(2) << iacb
	       << setw(9) << setprecision(4) << at.mass(i);
	cout << setw(9) << setprecision(3) << cgb;
	if(at.charge(i)!=0.0)
	  cout << "     1     1\n";
	else
	  cout << "     1     0\n";
      }
      cout << "END\n";
      
      cout << "PERTATOMPAIR\n   0\nEND\n";
      cout << "PERTBONDH\n   0\nEND\n";
      cout << "PERTBOND\n   0\nEND\n";
      cout << "PERTBANGLEH\n   0\nEND\n";
      cout << "PERTBANGLE\n   0\nEND\n";
      cout << "PERTIMPDIHEDRALH\n   0\nEND\n";
      cout << "PERTIMPDIHEDRAL\n   0\nEND\n";
      cout << "PERTDIHEDRALH\n   0\nEND\n";
      cout << "PERTDIHEDRAL\n   0\nEND\n";
    }
    
    return 0;
  }
  catch(gromos::Exception e){
    cerr << e.what() << endl;
    return 1;
  }
}





