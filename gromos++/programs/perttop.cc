#include <iostream>
#include <iomanip>
#include "../src/args/Arguments.h"
#include "../src/gio/InTopology.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/System.h"

using namespace gcore;
using namespace gio;
using namespace args;


int main(int argc, char *argv[]){

  char *knowns[] = {"topo", "atoms", "first"};
  int nknowns = 3;
  
  string usage = argv[0];
  usage += "\n\t@topo <topology>\n";
  usage += "\t@atoms <first and last atom in new topology/all>\n";
  usage += "\t@first <IACB of the first perturbed atom>\n";
   
  try{
    int atoff, atend;
    
    Arguments args(argc, argv, nknowns, knowns, usage);

    args.check("atoms",1);
    InTopology it(args["topo"]);
    System sys(it.system());
    MoleculeTopology mti=sys.mol(0).topology();
    AtomTopology ato;
    if (args["atoms"]=="all") {
      atoff=0;
      atend=mti.numAtoms();
    } else{
      args.check("atoms",2);
      Arguments::const_iterator iter=args.lower_bound("atoms");
      atoff=atoi(iter++->second.c_str())-1;
      atend=atoi(iter->second.c_str());
      if(atoff<0||atoff>mti.numAtoms())
        throw Arguments::Exception("illegal first atom value");
      if(atend<atoff||atend>mti.numAtoms())
        throw Arguments::Exception("illegal last atom value");
    }
    args.check("first",1);
    Arguments::const_iterator iter=args.lower_bound("first");
    int iacb1=atoi(iter->second.c_str());
    
    int offset=1;
    
    int iacb=19;
    double cgb=0.0;
    
    cout << "TITLE\n";
    cout << "Perturbation of "<< mti.resName(0) <<": ";
    
    if(iacb!=iacb1)
      cout << "atom "<< atoff+offset++ << " into IACB = " << iacb1 << ";\n";
    if(atoff+offset<=atend) {
       cout << "atom " << atoff+offset;
      if(atoff+offset<atend) cout << " to " << atend;
      cout << " into IACB = " << iacb << "; all neutral\n";
    }
    
    cout << "From topology: " << args["topo"] << endl;
    cout << "END\n";
    cout << "PERTATOM\n";
    cout << atend-atoff << endl;
    cout << "# JLA RESNR ATNAM     IACB      WMB      CGB ISCLJ  ISCC\n";
    for(int i=atoff;i<atend;i++){
      cout.setf(ios::fixed);
      
      cout << setw(5) << i+1
           << setw(6) << mti.resNum(i)+1
           << " " << setw(4) << mti.atom(i).name() <<"\t";
      if(i==atoff&&iacb1!=iacb) {
        double mass=mti.atom(i).mass()+1.008*(iacb1-mti.atom(i).iac()-1);
	cout << setw(2) << iacb1
             << setw(9) << setprecision(4) << mass;
      }
      else 
        cout << setw(2) << iacb
             << setw(9) << setprecision(4) << mti.atom(i).mass();
      cout << setw(9) << setprecision(3) << cgb;
      if(mti.atom(i).charge()!=0.0)
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
    
    return 0;
  }
  catch(gromos::Exception e){
    cerr << e.what() << endl;
    return 1;
  }
}





