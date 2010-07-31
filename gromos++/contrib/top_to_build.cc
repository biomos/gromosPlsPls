#include <cassert>
#include <iostream>
#include <iomanip>
#include <map>
#include "../src/args/Arguments.h"
#include "../src/gio/InTopology.h"
#include "../src/gio/InBuildingBlock.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/LJException.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/BuildingBlock.h"
#include "../src/gcore/Exclusion.h"
#include "../src/gcore/BbSolute.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gcore/Bond.h"
#include "../src/gcore/Angle.h"
#include "../src/gcore/Improper.h"
#include "../src/gcore/Dihedral.h"

using namespace gcore;
using namespace gio;
using namespace args;
using namespace std;


int main(int argc, char *argv[]){

  Argument_List knowns;
  knowns << "topo";
  
  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo <topology>\n";
  
  try{
    Arguments args(argc, argv, knowns, usage);


    //create a map to translate a mass back to an index.
    std::map<double, int> mass;
    mass[1.008]=1;
    mass[13.019]=3;
    mass[14.027]=4;
    mass[15.035]=5;
    mass[16.043]=6;
    mass[12.011]=12;
    mass[14.0067]=14;
    mass[15.9994]=16;
    mass[18.9984]=19;
    mass[22.9898]=23;
    mass[24.305]=24;
    mass[28.08]=28;
    mass[30.9738]=31;
    mass[32.06]=32;
    mass[35.453]=35;
    mass[39.948]=39;
    mass[40.08]=40;
    mass[55.847]=56;
    mass[63.546]=63;
    mass[65.37]=65;
    mass[79.904]=80;
    
    InTopology it(args["topo"]);

    System sys(it.system());
    if(sys.numMolecules()!=1){
      throw gromos::Exception("top_to_build", "can only handle topologies with one molecule (use red_top first)");
    }
    
    MoleculeTopology nbb=sys.mol(0).topology();
    

    //Now we need the output, but there is no outBuildingBlock yet!
    cout.precision(5);
    cout << "# Building block converted from\n";
    cout << "# " << args["topo"] << endl;
    cout << "# no previous exclusions taken into account\n";
    
    cout << "#\n";
    cout << "MTBUILDBLSOLUTE" << endl;
    cout << "# RNME" << endl;
    cout << nbb.resName(0) << endl;
    cout << "# NMAT, NLIN" << endl;
    cout << setw(6) << nbb.numAtoms()
	 << setw(6) << 0 << endl;

    cout << "#ATOM ANM  IACM MASS        CGMICGM MAE MSAE" << endl;
    for(int i=0; i<nbb.numAtoms(); i++){
      cout << setw(5) << i+1
	   << setw(5) << nbb.atom(i).name()
	   << setw(5) << nbb.atom(i).iac()+1
	   << setw(5) << mass[nbb.atom(i).mass()]
	   << setw(11) << nbb.atom(i).charge()
	   << setw(4) << nbb.atom(i).chargeGroup();
      
      cout << setw(4) << nbb.atom(i).exclusion().size();
      for(int j=0; j< nbb.atom(i).exclusion().size(); j++)
	cout << setw(5) << nbb.atom(i).exclusion().atom(j)+1;
      cout << endl;
    }
    
    cout << "#  NB" << endl;
    cout << setw(5) << nbb.numBonds() << endl;
    cout << "#  IB   JB  MCB" << endl;
    BondIterator bi(nbb);
    for(;bi;++bi)
      cout << setw(5) << bi()[0]+1
	   << setw(5) << bi()[1]+1
	   << setw(5) << bi().type()+1 << endl;
    cout << "# NBA" << endl;
    cout << setw(5) << nbb.numAngles() << endl;
    cout << "#  IB   JB   KB  MCB" << endl;
    AngleIterator ai(nbb);
    for(;ai;++ai)
      cout << setw(5) << ai()[0]+1
	   << setw(5) << ai()[1]+1
	   << setw(5) << ai()[2]+1
	   << setw(5) << ai().type()+1 << endl;
    cout << "# NIDA" << endl;
    cout << setw(5) << nbb.numImpropers() << endl;
    cout << "#  IB   JB   KB   LB  MCB" << endl;
    ImproperIterator ii(nbb);
    for(;ii;++ii)
      cout <<  setw(5) << ii()[0]+1
	   <<  setw(5) << ii()[1]+1
	   <<  setw(5) << ii()[2]+1
	   <<  setw(5) << ii()[3]+1
	   <<  setw(5) << ii().type()+1 << endl;
    cout << "# NDA" << endl;
    cout << setw(5) << nbb.numDihedrals() << endl;
    cout << "#  IB   JB   KB   LB  MCB" << endl;
    DihedralIterator di(nbb);
    for(;di;++di)
      cout <<  setw(5) << di()[0]+1
	   <<  setw(5) << di()[1]+1
	   <<  setw(5) << di()[2]+1
	   <<  setw(5) << di()[3]+1
	   <<  setw(5) << di().type()+1 << endl;   
    cout << "END" << endl;
    
    

  }
  catch(gromos::Exception e){
    cerr << e.what() << endl;
    return 1;
  }
}





