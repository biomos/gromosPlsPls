#include <cassert>
#include <iostream>
#include <iomanip>
#include <set>
#include "../src/args/Arguments.h"
#include "../src/gio/InBuildingBlock.h"
#include "../src/gcore/BuildingBlock.h"
#include "../src/gcore/BbSolute.h"
#include "../src/gcore/SolventTopology.h"


#include "../src/gcore/Bond.h"
#include "../src/gcore/Angle.h"
#include "../src/gcore/Dihedral.h"
#include "../src/gcore/Improper.h"
#include "../src/gcore/Exclusion.h"
#include "../src/gcore/AtomTopology.h"

#include "../src/gio/Ginstream.h"

using namespace gcore;
using namespace gio;
using namespace args;
using namespace std;


int main(int argc, char *argv[]){

  char *knowns[] = {"build"};
  int nknowns = 1;
  
  string usage = argv[0];
  usage += "\n\t@build <mtb-file>\n";

  try{
    Arguments args(argc, argv, nknowns, knowns, usage);

    // read in the building block file
    InBuildingBlock ibb(args["build"]);
    BuildingBlock mtb(ibb.building());

    // keep track of the number of atom, bond, angle etc. types required
    int atomtype=0, bondtype=0, angletype=0, dihedraltype=0, impropertype=0;
    string atomname, bondname, anglename, dihedralname, impropername;
    set<int> atoms, bonds, angles, dihedrals, impropers;

    // print out a list of the buildingblocks

    cout << "Number of Solute Building Blocks (MTBUILDBLSOLUTE): "
	 << mtb.numBbSolutes() << endl;
    for(int i=0; i<mtb.numBbSolutes(); i++){
      cout << "  " << mtb.bb(i).resName() << endl;
      for(int j=0; j<mtb.bb(i).numAtoms(); j++){
	atoms.insert(mtb.bb(i).atom(j).iac());
	if(mtb.bb(i).atom(j).iac() > atomtype) {
	  atomname=mtb.bb(i).resName();
	  atomtype=mtb.bb(i).atom(j).iac();
	}
      }
      
      BondIterator bi(mtb.bb(i));
      for(;bi;++bi){
	bonds.insert(bi().type());
	if(bi().type()>bondtype){
	  bondname= mtb.bb(i).resName();
	  bondtype=bi().type();
	}
      }
      
      AngleIterator ai(mtb.bb(i));
      for(; ai; ++ai){
	angles.insert(ai().type());
	if(ai().type()>angletype){
	  anglename=mtb.bb(i).resName();
	  angletype=ai().type();
	}
      }
      
      ImproperIterator ii(mtb.bb(i));
      for(; ii; ++ii){
	impropers.insert(ii().type());
	if(ii().type()>impropertype){
	  impropername=mtb.bb(i).resName();
	  impropertype=ii().type();
	}
      }
      
      DihedralIterator di(mtb.bb(i));
      for(; di; ++di){
	dihedrals.insert(di().type());
	if(di().type()>dihedraltype){
	  dihedralname=mtb.bb(i).resName();
	  dihedraltype=di().type();
	}
      }
    }
    cout << endl;
    cout << "Number of Solvent Building Blocks (MTBUILDBLSOLVENT): "
	 << mtb.numBbSolvents() << endl;
    for(int i=0; i<mtb.numBbSolvents(); i++){
      cout << "  " << mtb.bs(i).solvName() << endl;
      for(int j=0; j<mtb.bs(i).numAtoms(); j++){
	atoms.insert(mtb.bs(i).atom(j).iac());
	if(mtb.bs(i).atom(j).iac() > atomtype){
	  atomname=mtb.bs(i).solvName();
	  atomtype=mtb.bs(i).atom(j).iac();
	}
      }
    }
    cout << endl;
    
    
    cout << "Number of End-group Building Blocks (MTBUILDBLEND): "
	 << mtb.numBbEnds() << endl;
    for(int i=0; i<mtb.numBbEnds(); i++){
	
      cout << "  " << mtb.be(i).resName() << endl;
      for(int j=0; j<mtb.be(i).numAtoms(); j++){
	
	atoms.insert(mtb.be(i).atom(j).iac());
	if(mtb.be(i).atom(j).iac() > atomtype) {
	  atomname=mtb.be(i).resName();
	  atomtype=mtb.be(i).atom(j).iac();
	}
      }
      
      BondIterator bi(mtb.be(i));
      for(;bi;++bi){
	bonds.insert(bi().type());
	if(bi().type()>bondtype){
	  bondname=mtb.be(i).resName();
	  bondtype=bi().type();
	}
      }
      
      AngleIterator ai(mtb.be(i));
      for(; ai; ++ai){
	angles.insert(ai().type());
	if(ai().type()>angletype){
	  anglename=mtb.be(i).resName();
	  angletype=ai().type();
	}
      }
      
      ImproperIterator ii(mtb.be(i));
      for(; ii; ++ii){
	impropers.insert(ii().type());
	if(ii().type()>impropertype){
	  impropername=mtb.be(i).resName();
	  impropertype=ii().type();
	}
      }
      
      DihedralIterator di(mtb.be(i));
      for(; di; ++di){
	dihedrals.insert(di().type());
	if(di().type()>dihedraltype){
	  dihedralname=mtb.be(i).resName();
	  dihedraltype=di().type();
	}
      }
    }
    cout << endl;
    cout << "Highest types encountered for:   (in building block)" << endl;
    cout << "  atoms     : " << setw(5) << atomtype+1 
	 << "  (" << atomname << ")" << endl;
    cout << "  bonds     : " << setw(5) << bondtype +1 
	 << "  (" << bondname << ")" << endl;
    cout << "  angles    : " << setw(5) << angletype +1 
	 << "  (" << anglename << ")" << endl;
    cout << "  impropers : " << setw(5) << impropertype +1 
	 << "  (" << impropername << ")"<< endl;
    cout << "  dihedrals : " << setw(5) << dihedraltype +1 
	 << "  (" << dihedralname << ")" << endl;
    cout << endl;
    cout << "Types that were not seen below these values: " << endl;
    cout << "  atoms     : ";
    for(int i=0; i<atomtype; i++)
      if(atoms.count(i)==0) cout << setw(5) << i+1 << " ";
    cout << endl;
    
    cout << "  bonds     : ";
    for(int i=0; i<bondtype; i++)
      if(bonds.count(i)==0) cout << setw(5) <<i+1 << " ";
    cout << endl;
    
    cout << "  angles    : ";
    for(int i=0; i<angletype; i++)
      if(angles.count(i)==0) cout << setw(5) <<i+1 << " ";
    cout << endl;
    
    cout << "  impropers : ";
    for(int i=0; i<impropertype; i++)
      if(impropers.count(i)==0) cout << setw(5) <<i+1 << " ";
    cout << endl;
    
    cout << "  dihedrals : ";
    for(int i=0; i<dihedraltype; i++)
      if(dihedrals.count(i)==0) cout << setw(5) <<i+1 << " ";
    cout << endl;
    
  }
  catch(gromos::Exception e){
    cerr << e.what() << endl;
    return 1;
  }
}




