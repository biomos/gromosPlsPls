#include <iostream>
#include "InBuildingBlock.h"
#include "../gcore/BuildingBlock.h"
#include "../gcore/BbSolute.h"
#include "../gcore/SolventTopology.h"

using namespace gcore;
using namespace gio;

int main(int argc, char *argv[]){
  if(argc !=2){
    cerr << "Usage: " + string(argv[0]) + " <Building Block file>\n";
    exit(1);
  }
  try{
    InBuildingBlock ibb(argv[1]);
    BuildingBlock bld(ibb.building());

    cout << "Number of BbSolutes: " << bld.numBbSolutes() << endl;
    cout << " BbSolute number 24: " << bld.bb(23).resName() << endl;
    cout << "  has " << bld.bb(23).numAtoms() << " atoms and "
         << bld.bb(23).numPexcl() << " preceding exclusions" << endl;
    cout << "Number of BbSolvents: " << bld.numBbSolvents() << endl;
    cout << "  BbSolvent number 3: " << bld.bs(2).solvName() << endl;
    cout << "  has " << bld.bs(2).numAtoms() << " atoms" << endl;
    cout << "Number of BbEnds: " << bld.numBbEnds() << endl;
    cout << "  BbEnd number 2: " << bld.be(2).resName() << endl;
    cout << "  will replace " << bld.be(2).rep() << " atoms" << endl;

    int index=bld.findBb("DADE");
    BbDihIt di(bld.bb(index-1));
    int count =0;
    for(;di;++di) count++; 
    cout << "DADE has " << count << endl;
  return 0;
  }
  catch(gromos::Exception e){
    cerr << e.what() << endl;
    return 1;
  }
}





