#include <iostream>
#include "../gio/InTopology.h"
#include "PropertyContainer.h"
#include "../gcore/System.h"
#include "../gcore/MoleculeTopology.h"
#include "../gcore/Molecule.h"

using namespace gcore;
using namespace gio;
using namespace utils;

int main(int argc, char *argv[]){
  if(argc !=3){
    cerr << "Usage: " + string(argv[0]) + " <Topology> <propertyspecifier>\n";
    exit(1);
  }
  try{
    InTopology it(argv[1]);
    System sys(it.system());
    string s=argv[2];
    
    // PropertyContainer as(sys, s);
    PropertyContainer bs(sys);    
    bs.addSpecifier(argv[2]);

    cout << bs << endl;
    
    return 0;
  }
  catch(gromos::Exception e){
    cerr << e.what() << endl;
    return 1;
  }
}





