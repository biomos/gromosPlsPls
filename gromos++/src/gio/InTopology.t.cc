#include <iostream>
#include "InTopology.h"
#include "OutTopology.h"
#include "../gcore/GromosForceField.h"
#include "../gcore/System.h"

using namespace gcore;
using namespace gio;

int main(int argc, char *argv[]){
  if(argc !=2){
    cerr << "Usage: " + string(argv[0]) + " <Topology>\n";
    exit(1);
  }
  try{
    InTopology it(argv[1]);
  System sys(it.system());
  GromosForceField gff(it.forceField());
  cout << sys.numMolecules();
  OutTopology ot(cout);
  ot.setTitle(it.title());
  ot.write(sys,gff);

  return 0;
  }
  catch(gromos::Exception e){
    cerr << e.what() << endl;
    return 1;
  }
}





