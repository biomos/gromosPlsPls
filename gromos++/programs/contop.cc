#include <iostream>
#include "../src/args/Arguments.h"
#include "../src/gio/InTopology.h"
#include "../src/gio/InParameter.h"
#include "../src/gio/OutTopology.h"
#include "../src/gcore/GromosForceField.h"
#include "../src/gcore/System.h"

using namespace gcore;
using namespace gio;
using namespace args;


int main(int argc, char *argv[]){

  char *knowns[] = {"topo", "param"};
  int nknowns = 2;
  
  string usage = argv[0];
  usage += "\n\t@topo <topology>\n";
  usage += "\t@param <gromos parameter file/topo>\n";
  
  try{
    Arguments args(argc, argv, nknowns, knowns, usage);

    InTopology it(args["topo"]);


    System sys(it.system());

    
    OutTopology ot(cout);
    string addtitle;
    addtitle+="\nCONTOP parameters: "+args["param"];
    ot.setTitle(it.title()+addtitle);

    InParameter ip(args["param"]);
    GromosForceField gff(ip.forceField());
    gff.setFpepsi(it.forceField().fpepsi());
    gff.setHbar(it.forceField().hbar());
     
    ot.write(it.system(),gff);
    
  }
  catch(gromos::Exception e){
    cerr << e.what() << endl;
    return 1;
  }
}





