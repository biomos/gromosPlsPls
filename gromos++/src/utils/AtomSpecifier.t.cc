#include <iostream>
#include "../gio/InTopology.h"
#include "AtomSpecifier.h"
#include "../gcore/System.h"
#include "../gcore/MoleculeTopology.h"
#include "../gcore/Molecule.h"

using namespace gcore;
using namespace gio;
using namespace utils;

int main(int argc, char *argv[]){
  if(argc !=3){
    cerr << "Usage: " + string(argv[0]) + " <Topology> <atomspecifier>\n";
    exit(1);
  }
  try{
    InTopology it(argv[1]);
    System sys(it.system());
    string s=argv[2];
    
    AtomSpecifier as(sys, s);
    AtomSpecifier bs(sys);
    
    cout << s << " consists of " << as.size() << " atoms:\n";
    for(int i=0; i< as.size();i++)
      cout << as.mol(i) << " : " << as.atom(i) << endl;
    
    
    bs.addSpecifier("2:1-3");
    
    AtomSpecifier cs(sys);
    cs.addAtom(1,4);
    
    cs=as+bs;
    cout << endl;
    
    cout << endl <<  cs.size() << endl;
    for(int i=0; i< cs.size();i++)
      cout << cs.mol(i) << " : " << cs.atom(i) << endl;
    cs.removeAtom(1,1);
    
    cout << endl <<  cs.size() << endl;
    for(int i=0; i< cs.size();i++)
      cout << cs.mol(i) << " : " << cs.atom(i) << endl;

    cs.addAtom(1,6);
    cs.addType("CA");
    
    cout << endl <<  cs.size() << endl;
    for(int i=0; i< cs.size();i++)
      cout << cs.mol(i) << " : " << cs.atom(i) << endl;

    

    
  return 0;
  }
  catch(gromos::Exception e){
    cerr << e.what() << endl;
    return 1;
  }
}





