#include "VirtualAtom.cc"
#include "../gcore/System.h"
#include "../gio/InTopology.h"
#include "../gio/InG96.h"
#include <string>
#include <iostream>

using namespace gcore;
using namespace gio;
using namespace utils;

ostream &operator<<(ostream &os, const gmath::Vec &v)
{
  os << v[0] << ' ' << v[1] << ' ' << v[2];
  return os;
}


int main(int argc, char *argv[]){
  if(argc != 6){
    cerr << "Usage: " << argv[0] << " <topology> <coordinates> <mol nr.> <atom nr.> <type>" << endl;
    exit(1);
  }
  string top = argv[1];
  string coord = argv[2];
  
  int mol = atoi(argv[3])-1;
  int atom = atoi(argv[4])-1;
  int type = atoi(argv[5]);
  try{
    
    InTopology it(top);
    System sys=it.system();
    
    InG96 ic(coord);
    ic >> sys;
    
    
    VirtualAtom virt(sys,mol,atom,type);
    
    cout << "Virtual atom created with configuration\n";
    for(int i=0;i<4;++i)
      cout << virt[i] << ' ';
    cout << virt.type() << endl;
    
    cout << "Position calculated as: ";
    cout << virt.pos() << endl;
    
  }
  catch(gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  
  return 0;
}
