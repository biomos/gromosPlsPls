// args_Arguments.t.cc

#include "Arguments.h"
#include <iostream>

using namespace std;
using args::Arguments;
//using namespace std;

namespace std{
  std::ostream &operator<<(std::ostream &os, const Arguments &args){
    for(Arguments::const_iterator iter=args.begin();
	iter!=args.end();++iter){
      os << iter->first << ' ' << iter->second << endl;
    }
    return os;
  }
}

int main(int argc, char **argv){

  char *knowns[] = {"bla", "gug"};
  string usage = argv[0];
  usage+=" @bla <testargs> @gug <testargs>";
  
 try{
 Arguments args(argc, argv, 2, knowns, usage) ;
 cout << args.count("gug")  << endl;
 
 cout << args;


 }
 catch(const gromos::Exception &e){
   cerr << e.what() << endl;
   exit(1);
 }

  return 0;
}
