// args_Arguments.t.cc

#include <cassert>
#include <cstdlib>
#include "Arguments.h"
#include <iostream>

using namespace std;
using namespace args;

int debug_level = 0;

namespace std {

  std::ostream & operator<<(std::ostream &os, const Arguments &args) {
    for (Arguments::const_iterator iter = args.begin();
            iter != args.end(); ++iter) {
      os << iter->first << ' ' << iter->second << endl;
    }
    return os;
  }
}

int main(int argc, char **argv) {
  Argument_List knowns;
  knowns << "bla" << "gug";
  string usage = argv[0];
  usage += " @bla <testargs> @gug <testargs>";

  try {
    Arguments args(argc, argv, knowns, usage);
    cout << args.count("gug") << endl;

    cout << args;


  } catch (const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }

  return 0;
}
