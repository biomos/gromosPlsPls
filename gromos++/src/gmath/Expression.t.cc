#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <cstdlib>
#include "Expression.h"
#include "../gromos/Exception.h"

using gmath::Expression;

using namespace std;

int main(int argc, char *argv[]) {
  try {
    if (argc < 2) {
      cerr << "Usage: " + string(argv[0]) + " <expression>" << endl;
      exit(1);
    }
    ostringstream os;
    for (int i = 1; i < argc; i++)
      os << argv[i] << " ";

    string s = os.str();
    Expression e(s);
    e.writeExpression(cout);
    cout << e.value() << endl;
  }  catch (gromos::Exception e) {
    cerr << e.what() << endl;
    return 1;
  }
}
