#include "Exception.h"
#include <cstdlib>
#include <iostream>

using namespace std;

struct my_Exception : public gromos::Exception {

  my_Exception(const string & str) :
  gromos::Exception("my_Class", str) {
  }
};

int main() {
  my_Exception e("SOME ERROR!");
  cout << e.what() << endl;
  try {
    throw;
  } catch (const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}

/* OUTPUT:
my_Class: SOME ERROR!
 */
