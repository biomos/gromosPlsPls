// gcore_Exclusion.t.cc

#include "Exclusion.h"
#include <iostream>

using gcore::Exclusion;
using namespace std;

static void output(const Exclusion &e) {
  cout << "size: " << e.size() << endl
          << "Exclusions: ";

  for (int i = 0; i < e.size(); ++i)
    cout << e.atom(i) << ' ';
  cout << endl;
}

int main() {
  Exclusion e;
  e.insert(3);
  e.insert(2);
  e.insert(5);
  e.insert(4);
  e.insert(4);

  output(e);

  Exclusion g = e;
  output(g);

  Exclusion f;
  f = g;
  output(f);
  f.erase(3);
  output(f);
  return 0;
}


/* Output:
size: 4
Exclusions: 2 3 4 5 
size: 4
Exclusions: 2 3 4 5 
size: 4
Exclusions: 2 3 4 5 
size: 3
Exclusions: 2 4 5 
 */
