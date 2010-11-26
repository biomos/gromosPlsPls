#include "Distribution.h"
#include <iostream>

using gmath::Distribution;

using namespace std;

int main() {
  try {

    Distribution di(0, 10, 5);
    cout << "distribution\n";
    di.write(cout);
    di.add(3.2);
    di.add(0.4);
    di.add(4.0);
    di.add(8.0);
    double test = 12;
    if (test != di.add(test)) cout << "value " << test
            << " out of range, not added\n";

    cout << "\ndistribution after adding elements\n";

    di.write(cout);
    cout << "average " << di.ave() << endl;
    cout << "rmsd " << di.rmsd() << endl;
    cout << "number of elements in section 4: " << di[3] << endl;
    cout << "middle value of section 4: " << di.value(3) << endl;

    return 0;
  }  catch (gromos::Exception e) {
    cerr << e.what() << endl;
    return 1;
  }
}
