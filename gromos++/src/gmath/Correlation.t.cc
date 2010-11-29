#include "Correlation.h"
#include "../gromos/Exception.h"
#include <cstdlib>
#include <iostream>
#include <vector>
#include <string>

using gmath::Correlation;

using namespace std;

int main() {
  try {

    vector<double> a;
    vector<double> b;
    srand(23123);

    for (int j = 0; j < 1000; j++) {
      a.push_back(rand() * 1000.0 / RAND_MAX);
      b.push_back(j);
    }
    Correlation c(a, b);
    Correlation d(a, b);
    c.calc_direct();
    //for(int i=0; i< c.size(); i++){
    //	cout << i << "\t" << c[i] << endl;
    //  }
    d.calc_fft();

    for (int i = 0; i < c.size(); i++) {
      cout << i << "\t" << c[i] << "\t" << d[i] << endl;
    }

    vector<double> w(1000);
    vector<double> s(1000);
    d.spectrum(w, s, 0.002, 0.5);
    for (size_t i = 0; i < w.size(); i++)
      cout << w[i] << "\t" << s[i] << endl;

    return 0;
  }  catch (gromos::Exception e) {
    cerr << e.what() << endl;
    return 1;
  }
}
