#include <vector>

#include "Vec.h"
#include "../gcore/Box.h"
#include "Mesh.h"
#include <iostream>

using namespace std;
using namespace gcore;
using namespace gmath;

int main(int argc, char** argv) {
  Box b(Vec(10.0,  0.0,  0.0),
        Vec( 0.0, 10.0,  0.0),
        Vec( 2.0,  2.0, 10.0));

  Vec centre = b.K() + b.L() + b.M();
  centre *= 0.5;

  Mesh<double> m;
  m.setBox(b);
  m.setTitle("Nathans Mesh");
  m.resize(100, 60, 40);
  m = 0.0;

  for (int i = 0; i < m.size()[0]; ++i) {
    for (int j = 0; j < m.size()[1]; ++j) {
      for (int k = 0; k < m.size()[2]; ++k) {
        m(i, j, k) = (centre - m.pos(i, j, k)).abs();
      }
    }
  }

  m.write(cout);
  return 0;
}
