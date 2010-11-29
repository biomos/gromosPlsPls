// utils_Noe.t.cc

#include <cassert>
#include <cstdlib>
#include "Noe.h"
#include "../gcore/System.h"
#include "../gio/InTopology.h"
#include "../gio/InG96.h"
#include <string>
#include <iostream>

using namespace utils;
using namespace gcore;
using namespace gio;

using namespace std;

int main(int argc, char *argv[]) {
  if (argc != 4) {
    cerr << "Usage: " << argv[0] << " <topology> <coordinates> <Noe line>" << endl;
    exit(1);
  }
  string top = argv[1];
  string coord = argv[2];

  string line = argv[3];

  try {

    InTopology it(top);
    System sys = it.system();

    InG96 ic(coord);
    ic >> sys;


    Noe noe(sys, line, 0.1, 0.15);

    cout << "The distance restraints are \n";
    for (int i = 0; i < noe.numDistances(); ++i)
      cout << noe.distRes(i) << endl;

    cout << "\nThe distances calculated are:\n";
    for (int i = 0; i < noe.numDistances(); ++i)
      cout << noe.distance(i) << ' ';

    cout << "\n\nThe following reference distances were found. \nUncorrected: ";
    for (int i = 0; i < noe.numReferences(); ++i)
      cout << noe.reference(i) << ' ';
    cout << "\nCorrected: ";
    for (int i = 0; i < noe.numReferences(); ++i)
      cout << noe.correctedReference(i) << ' ';
    cout << endl;


  } catch (gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }

  return 0;
}
