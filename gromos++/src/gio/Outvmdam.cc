// gio_Outvmdam.cc

#include <cassert>
#include <iostream>
#include <iomanip>
#include <set>
#include <sstream>
#include "Outvmdam.h"
#include "../gcore/System.h"
#include "../gcore/Molecule.h"
#include "../gcore/LJException.h"
#include "../gcore/MoleculeTopology.h"
#include "../gcore/Solvent.h"
#include "../gcore/SolventTopology.h"
#include "../gcore/AtomTopology.h"
#include "../gmath/Vec.h"
#include "../gcore/Box.h"
#include "../utils/AtomSpecifier.h"

using gio::Outvmdam;
using namespace gcore;
using namespace std;
using namespace utils;

class gio::Outvmdam_i {
  friend class gio::Outvmdam;
  ostream &d_os;
  int d_count, d_resoff, d_switch;
  int count, eq;
  double d_factor;

  Outvmdam_i(ostream &os, double factor = 10.0) :
  d_os(os), d_count(0), d_switch(0), d_factor(10.0) {
    count = 0;
  }

  ~Outvmdam_i() {
  }
  void select(const string &thing);
  void writeSingleS(const Solvent &sol);
  void writeSingleM(const Molecule &mol);
  void writeAtomSpecifier(const AtomSpecifier & atoms);
};

Outvmdam::Outvmdam(ostream &os, double f) :
OutCoordinates(),
d_this(new Outvmdam_i(os, f)), factor(f) {
}

Outvmdam::Outvmdam(double factor) :
OutCoordinates(), factor(factor) {
  d_this = 0;
}

Outvmdam::~Outvmdam() {
  if (d_this)delete d_this;
}

void Outvmdam::writeTitle(const string &title) {
  d_this->d_os << "TITLE       " << title << "\n";
}

void Outvmdam::writeTimestep(const int step, const double time) {
  // not implemented
}

void Outvmdam::select(const string &thing) {
  if (thing == "ALL") {
    d_this->d_switch = 1;
  } else if (thing == "SOLVENT") {
    d_this->d_switch = 2;
  } else {
    d_this->d_switch = 0;
  }
}

void Outvmdam::open(ostream &os) {
  if (d_this) {
    delete d_this;
  }
  d_this = new Outvmdam_i(os, factor);
}

void Outvmdam::close() {
  if (d_this)delete d_this;
  d_this = 0;
}

Outvmdam &Outvmdam::operator<<(const gcore::System &sys) {
  d_this->d_count = 0;
  d_this->d_resoff = 1;
  d_this->eq = 10;
  d_this->count = 0;
  if (d_this->d_switch <= 1)
    for (int i = 0; i < sys.numMolecules(); ++i)
      d_this->writeSingleM(sys.mol(i));
  if (d_this->d_switch >= 1)
    for (int i = 0; i < sys.numSolvents(); ++i)
      d_this->writeSingleS(sys.sol(i));

  d_this->d_os << "\n";

  return *this;
}

Outvmdam & Outvmdam::operator<<(const AtomSpecifier & atoms) {
  d_this->d_count = 0;
  d_this->d_resoff = 1;
  d_this->eq = 10;
  d_this->count = 0;
  d_this->writeAtomSpecifier(atoms);
  d_this->d_os << "\n";


  return *this;
}


void gio::Outvmdam_i::writeSingleM(const Molecule &mol) {
  d_os.setf(ios::fixed, ios::floatfield);
  d_os.setf(ios::unitbuf);
  d_os.precision(3);

  for (int i = 0; i < mol.numAtoms(); ++i) {
    ++d_count;


    d_os.setf(ios::right, ios::adjustfield);

    if (count && count % eq == 0) d_os << endl;
    d_os << setw(8) << mol.pos(i)[0]*d_factor;
    count++;
    if (count && count % eq == 0) d_os << endl;
    d_os << setw(8) << mol.pos(i)[1]*d_factor;
    count++;
    if (count && count % eq == 0) d_os << endl;
    d_os << setw(8) << mol.pos(i)[2]*d_factor;
    count++;

  }
}

void gio::Outvmdam_i::writeSingleS(const Solvent &sol) {
  d_os.setf(ios::fixed, ios::floatfield);
  d_os.setf(ios::unitbuf);
  d_os.precision(3);

  for (int i = 0; i < sol.numPos(); ++i) {
    ++d_count;


    d_os.setf(ios::right, ios::adjustfield);

    if (count && count % eq == 0) d_os << endl;
    d_os << setw(8) << sol.pos(i)[0]*d_factor;
    count++;
    if (count && count % eq == 0) d_os << endl;
    d_os << setw(8) << sol.pos(i)[1]*d_factor;
    count++;
    if (count && count % eq == 0) d_os << endl;
    d_os << setw(8) << sol.pos(i)[2]*d_factor;
    count++;

  }
}

void gio::Outvmdam_i::writeAtomSpecifier(const AtomSpecifier& atoms) {
  d_os.setf(ios::fixed, ios::floatfield);
  d_os.setf(ios::unitbuf);
  d_os.precision(3);

  for (int i = 0; i < atoms.size(); ++i) {
    ++d_count;

    d_os.setf(ios::right, ios::adjustfield);
    if (count && count % eq == 0) d_os << endl;
    d_os << setw(8) << atoms.pos(i)[0]*d_factor;
    count++;
    if (count && count % eq == 0) d_os << endl;
    d_os << setw(8) << atoms.pos(i)[1]*d_factor;
    count++;
    if (count && count % eq == 0) d_os << endl;
    d_os << setw(8) << atoms.pos(i)[2]*d_factor;
    count++;
  }
}
