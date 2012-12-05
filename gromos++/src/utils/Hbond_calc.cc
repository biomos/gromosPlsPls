#include <cassert>
#include "../args/Arguments.h"
#include "../bound/Boundary.h"
#include "../args/BoundaryParser.h"
#include "../gio/InG96.h"
#include "../gio/Ginstream.h"
#include "Neighbours.h"

#include "Hbond_calc.h"
#include "Hbond.h"

using namespace bound;
using namespace gio;
using namespace args;

using utils::HB_calc;
using args::Arguments;

void HB_calc::readinmasses(std::string filename) {
  //new read in stuff, using the Ginstream...
  Ginstream nf(filename);
  vector<string> buffer;
  nf.getblock(buffer);

  if (buffer[0] != "HYDROGENMASS")
    throw gromos::Exception("Hbondcalc", "Mass file does not contain a HYDROGENMASS block!");

  istringstream is;
  double mass = 0;

  //read in the hydrogen masses...
  for (unsigned int j = 1; j < buffer.size() - 1; j++) {
    is.clear();
    is.str(buffer[j]);
    is >> mass;
    mass_hydrogens.push_back(mass);
  }

  //get the ACCEPTORMASS block
  nf.getblock(buffer);

  if (buffer[0] != "ACCEPTORMASS")
    throw gromos::Exception("Hbondcalc", "Mass file does not contain a ACCEPTORMASS block!");

  //read in the acceptor masses...
  for (unsigned int j = 1; j < buffer.size() - 1; j++) {
    is.clear();
    is.str(buffer[j]);
    is >> mass;
    mass_acceptors.push_back(mass);
  }
}//end HB_calc::readinmasses()

void HB_calc::determineAtoms() {
  // we allow for the definition of four groups of atoms
  // 1. Donor atoms of group A
  // 2. Acceptor atoms of group A
  // 3. Donor atoms of group B
  // 4. Acceptor atoms of group B

  // everything can contain solvent. So to find out how many of those we have,
  // read a frame
  readframe();

  Arguments::const_iterator iter = args -> lower_bound("DonorAtomsA");
  Arguments::const_iterator to = args -> upper_bound("DonorAtomsA");

  for (; iter != to; iter++) {
    string spec = iter->second.c_str();
    donors.addSpecifier(spec);
  }

  iter = args -> lower_bound("AcceptorAtomsA");
  to = args -> upper_bound("AcceptorAtomsA");

  for (; iter != to; iter++) {
    string spec = iter->second.c_str();
    acceptors.addSpecifier(spec);
  }

  //sort them and find the atoms bound to the donor
  donors.sort();
  acceptors.sort();

  int m, a;
  for (int i = 0; i < donors.size(); ++i) {
    m = donors.mol(i);
    a = donors.atom(i);
    if (m < 0) {
      int j = a % sys->sol(0).topology().numAtoms();
      Neighbours neigh(*sys, 0, j, 0);
      bound.addAtomStrict(-1, a - j + neigh[0]);
    } else {
      Neighbours neigh(*sys, m, a);
      bound.addAtomStrict(m, neigh[0]);
    }
  }
  // store how many acceptors and donor we have in A
  num_A_donors = donors.size();
  num_A_acceptors = acceptors.size();

  // if there is no B specified, we take A
  if (args->count("DonorAtomsB") <= 0 && args->count("AcceptorAtomsB") <= 0) {
    for (int i = 0; i < num_A_donors; i++) {
      donors.addAtomStrict(donors.mol(i), donors.atom(i));
      bound.addAtomStrict(bound.mol(i), bound.atom(i));
    }
    for (int i = 0; i < num_A_acceptors; i++) {
      acceptors.addAtomStrict(acceptors.mol(i), acceptors.atom(i));
    }
  } else {
    AtomSpecifier donor_B(*sys);
    AtomSpecifier bound_B(*sys);
    AtomSpecifier acceptor_B(*sys);

    iter = args -> lower_bound("DonorAtomsB");
    to = args -> upper_bound("DonorAtomsB");
    for (; iter != to; iter++) {
      string spec = iter->second.c_str();
      donor_B.addSpecifier(spec);
    }
    iter = args -> lower_bound("AcceptorAtomsB");
    to = args -> upper_bound("AcceptorAtomsB");
    for (; iter != to; iter++) {
      string spec = iter->second.c_str();
      acceptor_B.addSpecifier(spec);
    }



    // and sort these as well and find the hydrogens bound to donor_B
    donor_B.sort();
    acceptor_B.sort();
    for (int i = 0; i < donor_B.size(); ++i) {
      m = donor_B.mol(i);
      a = donor_B.atom(i);
      if (m < 0) {
        int j = a % sys->sol(0).topology().numAtoms();
        Neighbours neigh(*sys, 0, j, 0);
        bound_B.addAtomStrict(-1, a - j + neigh[0]);
      } else {
        Neighbours neigh(*sys, m, a);
        bound_B.addAtomStrict(m, neigh[0]);
      }
    }

    // copy them into the d_donors, d_bound, d_acceptors
    for (int i = 0; i < donor_B.size(); i++) {
      donors.addAtomStrict(donor_B.mol(i), donor_B.atom(i));
      bound.addAtomStrict(bound_B.mol(i), bound_B.atom(i));
    }
    for (int i = 0; i < acceptor_B.size(); i++) {
      acceptors.addAtomStrict(acceptor_B.mol(i), acceptor_B.atom(i));
    }

  }
} //end HB_calc::determineAtoms()

void HB_calc::determineAtomsbymass() {
  bool keep = false;

  //donors
  for (int i = 0; i < donors.size(); ++i) {
    keep = false;
    for (unsigned int j = 0; j < mass_hydrogens.size(); ++j) {
      if (donors.mass(i) == mass_hydrogens[j]) {
        keep = true;
        break;
      }
    }
    if (!keep) {
      donors.removeAtom(i);
      bound.removeAtom(i);
      if (i < num_A_donors) num_A_donors--;
      i--;
    }
  }

  // acceptors 
  for (int i = 0; i < acceptors.size(); ++i) {

    keep = false;
    for (unsigned int j = 0; j < mass_acceptors.size(); ++j) {
      if (acceptors.mass(i) == mass_acceptors[j]) {
        keep = true;
        break;
      }
    }
    if (!keep) {
      acceptors.removeAtom(i);
      if (i < num_A_acceptors) num_A_acceptors--;
      i--;
    }
  }
} //end HB_calc::determineAtomsbymass()

void HB_calc::readframe() {

  InG96 icc;

  try {
    args -> check("ref", 1);
    Arguments::const_iterator iterr = args -> lower_bound("ref");
    icc.open((iterr->second).c_str());
  } catch (const Arguments::Exception &) {
    args -> check("traj", 1);
    Arguments::const_iterator iterr = args -> lower_bound("traj");
    icc.open((iterr->second).c_str());
  }
  icc.select("ALL");
  icc >> *sys;
  icc.close();
}//end HB_calc::readframe()

void HB_calc::setval(gcore::System& sys, args::Arguments& args) {
  this->sys = &sys;
  this->args = &args;
  donors = AtomSpecifier(sys);
  bound = AtomSpecifier(sys);
  acceptors = AtomSpecifier(sys);
  time = 0.0;
  determineAtoms();
  //check for massfile
  string mfile;
  if (args.count("massfile") > 0) {
    Arguments::const_iterator iter = args.lower_bound("massfile");
    if (iter != args.upper_bound("massfile")) {
      mfile = (iter->second.c_str());
    }
    readinmasses(mfile);
    determineAtomsbymass();
  }
}//end HB_calc::setval()

void HB_calc::init() {
  pbc = BoundaryParser::boundary(*sys, *args);
  frames = 0;
}//end HB_calc::init()

bool HB_calc::neighbour(int i, int j) {
  return (bound.atom(i) != acceptors.atom(j) ||
          bound.mol(i) != acceptors.mol(j));
}//end HB_calc::neighbour()

bool HB_calc::distances(double &dist, gmath::Vec &tmpA) {
  dist = tmpA.abs2();
  return (max_distance2 >= dist);
}//end HB_calc::distances()

bool HB_calc::angle(int i, double &angles, gmath::Vec &bound_i, gmath::Vec &tmpA) {
  gmath::Vec tmpB;
  tmpB = bound_i - donors.pos(i);
  angles = acos((tmpA.dot(tmpB)) / (tmpA.abs() * tmpB.abs()))*180 / M_PI;
  return (min_angle <= angles);
}//end HB_calc::angle()

void HB_calc::settime(double times) {
  time = times;
}//end HB_calc::settime()

