// some helper classes to construct and write a perturbation topology
namespace utils {

using namespace std;

struct PertAtom {
  unsigned int gromosNum;
  string name;
  unsigned int resNum;
  int iac[2];
  double mass[2];
  double charge[2];
  double alphaLJ;
  double alphaCRF;
  
  PertAtom() : alphaLJ(1.51), alphaCRF(0.5) {}
  PertAtom(double lj, double crf) : alphaLJ(lj), alphaCRF(crf) {}
};

ostream & operator<<(ostream & os, const PertAtom &pa) {
  os.precision(5);
  os.setf(ios::fixed, ios::floatfield);
  os << setw(6) << pa.gromosNum+1 << ' ' << setw(4) << pa.resNum + 1 
     << ' ' << setw(4) << pa.name  << setw(4) << pa.iac[0]+1 << setw(9) 
     << pa.mass[0] << setw(9) << pa.charge[0]
     << "   " 
     << setw(4) << pa.iac[1]+1 << setw(9) << pa.mass[1] << setw(9) << pa.charge[1];
  os.precision(2);
  os << "   " << setw(6) << pa.alphaLJ << setw(6) << pa.alphaCRF << endl;
  return os;
}

struct PertAtomPair {
  unsigned int gromosNum[2];
  unsigned int interactionType[2];
  
  bool operator<(const PertAtomPair& pa) const {
    if (gromosNum[0] == pa.gromosNum[0])
      return gromosNum[1] < pa.gromosNum[1];
    else
      return gromosNum[0] < pa.gromosNum[0];
  }
};

ostream & operator<<(ostream & os, const PertAtomPair &pap) {
  os << setw(6) << pap.gromosNum[0]+1 << setw(6) << pap.gromosNum[1]+1 
     << setw(5) << pap.interactionType[0] << setw(5) << pap.interactionType[1]
     << endl;
  return os;
}

struct PertBond {
  unsigned int gromosNum[2];
  unsigned int type[2];
};

ostream & operator<<(ostream & os, const PertBond &pb) {
  os << setw(6) << pb.gromosNum[0]+1 << setw(6) << pb.gromosNum[1]+1 
     << setw(5) << pb.type[0]+1 << setw(5) << pb.type[1]+1 << endl;
  return os;
}

struct PertAngle {
  unsigned int gromosNum[3];
  unsigned int type[2];
};

ostream & operator<<(ostream & os, const PertAngle &pa) {
  os << setw(6) << pa.gromosNum[0]+1 << setw(6) << pa.gromosNum[1]+1 
     << setw(6) << pa.gromosNum[2]+1
     << setw(5) << pa.type[0]+1 << setw(5) << pa.type[1]+1 << endl;
  return os;
}

struct PertImproper {
  unsigned int gromosNum[4];
  unsigned int type[2];
};

ostream & operator<<(ostream & os, const PertImproper &pi) {
  os << setw(6) << pi.gromosNum[0]+1 << setw(6) << pi.gromosNum[1]+1 
     << setw(6) << pi.gromosNum[2]+1 << setw(6) << pi.gromosNum[3]+1
     << setw(5) << pi.type[0]+1 << setw(5) << pi.type[1]+1 << endl;
  return os;
}

struct PertDihedral {
  unsigned int gromosNum[4];
  unsigned int type[2];
};

ostream & operator<<(ostream & os, const PertDihedral &pd) {
  os << setw(6) << pd.gromosNum[0]+1 << setw(6) << pd.gromosNum[1]+1 
     << setw(6) << pd.gromosNum[2]+1 << setw(6) << pd.gromosNum[3]+1
     << setw(5) << pd.type[0]+1 << setw(5) << pd.type[1]+1 << endl;
  return os;
}

bool exclusionContains(const gcore::Exclusion& ex, int atom) {
  for(int i = 0; i < ex.size(); i++) {
    if (ex.atom(i) == atom)
      return true;
  }
  return false;
}

bool containsGromosNum(const utils::AtomSpecifier& spec, int gromosNum) {
  for(int i = 0; i < spec.size(); i++) {
    if (spec.gromosAtom(i) == gromosNum)
      return true;
  }
  return false;
}

bool containsGromosNum(const utils::AtomSpecifier& spec, const gcore::Bond& bond) {
  return containsGromosNum(spec, bond[0]) ||
         containsGromosNum(spec, bond[1]);
}

bool containsGromosNum(const utils::AtomSpecifier& spec, const gcore::Angle& angle) {
  return containsGromosNum(spec, angle[0]) ||
         containsGromosNum(spec, angle[1]) ||
         containsGromosNum(spec, angle[2]);
}

bool containsGromosNum(const utils::AtomSpecifier& spec, const gcore::Dihedral& dihedral) {
  return containsGromosNum(spec, dihedral[0]) ||
         containsGromosNum(spec, dihedral[1]) ||
         containsGromosNum(spec, dihedral[2]) ||
         containsGromosNum(spec, dihedral[3]);
}

bool containsGromosNum(const utils::AtomSpecifier& spec, const gcore::Improper& improper) {
  return containsGromosNum(spec, improper[0]) ||
         containsGromosNum(spec, improper[1]) ||
         containsGromosNum(spec, improper[2]) ||
         containsGromosNum(spec, improper[3]);
}

}


