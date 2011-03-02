#include <cassert>
#include <vector>
#include <map>
#include <set>
#include <iomanip>
#include <sstream>

#include "../gcore/LJException.h"
#include "AtomSpecifier.h"
#include "../args/Arguments.h"
#include "../gcore/System.h"
#include "RDF.h"
#include "NeutronScattering.h"

using namespace std;
using namespace args;
using namespace gcore;

namespace utils {

  class iNS {
  public:

    int d_grid;
    double d_Qmax;
    vector<RDF> d_rdf;
    vector<vector<double> > d_Sintra;
    vector<vector<double> > d_Sinter;
    vector<double> IofQ;
    System *d_sys;
    multimap<int, int> d_comb;
    map<int, double> d_scattLen;
    map<int, double> d_sigma;
    AtomSpecifier d_atoms;
    vector<int> d_weightIntra;
    vector<int> d_weightInter;
    map<int, double> d_afraction;

  };

  NS::NS(System *sys, args::Arguments::const_iterator firsttrj,
            args::Arguments::const_iterator lasttrj) {
    d_this = new iNS;
    setSystem(sys);
    setTrajectories(firsttrj, lasttrj);
    d_this->d_grid = 200;
  }

  NS::~NS(void) {
    if (d_this) {
      delete d_this;
    }
  }

  void NS::setGrid(int grid) {
    assert(d_this != NULL);
    assert(d_this->d_rdf.size() == d_this->d_Sinter.size() &&
            d_this->d_rdf.size() == d_this->d_Sintra.size());
    d_this->d_grid = grid;
    for (unsigned int i = 0; i < d_this->d_rdf.size(); ++i) {
      d_this->d_rdf[i].setGrid(grid);
      d_this->d_Sinter[i].resize(grid);
      d_this->d_Sintra[i].resize(grid);
    }
  }

  void NS::setCut(int cut) {
    assert(d_this != NULL);
    for (unsigned int i = 0; i < d_this->d_rdf.size(); ++i) {
      d_this->d_rdf[i].setCut(cut);
    }
  }

  void NS::setQmax(int Qmax) {
    assert(d_this != NULL);
    d_this->d_Qmax = Qmax;
  }

  int NS::addAtoms(string s) {
    assert(d_this != NULL);
    d_this->d_atoms.addSpecifier(s);
    d_this->d_atoms.sort();
    return d_this->d_atoms.size();
  }

  void NS::setSystem(System *sys) {
    assert(d_this != NULL);
    d_this->d_sys = sys;
    d_this->d_atoms.setSystem(*sys);
    for(unsigned int i = 0; i < d_this->d_rdf.size(); ++i) {
      d_this->d_rdf[i].setSystem(sys);
    }
  }

  int NS::getCombinations(void) {
    // make sure there are no (old) combinations in the list now
    d_this->d_comb.clear();
    for(int c = 0; c < d_this->d_atoms.size(); ++c) {
      for(int w = 0; w < d_this->d_atoms.size(); ++w) {
        int iacc = d_this->d_atoms.iac(c); // has to be in this loop since it
                                           // need to be reset in case of
                                           // iacw > iacc, which is switched
                                           // below
        int iacw = d_this->d_atoms.iac(w);
        // to make sure the centres are smaller than the with in the d_comb map...
        if(iacw < iacc) {
          int tmp = iacc;
          iacc = iacw;
          iacw = tmp;
        }
        bool found = false;
        multimap<int, int>::const_iterator start = d_this->d_comb.lower_bound(iacc);
        multimap<int, int>::const_iterator stop = d_this->d_comb.upper_bound(iacc);
        multimap<int, int>::const_iterator it;
        for (it = start; it != stop; ++it) {
          if (it->second == iacw) {
            found = true;
            break;
          }
        }
        if (!found) {
          d_this->d_comb.insert(pair<int, int>(iacc, iacw));
        }
      }
    }
    // sort the multimap of combinations
    map<int, set<int> > comb;
    {
      multimap<int, int>::iterator it;
      for(it = d_this->d_comb.begin(); it != d_this->d_comb.end(); ++it) {
        map<int, set<int> >::iterator iit = comb.find(it->first);
        if(iit != comb.end()) {
          iit->second.insert(it->second);
        } else {
          set<int> s;
          s.insert(it->second);
          comb.insert(pair<int, set<int> >(it->first, s));
        }
      }
      d_this->d_comb.clear();
      map<int, set<int> >::iterator itc;
      for(itc = comb.begin(); itc != comb.end(); itc++) {
        set<int>::iterator its;
        for(its = itc->second.begin(); its != itc->second.end(); ++its) {
          d_this->d_comb.insert(pair<int, int>(itc->first, *its));
        }
      }
    }
    // now resize the depending vector lengths
    d_this->d_Sinter.resize(d_this->d_comb.size());
    d_this->d_Sintra.resize(d_this->d_comb.size());
    d_this->d_rdf.resize(d_this->d_comb.size());
    d_this->d_weightInter.resize(d_this->d_comb.size());
    d_this->d_weightIntra.resize(d_this->d_comb.size());
    // now calculate the mole fractions
    d_this->d_afraction.clear();
    map<int, set<int> >::iterator it;
    for(it = comb.begin(); it != comb.end(); ++it) {
      int num = 0;
      for(int a = 0; a < d_this->d_atoms.size(); ++a) {
        if(it->first == d_this->d_atoms.iac(a)) {
          num++;
        }
      }
      d_this->d_afraction.insert(pair<int, double>(it->first, double(num)/double(d_this->d_atoms.size())));
    }
    return d_this->d_comb.size();
  }

  void NS::check(void) {
    if(!d_this) {
      stringstream msg;
      msg << "inertialisation of the implementation calss iNS failed";
      throw gromos::Exception("class utils::NS", msg.str());
    }
    if(!d_this->d_sys) {
      stringstream msg;
      msg << "no system (gcore::System) set";
      throw gromos::Exception("class utils::NS", msg.str());
    }
  }

  void NS::setTrajectories(args::Arguments::const_iterator firsttrj,
            args::Arguments::const_iterator lasttrj) {
    assert(d_this != NULL);
    for(unsigned int i = 0; i < d_this->d_rdf.size(); i++) {
      d_this->d_rdf[i].setTrajectories(firsttrj, lasttrj);
    }
  }

  void NS::setRDFatoms() {
    assert(d_this != NULL);
    assert(d_this->d_comb.size() == d_this->d_rdf.size() && d_this->d_comb.size() > 0);
    int i = 0;
    multimap<int, int>::iterator it;
    for(it = d_this->d_comb.begin(); it != d_this->d_comb.end(); ++it) {
      for(int asn = 0; asn < d_this->d_atoms.size(); ++asn) {
        if (d_this->d_atoms.iac(asn) == it->first) {
          int m = d_this->d_atoms.mol(asn);
          int a = d_this->d_atoms.atom(asn);
          d_this->d_rdf[i].addCentersAtom(m, a);
        }
        if (d_this->d_atoms.iac(asn) == it->second) {
          int m = d_this->d_atoms.mol(asn);
          int a = d_this->d_atoms.atom(asn);
          d_this->d_rdf[i].addWithAtom(m, a);
        }
      }
      ++i;
    }
  }

  void NS::calcRDFsInterAll() {
    assert(d_this != NULL);
    assert(d_this->d_comb.size() == d_this->d_rdf.size());
    for(unsigned int i = 0; i < d_this->d_rdf.size(); ++i) {
      d_this->d_rdf[i].calculateInter();
    }
  }

  void NS::printRDFs(std::ostream &os) {
    for(unsigned int i = 0; i < d_this->d_rdf.size(); ++i) {
      d_this->d_rdf[i].print(os);
      os << endl;
    }
  }

  void NS::getWeights(void) {
    assert(d_this != NULL);
    assert(d_this->d_atoms.size() != 0);
    int i = 0;
    multimap<int, int>::iterator it;
    for (it = d_this->d_comb.begin(); it != d_this->d_comb.end(); ++it) {
      // make sure there are no old weights
      d_this->d_weightIntra[i] = 0;
      d_this->d_weightInter[i] = 0;
      // get the intramolecular weights
      for (int c = 0; c < d_this->d_atoms.size(); ++c) {
        for (int w = 0; w < d_this->d_atoms.size(); ++w) {
          int iacc = d_this->d_atoms.iac(c); // has to be in this loop since it
                                             // need to be reset in case of
                                             // iacw > iacc, which is switched
                                             // below
          int iacw = d_this->d_atoms.iac(w);
          // make sure iacc <= iacw
          if (iacw < iacc) {
            int tmp = iacc;
            iacc = iacw;
            iacw = tmp;
          }
          // skip if centre and with atom are in different molecules
          if (d_this->d_atoms.mol(c) != d_this->d_atoms.mol(w)) {
            continue;
          }
          // skip if the same molecule but also the same atom
          if(d_this->d_atoms.mol(c) == d_this->d_atoms.mol(w) &&
                  d_this->d_atoms.atom(c) == d_this->d_atoms.atom(w)) {
            continue;
          }
          // skipt if centre or with atoms have not the right type (IAC)
          if(iacc != it->first || iacw != it->second) {
            continue;
          }
          d_this->d_weightIntra[i]++;
        }
      }
      // and now the inter-molecular weights
      for (int c = 0; c < d_this->d_atoms.size(); ++c) {
        int iacc = d_this->d_atoms.iac(c);
        for (int w = 0; w < d_this->d_atoms.size(); ++w) {
          int iacw = d_this->d_atoms.iac(w);
          // we only count if iacc >= iacw
          if(iacw < iacc) {
            continue;
          }
          // skip if centre and with atom are in the same molecule
          if (d_this->d_atoms.mol(c) == d_this->d_atoms.mol(w)) {
            continue;
          }
          // skipt if centre or with atoms have not the right type (IAC)
          if (iacc != it->first || iacw != it->second) {
            continue;
          }
          // set the f constant to 1 if centre and with of the same type
          // (to avoid double counting)
          double f = 0.0;
          if (iacc == iacw) {
            f = 1.0;
          }
          d_this->d_weightInter[i] += (2.0 - f);
        }
      }
      ++i;
    }
  }

  /**
   * Prints the combination of centre to with IAC numbers.
   */
  void NS::print(ostream &os) {

    os << "NEUTRONSCATTERING\n";
    os << "# IACI: IAC number of atom i\n";
    os << "# IACJ: IAC number of atom j\n";
    os << "# WIJM: weight of the corresponding intra-molecular structure factor\n";
    os << "# WIJD: weight of the corresponding inter-molecular structure factor\n";
    os << "# NC  : number of IAC combinations (centre-to-with atoms)\n";
    os << "#\n";
    os << "#" << setw(7) << "NC" << endl;
    os << setw(8) << d_this->d_comb.size() << endl;
    os << "#" << setw(7) << "IACI" << setw(8) << "IACJ"
            << setw(20) << "WIJM" << setw(20) << "WIJD" << endl;
    os.precision(9);
    int i = 0;
    for (multimap<int, int>::iterator it = d_this->d_comb.begin();
            it != d_this->d_comb.end(); ++it) {
      os << setw(8) << (it->first) + 1 << setw(8) << (it->second) + 1<< scientific
              << setw(20) << d_this->d_weightIntra[i]
              << setw(20) << d_this->d_weightInter[i] << endl;
      ++i;
    }
    os << "# NDAT: Number of different atom types\n";
    os << "# IAC:  Integer atom code (atom tupe)\n";
    os << "# AFR:  atom fraction (AFR = NATT / TNA)\n";
    os << "#\n";
    os << "#" << setw(14) << "NDAT" << setw(15) << endl;
    os << setw(15) << d_this->d_afraction.size() << endl;
    os << "#" << setw(14) << "IAC" << setw(15) << setw(20) << "AFR" << endl;
    for(map<int, double>::iterator it = d_this->d_afraction.begin();
            it != d_this->d_afraction.end(); ++it) {
      os << setw(15) << (it->first) + 1
              << setw(20) << scientific << it->second << endl;
    }
    os << "END\n";
  }

}
