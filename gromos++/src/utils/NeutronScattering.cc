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
    double d_cut;
    double d_Qmax;
    vector<RDF> d_rdf;
    vector<vector<double> > d_Sintra;
    vector<vector<double> > d_Sinter;
    vector<double> IofQ;
    System *d_sys;
    multimap<int, int> d_comb;
    map<int, double> d_scattLen;
    map<int, double> d_sigma;
    AtomSpecifier d_centre;
    AtomSpecifier d_with;

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

  int NS::addCenters(string s) {
    assert(d_this != NULL);
    for (unsigned int i = 0; i < d_this->d_rdf.size(); ++i) {
      d_this->d_rdf[i].addCenters(s);
    }
    return d_this->d_centre.addSpecifier(s);
  }

  int NS::addWiths(string s) {
    assert(d_this != NULL);
    for (unsigned int i = 0; i < d_this->d_rdf.size(); ++i) {
      d_this->d_rdf[i].addWiths(s);
    }
    return d_this->d_with.addSpecifier(s);
  }

  void NS::setSystem(System *sys) {
    assert(d_this != NULL);
    d_this->d_sys = sys;
    d_this->d_centre.setSystem(*sys);
    d_this->d_with.setSystem(*sys);
    for(unsigned int i = 0; i < d_this->d_rdf.size(); ++i) {
      d_this->d_rdf[i].setSystem(sys);
    }
  }

  int NS::getCombinations(void) {
    set<int> centres;
    set<int> withs;
    for (int c = 0; c < d_this->d_centre.size(); ++c) {
      centres.insert(d_this->d_centre.iac(c));
    }
    for (int w = 0; w < d_this->d_with.size(); ++w) {
      withs.insert(d_this->d_with.iac(w));
    }
    for(set<int>::const_iterator cit = centres.begin(); cit != centres.end(); cit++) {
      for(set<int>::const_iterator wit = withs.begin(); wit != withs.end(); wit++) {
        bool found = false;
        multimap<int, int>::const_iterator start = d_this->d_comb.lower_bound(*cit);
        multimap<int, int>::const_iterator stop = d_this->d_comb.upper_bound(*cit);
        multimap<int, int>::const_iterator it;
        for (it = start; it != stop; it++) {
          if (it->second == *wit) {
            found = true;
            break;
          }
        }
        start = d_this->d_comb.lower_bound(*wit);
        stop = d_this->d_comb.upper_bound(*wit);
        for (it = start; it != stop; it++) {
          if (it->second == *cit) {
            found = true;
            break;
          }
        }
        if (!found) {
          d_this->d_comb.insert(pair<int, int>(*cit, *wit));
        }
      }
    }
    // sort the combinations to fulfill c <= w
    multimap<int, int>::iterator it;
    for(it = d_this->d_comb.begin(); it != d_this->d_comb.end(); it++) {
      if(it->second < it->first) {
        int iaci = it->second;
        int iacj = it->first;
        d_this->d_comb.erase(it);
        d_this->d_comb.insert(pair<int, int>(iaci, iacj));
      }
    }
    // now resize the depending vector lengths
    d_this->d_Sinter.resize(d_this->d_comb.size());
    d_this->d_Sintra.resize(d_this->d_comb.size());
    d_this->d_rdf.resize(d_this->d_comb.size());
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

  /**
   * Prints the combination of centre to with IAC numbers.
   */
  void NS::printComb(ostream &os) {

    os << "IACCOMBINATIONS\n";
    os << "# IACI: IAC number of atom i\n";
    os << "# IACJ: IAC number of atom j\n";
    os << "# NC  : number of IAC combinations (centre-to-with atoms)\n";
    os << "#\n";
    os << "#" << setw(7) << "NC" << endl;
    os << setw(8) << d_this->d_comb.size() << endl;
      os << "#" << setw(7) << "IACI" << setw(8) << "IACJ" << endl;
    for(multimap<int, int>::const_iterator it = d_this->d_comb.begin();
            it != d_this->d_comb.end(); it++) {
      os << setw(8) << it->first << setw(8) << it->second << endl;
    }
    os << "END\n";
  }

}
