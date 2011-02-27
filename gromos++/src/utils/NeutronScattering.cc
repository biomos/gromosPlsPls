#include <cassert>
#include <vector>
#include <map>
#include <set>

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
    System *d_sys;
    multimap<int, int> d_comb;
    map<int, double> d_scattLen;
    map<int, double> d_sigma;
    AtomSpecifier d_centre;
    AtomSpecifier d_with;

  };

  NS::NS(System *sys) {
    d_this = new iNS;
    d_this->d_sys = sys;
    d_this->d_centre.setSystem(*d_this->d_sys);
    d_this->d_with.setSystem(*d_this->d_sys);
  }

  NS::~NS(void) {
    if (d_this) {
      delete d_this;
    }
  }

  void NS::setGrid(int grid) {
    assert(d_this != NULL);
    d_this->d_grid = grid;
    for (unsigned int i = 0; i < d_this->d_rdf.size(); ++i) {
      d_this->d_rdf[i].setGrid(grid);
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
    return d_this->d_comb.size();
  }

  /**
   * Prints the combination of centre to with IAC numbers.
   */
  void NS::printComb() {
    for(multimap<int, int>::const_iterator it = d_this->d_comb.begin();
            it != d_this->d_comb.end(); it++) {
      cout << it->first << "\t" << it->second << endl;
    }
  }

}
