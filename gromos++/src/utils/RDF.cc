
#include <cassert>
#include <set>
#include <string>
#include <iostream>
#include <vector>

#include "../gcore/System.h"
#include "AtomSpecifier.h"
#include "../gio/InG96.h"
#include "../args/Arguments.h"
#include "../gcore/Box.h"
#include "../gmath/Vec.h"
#include "../bound/Boundary.h"
#include "../args/BoundaryParser.h"
#include "../gmath/Distribution.h"

#include <RDF.h>

using namespace std;

namespace utils {

  class iRDF {
  public:
    std::vector<double> d_rdf;
    gcore::System *d_sys;
    args::Arguments::const_iterator d_firsttrj;
    args::Arguments::const_iterator d_lasttrj;
    AtomSpecifier d_centre;
    AtomSpecifier d_with;
    unsigned int d_grid;
    double d_cut;

  };

  // definition of the (standard) construcor

  RDF::RDF(gcore::System *sys,
          args::Arguments::const_iterator firsttrj, args::Arguments::const_iterator lasttrj) {
    d_this = new iRDF;
    d_this->d_sys = sys;
    d_this->d_grid = 200;
    d_this->d_cut = 1.5;
    d_this->d_firsttrj = firsttrj;
    d_this->d_lasttrj = lasttrj;
    d_this->d_rdf.resize(d_this->d_grid);
    d_this->d_centre.setSystem(*d_this->d_sys);
    d_this->d_with.setSystem(*d_this->d_sys);
  }

  // definition of the copy constructor

  RDF::RDF(RDF &rdf) {
    d_this = new iRDF;
    *d_this = *rdf.d_this;
  }

  // definition of the destructor

  RDF::~RDF(void) {
    assert(d_this != NULL);
    delete d_this;
  }

  int RDF::addCenters(string s) {
    assert(d_this != NULL);
    return d_this->d_centre.addSpecifier(s);
  }

  void RDF::clearCenters(void) {
    assert(d_this != NULL);
    d_this->d_centre.clear();
  }

  int RDF::addWiths(string s) {
    assert(d_this != NULL);
    return d_this->d_with.addSpecifier(s);
  }

  void RDF::clearWiths(void) {
    assert(d_this != NULL);
    d_this->d_with.clear();
  }

  void RDF::setGrid(unsigned int grid) {
    assert(d_this != NULL);
    d_this->d_grid = grid;
    d_this->d_rdf.resize(d_this->d_grid);
  }

  void RDF::setCut(double cut) {
    assert(d_this != NULL);
    d_this->d_cut = cut;
  }

  void RDF::clearRDF(void) {
    assert(d_this != NULL);
    for(unsigned int i = 0; i < d_this->d_rdf.size(); ++i) {
      d_this->d_rdf[i] = 0.0;
    }
  }

  void RDF::calculateAll(void) {
    assert(d_this != NULL && d_this->d_sys != NULL);

    clearRDF();

    gio::InG96 ic;
    unsigned int count_frame = 0;

    double correct=4*acos(-1.0)*d_this->d_cut/double(d_this->d_grid);

    // loop over the different trajectory files
    for (args::Arguments::const_iterator trj = d_this->d_firsttrj; trj != d_this->d_lasttrj; trj++) {
      
      // open the trajectory file for reading the boxshape only
      // it is faster to do it here, close the file and reopen it later again
      // for the calculation
      ic.open(trj->second.c_str());
      ic >> *(d_this->d_sys);
      // get the boundary from the read box format
      bound::Boundary *pbc = args::BoundaryParser::boundary(*d_this->d_sys);
      ic.close();

      // reopen the same file for the calculation
      ic.open(trj->second.c_str());
      ic.select("ALL");

      // loop over frames of the current trajectory file
      while (!ic.eof()) {

        // read the next frame
        ic >> *(d_this->d_sys);
        count_frame++;

        // calculate the volume
        double vol_corr = 1;
        if(pbc->type()=='t') vol_corr=0.5;
        double vol = d_this->d_sys->box().K_L_M() * vol_corr;

        // loop over the centre atoms
#ifdef OMP
#pragma omp parallel for
#endif
        for (int c = 0; c < d_this->d_centre.size(); c++) {

          // the distribution array
          gmath::Distribution dist(0, d_this->d_cut, d_this->d_grid);

          // the coordinates of the centre atom
          const gmath::Vec & centre_coord = *(d_this->d_centre.coord(c));

          // to know if this atom is also in the with set.
          int inwith = 0;
          if (d_this->d_with.findAtom(d_this->d_centre.mol(c), d_this->d_centre.atom(c))>-1) inwith = 1;

          // loop over the width atoms
          for (int w = 0; w < d_this->d_with.size(); w++) {

            // only do the calculations if the centre and with atom are not identical
            if (!(d_this->d_with.mol(w) == d_this->d_centre.mol(c) && d_this->d_with.atom(w) == d_this->d_centre.atom(c))) {
              const Vec & tmp = pbc->nearestImage(centre_coord, *(d_this->d_with.coord(w)), d_this->d_sys->box());
              dist.add((tmp - centre_coord).abs());
            }

          } /* end of loop over with atoms */

          const double dens = (d_this->d_with.size() - inwith) / vol;
          for (unsigned int k = 0; k < d_this->d_grid; k++) {
            const double r = dist.value(k);
            const double rdf_val = double(dist[k]) / (dens * correct * r * r);
#ifdef OMP
#pragma omp critical
#endif
            {
              d_this->d_rdf[k] += rdf_val;
            }
          }

        } /* end of loop over centre atoms */

      } /* end of loop over frames */

      // close the current trajectory file
      ic.close();

    } /* end of loop over trajectory files */

    // correct the distribution for the number of frames and the number of centre atoms
    int divide = count_frame * d_this->d_centre.size();
#ifdef OMP
#pragma omp parallel for
#endif 
    for (unsigned int i = 0; i < d_this->d_grid; i++) {
      d_this->d_rdf[i] /= double(divide);
    }

  } /* end of RDF::calculateAll() */

  void RDF::calculateInter(void) {
    assert(d_this != NULL && d_this->d_sys != NULL);

    clearRDF();

    gio::InG96 ic;
    unsigned int count_frame = 0;

    double correct=4*acos(-1.0)*d_this->d_cut/double(d_this->d_grid);

    // loop over the different trajectory files
    for (args::Arguments::const_iterator trj = d_this->d_firsttrj; trj != d_this->d_lasttrj; trj++) {

      // open the trajectory file for reading the boxshape only
      // it is faster to do it here, close the file and reopen it later again
      // for the calculation
      ic.open(trj->second.c_str());
      ic >> *(d_this->d_sys);
      // get the boundary from the read box format
      bound::Boundary *pbc = args::BoundaryParser::boundary(*d_this->d_sys);
      ic.close();

      // reopen the same file for the calculation
      ic.open(trj->second.c_str());
      ic.select("ALL");

      // loop over frames of the current trajectory file
      while (!ic.eof()) {

        // read the next frame
        ic >> *(d_this->d_sys);
        count_frame++;

        // calculate the volume
        double vol_corr = 1;
        if(pbc->type()=='t') vol_corr=0.5;
        double vol = d_this->d_sys->box().K_L_M() * vol_corr;

        // loop over the centre atoms
#ifdef OMP
#pragma omp parallel for
#endif
        for (int c = 0; c < d_this->d_centre.size(); c++) {

          // the distribution array
          gmath::Distribution dist(0, d_this->d_cut, d_this->d_grid);

          // the coordinates of the centre atom
          const gmath::Vec & centre_coord = *(d_this->d_centre.coord(c));

          // to know if this atom is also in the with set.
          int inwith = 0;
          if (d_this->d_with.findAtom(d_this->d_centre.mol(c), d_this->d_centre.atom(c))>-1) inwith = 1;

          // loop over the width atoms
          for (int w = 0; w < d_this->d_with.size(); w++) {

            // only do the calculations if the centre and with atom are within different molecules
            if (!(d_this->d_with.mol(w) == d_this->d_centre.mol(c))) {
              const Vec & tmp = pbc->nearestImage(centre_coord, *(d_this->d_with.coord(w)), d_this->d_sys->box());
              dist.add((tmp - centre_coord).abs());
            }

          } /* end of loop over with atoms */

          const double dens = (d_this->d_with.size() - inwith) / vol;
          for (unsigned int k = 0; k < d_this->d_grid; k++) {
            const double r = dist.value(k);
            const double rdf_val = double(dist[k]) / (dens * correct * r * r);
#ifdef OMP
#pragma omp critical
#endif
            {
              d_this->d_rdf[k] += rdf_val;
            }
          }

        } /* end of loop over centre atoms */

      } /* end of loop over frames */

      // close the current trajectory file
      ic.close();

    } /* end of loop over trajectory files */

    // correct the distribution for the number of frames and the number of centre atoms
    int divide = count_frame * d_this->d_centre.size();
#ifdef OMP
#pragma omp parallel for
#endif
    for (unsigned int i = 0; i < d_this->d_grid; i++) {
      d_this->d_rdf[i] /= double(divide);
    }

  } /* end of RDF::calculateInter() */

  void RDF::print(std::ostream &os) {
    for (unsigned int i = 0; i < d_this->d_grid; i++) {
      double r = (double(i) + 0.5) * d_this->d_cut / d_this->d_grid;
      os << r << "\t" << d_this->d_rdf[i] << endl;
    }
  }


} /* end of namespace utils */
