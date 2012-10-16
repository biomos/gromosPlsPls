#include <cassert>
#include <set>
#include <string>
#include <iostream>
#include <vector>
#include <iomanip>

#include "../gcore/System.h"
#include "AtomSpecifier.h"
#include "../gio/InG96.h"
#include "../args/Arguments.h"
#include "../gcore/Box.h"
#include "../gmath/Vec.h"
#include "../bound/Boundary.h"
#include "../args/BoundaryParser.h"
#include "../gmath/Distribution.h"
#include "../bound/Boundary.h"
#include "../args/GatherParser.h"

#include <RDF.h>

using namespace std;

namespace utils {

  /**
   * Class iRDF is the implementation class of the class RDF, i.e. it stores
   * all data elements of the RDF class.
   *
   * @class iDRF
   * @ingroup utils
   * @author A. Eichenberger
   */
  class iRDF {
  public:
    /**
     * a pointer to GROMOS++like arguments
     */
    const args::Arguments *d_args;
    /**
     * Stores the radial distribution function.
     */
    std::vector<double> d_rdf;
    /**
     * Stores the local number of "center" with "with"
     */
    std::vector<double> d_local_mix;
    /**
     * Stores the local number of "center" with "center"
     */
    std::vector<double> d_local_self;
    /**
     * A pointer to the system for which the radial distribution function is
     * calculated.
     */
    gcore::System *d_sys;
    /**
     * An atom specifier to store the centre atoms for the rdf calculation
     */
    AtomSpecifier d_centre;
    /**
     * An atom specifier to store the with atoms for the rdf calculation
     */
    AtomSpecifier d_with;
    /**
     * The grid number, i.e. the resolution of the calculated rdf function
     */
    unsigned int d_grid;
    /**
     * The cutoff distance: the maximum distance to be considered for the radial
     * distribution function calculation.
     */
    double d_cut;

  };

  RDF::RDF() {
    d_this = new iRDF;
    d_this->d_grid = 200;
    d_this->d_cut = 1.5;
    d_this->d_rdf.resize(d_this->d_grid);
    d_this->d_local_mix.resize(d_this->d_grid);
    d_this->d_local_self.resize(d_this->d_grid);
  }

  RDF::RDF(gcore::System *sys, const args::Arguments *args) {
    d_this = new iRDF;
    d_this->d_sys = sys;
    d_this->d_grid = 200;
    d_this->d_cut = 1.5;
    d_this->d_rdf.resize(d_this->d_grid);
    d_this->d_local_mix.resize(d_this->d_grid);
    d_this->d_local_self.resize(d_this->d_grid);
    d_this->d_centre.setSystem(*d_this->d_sys);
    d_this->d_with.setSystem(*d_this->d_sys);
    d_this->d_args = args;
  }

  RDF::RDF(const RDF &rdf) {
    d_this = new iRDF;
    *d_this = *rdf.d_this;
  }

  RDF::~RDF(void) {
    assert(d_this != NULL);
    delete d_this;
  }

  int RDF::addCenters(string s) {
    assert(d_this != NULL);
    return d_this->d_centre.addSpecifier(s);
  }

  void RDF::addCentersAtom(int m, int a) {
    assert(d_this != NULL);
    d_this->d_centre.addAtom(m, a);
  }

  void RDF::addWithAtom(int m, int a) {
    assert(d_this != NULL);
    d_this->d_with.addAtom(m, a);
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

  void RDF::clearLocal(void) {
    assert(d_this != NULL);
    for(unsigned int i = 0; i < d_this->d_local_mix.size(); ++i) {
      d_this->d_local_mix[i] = 0.0;
      d_this->d_local_self[i] = 0.0;
    }
  }

  void RDF::setSystem(gcore::System *sys) {
    assert(d_this != NULL);
    d_this->d_sys = sys;
    d_this->d_centre.setSystem(*d_this->d_sys);
    d_this->d_with.setSystem(*d_this->d_sys);
  }

  void RDF::calculateAll(void) {
    assert(d_this != NULL && d_this->d_sys != NULL);

    clearRDF();

    gio::InG96 ic;
    unsigned int count_frame = 0;

    double correct=4*acos(-1.0)*d_this->d_cut/double(d_this->d_grid);
    
    // boundary conditions
    bound::Boundary *pbc;
    // gather method
    bound::Boundary::MemPtr gathmethod;
    
    // loop over the different trajectory files
    for (args::Arguments::const_iterator trj = d_this->d_args->lower_bound("traj"); trj != d_this->d_args->upper_bound("traj"); trj++) {
      
      // reading the boundary shape and gathering method
      ic.open(trj->second.c_str());
      ic.select("ALL");
      ic >> *(d_this->d_sys);
      // here we have to check whether we really got the atoms we want
      // maybe the solvent is missing.
      if (d_this->d_centre.size() == 0 || d_this->d_with.size() == 0) {
        string argument = d_this->d_centre.size() == 0 ? "centre" : "width";
        throw gromos::Exception("Rdf.cc", "No atoms specified for " + argument + " atoms!");
      }
      // parse boundary conditions
      pbc = args::BoundaryParser::boundary(*d_this->d_sys, *d_this->d_args);
      //parse gather method
      gathmethod = args::GatherParser::parse(*d_this->d_sys, *d_this->d_sys, *d_this->d_args);
      ic.close();

      // reopen the same file for the calculation
      ic.open(trj->second.c_str());
      ic.select("ALL");

      // loop over frames of the current trajectory file
      while (!ic.eof()) {

        // read the next frame
        ic >> *(d_this->d_sys);
        count_frame++;
        // gather the system
        (*pbc.*gathmethod)();

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

    // boundary conditions
    bound::Boundary *pbc;
    // gather method
    bound::Boundary::MemPtr gathmethod;
    
    // loop over the different trajectory files
    for (args::Arguments::const_iterator trj = d_this->d_args->lower_bound("traj"); trj != d_this->d_args->upper_bound("traj"); trj++) {

      ic.open(trj->second.c_str());
      ic.select("ALL");
      ic >> *(d_this->d_sys);
      // here we have to check whether we really got the atoms we want
      // maybe the solvent is missing.
      if (d_this->d_centre.size() == 0 || d_this->d_with.size() == 0) {
        string argument = d_this->d_centre.size() == 0 ? "centre" : "width";
        throw gromos::Exception("Rdf.cc", "No atoms specified for " + argument + " atoms!");
      }
      // parse boundary conditions
      pbc = args::BoundaryParser::boundary(*d_this->d_sys, *d_this->d_args);
      //parse gather method
      gathmethod = args::GatherParser::parse(*d_this->d_sys, *d_this->d_sys, *d_this->d_args);
      ic.close();

      // reopen the same file for the calculation
      ic.open(trj->second.c_str());
      ic.select("ALL");

      // loop over frames of the current trajectory file
      while (!ic.eof()) {

        // read the next frame
        ic >> *(d_this->d_sys);
        count_frame++;
        // gather the system
        (*pbc.*gathmethod)();

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

  void RDF::calculateLocal() {
    assert(d_this != NULL && d_this->d_sys != NULL);

    clearLocal();
    
    gio::InG96 ic;
    unsigned int count_frame = 0;

    // loop over the different trajectory files
    for (args::Arguments::const_iterator trj = d_this->d_args->lower_bound("traj"); trj != d_this->d_args->upper_bound("traj"); trj++) {
      
      // open the trajectory file for reading the boxshape only
      // it is faster to do it here, close the file and reopen it later again
      // for the calculation
      ic.open(trj->second.c_str());
      ic.select("ALL");
      ic >> *(d_this->d_sys);
      // here we have to check whether we really got the atoms we want
      // maybe the solvent is missing.
      if (d_this->d_centre.size() == 0 || d_this->d_with.size() == 0) {
        string argument = d_this->d_centre.size() == 0 ? "centre" : "width";
        throw gromos::Exception("Rdf.cc", "No atoms specified for " + argument + " atoms!");
      }
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

        // loop over the center atoms
#ifdef OMP
#pragma omp parallel for
#endif
        for (int c = 0; c < d_this->d_centre.size(); c++) {

          // the distribution array
          gmath::Distribution dist(0, d_this->d_cut, d_this->d_grid);
          gmath::Distribution dist2(0, d_this->d_cut, d_this->d_grid);

          // the coordinates of the center atom
          const gmath::Vec & centre_coord = *(d_this->d_centre.coord(c));

          // loop over the with atoms
          for (int w = 0; w < d_this->d_with.size(); w++) {
            // only do the calculations if the center and with atom are not identical
            if (!(d_this->d_with.mol(w) == d_this->d_centre.mol(c) && d_this->d_with.atom(w) == d_this->d_centre.atom(c))) {
              const Vec & tmp = pbc->nearestImage(centre_coord, *(d_this->d_with.coord(w)), d_this->d_sys->box());
              dist.add((tmp - centre_coord).abs());
            }
          } /* end of loop over with atoms */
          
          // loop over the center atoms
          for (int w = 0; w < d_this->d_centre.size(); w++) {
            // only do the calculations if the center and with atom are not identical
            if (!(d_this->d_centre.mol(w) == d_this->d_centre.mol(c) && d_this->d_centre.atom(w) == d_this->d_centre.atom(c))) {
              const Vec & tmp = pbc->nearestImage(centre_coord, *(d_this->d_centre.coord(w)), d_this->d_sys->box());
              dist2.add((tmp - centre_coord).abs());
            }
          } /* end of loop over center atoms */
          
          for (unsigned int k = 0; k < d_this->d_grid; k++) {
#ifdef OMP
#pragma omp critical
#endif
            {
              d_this->d_local_mix[k] += double(dist[k]);
              d_this->d_local_self[k] += double(dist2[k]);
            }
          }
        } /* end of loop over center atoms */
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
      d_this->d_local_mix[i] /= double(divide);
      d_this->d_local_self[i] /= double(divide);
    }
    
  } /* end of RDF::calculateLocal() */

  double RDF::calculateInterPDens(unsigned int numatoms) {
    assert(d_this != NULL && d_this->d_sys != NULL);

    // remember the average centre to with distance, which is needed for
    // the calculation of the inelastic (damped) neutron scattering intensities
    // (term of the Debye-Waller factor)
    double centre2with = 0.0;
    int countC2W = 0;

    clearRDF();

    gio::InG96 ic;
    unsigned int count_frame = 0;

    double correct=4*acos(-1.0)*d_this->d_cut/double(d_this->d_grid);

    // boundary conditions
    bound::Boundary *pbc;
    // gather method
    bound::Boundary::MemPtr gathmethod;
    
    // loop over the different trajectory files
    for (args::Arguments::const_iterator trj = d_this->d_args->lower_bound("traj"); trj != d_this->d_args->upper_bound("traj"); trj++) {

      ic.open(trj->second.c_str());
      ic >> *(d_this->d_sys);
      // here we have to check whether we really got the atoms we want
      // maybe the solvent is missing.
      if (d_this->d_centre.size() == 0 || d_this->d_with.size() == 0) {
        string argument = d_this->d_centre.size() == 0 ? "centre" : "width";
        throw gromos::Exception("Rdf.cc", "No atoms specified for " + argument + " atoms!");
      }
      // parse boundary conditions
      pbc = args::BoundaryParser::boundary(*d_this->d_sys, *d_this->d_args);
      //parse gather method
      gathmethod = args::GatherParser::parse(*d_this->d_sys, *d_this->d_sys, *d_this->d_args);
      ic.close();

      // reopen the same file for the calculation
      ic.open(trj->second.c_str());
      ic.select("ALL");

      // loop over frames of the current trajectory file
      while (!ic.eof()) {

        // read the next frame
        ic >> *(d_this->d_sys);
        count_frame++;
        // gather the system
        (*pbc.*gathmethod)();

        // calculate the volume
        double vol_corr = 1;
        if(pbc->type()=='t') vol_corr=0.5;
        double vol = d_this->d_sys->box().K_L_M() * vol_corr;
        double partDens = (double)numatoms / vol;

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
            } else {
              if(!(d_this->d_with.atom(w) == d_this->d_centre.atom(c))) {
                const Vec & tmp = pbc->nearestImage(centre_coord, *(d_this->d_with.coord(w)), d_this->d_sys->box());
                centre2with += (tmp - centre_coord).abs();
                countC2W++;
              }
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
              d_this->d_rdf[k] += (rdf_val - 1.0) * partDens;
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

    return centre2with/countC2W;

  } /* end of RDF::calculateInterPDens() */

  void RDF::print(std::ostream &os) {
    os.precision(9);
    for (unsigned int i = 0; i < d_this->d_grid; i++) {
      double r = (double(i) + 0.5) * d_this->d_cut / d_this->d_grid;
      os << setw(15) << r << setw(15) << d_this->d_rdf[i] << endl;
    }
  }

  void RDF::printLocal(std::ostream &os) {
    os.precision(9);
    for (unsigned int i = 0; i < d_this->d_grid; i++) {
      double r = (double(i) + 0.5) * d_this->d_cut / d_this->d_grid;
      double frac = d_this->d_local_self[i] / (d_this->d_local_self[i] + d_this->d_local_mix[i]);
      if (frac >= 0.0)
        os << setw(15) << r << setw(15) << frac << endl;
      else 
        os << setw(15) << r << setw(15) << "0" << endl;
    }
  }

  double RDF::rdf(unsigned int i) {
    assert(d_this->d_rdf.size() > i);
    return d_this->d_rdf[i];
  }


} /* end of namespace utils */
