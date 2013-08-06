#include <cmath>
#include <sstream>

#include "../gromos/Exception.h"
#include "../args/Arguments.h"
#include "../src/gio/InTopology.h"
#include "../src/args/BoundaryParser.h"
#include "../src/bound/Boundary.h"
#include "../gcore/System.h"

#include "SoluteAverageDistance.h"

static const int fgIndex = 0;
static const int cgIndex = 1;

namespace utils {
  SoluteAverageDistance::SoluteAverageDistance(gcore::System& sys, args::Arguments& args) :
          _solute(sys),
          _fgSolvent(sys),
          _cgSolvent(sys),
          _withEnergy(false),
          _sys(sys)
  {
    if ( !(args.count("solute") > 0 &&
         ( args.count("fgsolv") > 0 ||
           args.count("cgsolv") > 0) ) ) {
      throw gromos::Exception("SoluteAverageDistance",
              "You must specify solute and at least fg or cg solvent");
    }
    
    std::string keywords[] = {"solute", "fgsolv", "cgsolv"};
    AtomSpecifier * const as[] = {&_solute, &_fgSolvent, &_cgSolvent};

    for (unsigned int i = 0; i < 3; i++) {
      for (args::Arguments::const_iterator iter = args.lower_bound(keywords[i]),
              to = args.upper_bound(keywords[i]);
              iter != to; ++iter) {
        as[i]->addSpecifier(iter->second);
      }
    }
    _pbc = args::BoundaryParser::boundary(sys, args);
    
    // Should we also calculate the energy?
    if (args.count("energy") == 4){
      args::Arguments::const_iterator iter = args.lower_bound("energy");
      for (unsigned int i = 0; i < 2; i++){
        double params[] = {0.0, 0.0};
        for (unsigned int j = 0; j < 2; j++, iter++){
          std::istringstream is(iter->second);
          is >> params[j];
        }
        SAD_Param p(params[0], params[1]); // first force constant, then cutoff
        _params.push_back(p);
      }
      std::cerr << "# FG: Force Constant = " << _params[fgIndex].forceConstant << std::endl
                << "#     Cut Off        = " << _params[fgIndex].cutoff << std::endl;
      std::cerr << "# CG: Force Constant = " << _params[cgIndex].forceConstant << std::endl
                << "#     Cut Off        = " << _params[cgIndex].cutoff << std::endl;
      _withEnergy = true;
    } else if (args.count("energy") != -1) {
      std::cerr << "# There is probably something wrong with the numbers of "
              << "arguments for energy!" << std::endl;
      std::cerr << "# There were " << args.count("energy") << " given!\n";
    }
  }
  
  bool SoluteAverageDistance::withEnergy() const {
    return _withEnergy;
  }
  
  void SoluteAverageDistance::calculate() {
    
    double n = 6;
    
    AtomSpecifier * const as[] = {&_fgSolvent, &_cgSolvent};
    Distances * const ds[] = {&_fgDistances, &_cgDistances};
    for (unsigned int k = 0; k < 2; k++) {
      ds[k]->clear();
      for (int j = 0; j < as[k]->size(); j++) {
        double sum = 0.0;
        for (int i = 0; i < _solute.size(); i++) {
          gmath::Vec nim_sol = _pbc->nearestImage(as[k]->pos(j), _solute.pos(i),
                  _sys.box());
          gmath::Vec dist = as[k]->pos(j) - nim_sol;
          sum += pow(dist.abs(), -n);
        }
        ds[k]->push_back(pow(sum, -1.0 / n));
      }
    }
    return;
  }
  
  std::string SoluteAverageDistance::title() {
    std::string title("# Time   ");
    
    AtomSpecifier * const as[] = {&_fgSolvent, &_cgSolvent};
    for (unsigned int k = 0; k < 2; k++) {
      for (int i = 0; i < as[k]->size(); i++) {
        title += " ";
        title += as[k]->toString(i);
      }
    }
    title += "\n";
    return title;
  }
  
  void SoluteAverageDistance::distances(std::ostream &os) const {
    const Distances * const ds[] = {&_fgDistances, &_cgDistances};
    for (unsigned int k = 0; k < 2; k++) {
      Distances::const_iterator it = ds[k]->begin(),
              to = ds[k]->end();
      for (; it != to; it++) {
        os << " ";
        os << *it;
      }
    }
  }

  void SoluteAverageDistance::energies(std::ostream &os) const {
    const Distances * const ds[] = {&_fgDistances, &_cgDistances};
    for (unsigned int k = 0; k < 2; k++) {
      Distances::const_iterator it = ds[k]->begin(),
              to = ds[k]->end();
      for (; it != to; it++) {
        os << " ";
        const double dist = *it - _params[k].cutoff;
        if ((k == fgIndex && dist > 0.0) || (k == cgIndex && dist < 0.0)) {
          os << 0.5 * _params[k].forceConstant * dist * dist;
        } else {
          os << 0.0;
        }
      }
    }
  }
  
  std::ostream &operator<<(std::ostream &os, SoluteAverageDistance const & sad){
    sad.distances(os);
    return os;
  }
}
