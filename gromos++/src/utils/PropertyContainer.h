/* 	$Id$	 */

#ifndef MAXFLOAT
#define	MAXFLOAT	((float)3.40282346638528860e+38)
#endif

#ifndef INCLUDED_UTILS_PROPERTYCONTAINER
#define INCLUDED_UTILS_PROPERTYCONTAINER

#ifndef INCLUDED_VECTOR
#include <vector>
#define INCLUDED_VECTOR
#endif

#ifndef INCLUDED_STRING
#include <string>
#define INCLUDED_STRING
#endif

#ifndef INCLUDED_UTILS_PROPERTY
#include "Property.h"
#endif

#ifndef INCLUDED_IOSTREAM
#include <iostream>
#define INCLUDED_IOSTREAM
#endif

#ifndef INCLUDED_GROMOS_EXCEPTION
#include "../gromos/Exception.h"
#endif

#ifndef INCLUDED_GMATH_DISTRIBUTION
#include "../gmath/Distribution.h"
#endif

namespace gcore
{
  class System;
  class Molecule;
  class MoleculeTopology;
}

namespace gmath
{
  class Distribution;
}

namespace utils
{
  /**
   * Class PropertyContainer
   * Purpose: Implements a container for properties.
   *
   * Description:
   * This class implements a container for properties. It parses the
   * arguments and constructs the single properties. It loops over these
   * properties to calculate their values, averages or to print out the
   * information.
   *
   * @class PropertyContainer
   * @version Wed Jul 31 2002
   * @author M. Christen
   * @ingroup utils
   * @sa utils::Property
   */
  
  class PropertyContainer: public std::vector<Property *>
    {
    public:
      /**
       * Constructor.
       */
      PropertyContainer();
      /**
       * Constructor.
       * Most properties need a reference to the system.
       */
      PropertyContainer(gcore::System &sys);
      /**
       * Destructor.
       */
      virtual ~PropertyContainer();
      
      // accessors
      /**
       * All values of the properties can be added to a distribution.
       * This is in general only sensible, if all properties are of the
       * same type. See for an example the dist program.
       */
      void addDistribution(gmath::Distribution &distribution);
      /**
       * Remove the distribution.
       */
      void removeDistribution();
      /**
       * Accessor for the distribution.
       */
      gmath::Distribution &getDistribution();
      
      // methods
      /**
       * Add a property specifier (thereby constructing a property class).
       */
      int addSpecifier(std::string s);
      /**
       * Get title string.
       */
      std::string toTitle();
      /**
       * Get results of a calculation.
       */
      std::string toString();
      /**
       * Calculate all properties in the container.
       */
      void calc();
      /**
       * Check the boundaries of all the properties in the container.
       */
      std::string checkBounds();
      /**
       * Get the average over all the calc() calls. Probably only useful,
       * if all properties are of the same type.
       * Again: see the dist program for an example.
       */
      std::string averageOverRun();
      /**
       * Get the average over all the properties calculated in one calc()
       * call.
       */
      std::string averageOverProperties();
      /**
       * Or if you rather need the values than a string.
       * get average, rmsd from average, rmsd from z-value,
       * lower and upper bound of all properties after each calculation.
       */
      void averageOverProperties(double &av, double &rmsd, 
				 double &zrmsd, double &lb, 
				 double &ub, int &lp, int &up);
      /**
       * @struct Exception
       * PropertyContainer exception.
       */
      struct Exception: public gromos::Exception
      {
	/**
	 * Constructor.
	 */
	Exception(const std::string &what) : 
	  gromos::Exception("PropertyContainer", what) {}
      };

    protected:
      // internal functions
      /**
       * Parse the arguments string into property type and the rest of
       * the arguments.
       */
      void parse(std::string s);
      /**
       * Creates the specified property. Override if you want to add user
       * defined properties (ie properties derived from Property, but not
       * DistanceProperty, AngleProperty or TorsionProperty).
       * The remaining arguments are passed to the property to parse.
       */      
      virtual Property *createProperty(std::string type, int count,
				       std::string arguments[]);

      /**
       * get the argument number that contains the atoms (to check for
       * all molecules.
       * return -1 if no all molecules should be done.
       */
      virtual int all_mol_arg(std::string type);
      
      /**
       * Reference to the system. Often needed to construct a property.
       */
      gcore::System *d_sys;
      /**
       * Distribution that holds the data over all calls to calc() for all
       * the properties.
       * See program dist for an example.
       */
      gmath::Distribution *d_distribution;
     
      // friend for output
      /**
       * Prints all the values of the properties in the container.
       */
      friend std::ostream &operator<<(std::ostream &os, PropertyContainer &s);
    };  

  // inlines
  inline std::string PropertyContainer::checkBounds()
    {
      std::string s="";
      for(iterator it=begin(); it!=end(); ++it)
	s += (*it)->checkBounds();
      return s;
    }  
  
  inline void PropertyContainer::addDistribution(gmath::Distribution &distribution)
    {
      d_distribution = &distribution;
    }
  
  inline void PropertyContainer::removeDistribution()
    {
      d_distribution = NULL;
    }

  inline gmath::Distribution &PropertyContainer::getDistribution()
    {
      return *d_distribution;
    }
  
}

#endif

  
