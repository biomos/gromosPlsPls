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
  
  class PropertyContainer: public std::vector<Property *>
    {
    public:
      PropertyContainer();
      PropertyContainer(gcore::System &sys);
      virtual ~PropertyContainer();
      
      // accessors
      void addDistribution(gmath::Distribution &distribution);
      void removeDistribution();
      gmath::Distribution &getDistribution();
      
      // methods
      int addSpecifier(std::string s);
      std::string toTitle();
      std::string toString();

      void calc();
      std::string checkBounds();

      std::string averageOverRun();
      std::string averageOverProperties();
      
      struct Exception: public gromos::Exception
      {
	  Exception(const string &what) : gromos::Exception("PropertyContainer", what) {}
      };

    protected:
      // internal functions
      void parse(std::string s);
      virtual Property *createProperty(std::string type, int count,
				       std::string arguments[]);
      // member variables
      gcore::System *d_sys;
      gmath::Distribution *d_distribution;
     
      // friend for output
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

  
