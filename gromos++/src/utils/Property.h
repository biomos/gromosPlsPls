/* 	$Id$	 */

#ifndef MAXFLOAT
#define	MAXFLOAT	((float)3.40282346638528860e+38)
#endif

#ifndef INCLUDED_UTILS_PROPERTY
#define INCLUDED_UTILS_PROPERTY

#ifndef INCLUDED_VECTOR
#include <vector>
#define INCLUDED_VECTOR
#endif

#ifndef INCLUDED_STRING
#include <string>
#define INCLUDED_STRING
#endif

#ifndef INCLUDED_GROMOS_EXCEPTION
#include "../gromos/Exception.h"
#endif

namespace gcore
{
  class System;
  class Molecule;
  class MoleculeTopology;
}

namespace utils
{
  // general property class

  class Property;
  std::ostream &operator<<(std::ostream &os, Property &p);

  class Property
  {
  public:
    static const int MAXARGUMENTS = 10;

    Property(); // for user defined properties that do not need a system
    Property(gcore::System &sys);
    virtual ~Property();

    // accessor
    float getValue();
    float getZValue();
    
    // methods
    virtual float calc();
    virtual std::string checkBounds();
    virtual std::string toTitle();
    virtual std::string toString();
    virtual std::string average();

    struct Exception: public gromos::Exception
    {
      Exception(const string &what): gromos::Exception("Property", what) {}
    };
    
  protected:
    void parse(int count, std::string arguments[]);
    void parseAtoms(std::string atoms);                      // parses AtomSpecifier strings
    void _parseAtomsHelper(std::string substring, int &mol);

  protected:
    // member variables
    int REQUIREDARGUMENTS; 
    std::string d_type;    // property name (usually first part of class name)

    // bounds for 'boundary violation check'
    float d_ubound, d_lbound;
    // zero value (ie equilibrium/zero energy/standard value of property)
    float d_zvalue;
    // result from calc(), for printing,...
    // if this cannot be used (for a user defined property) the virtual 
    // functions have to be overwritten (to produce sensible output)
    float d_value;
    
    std::vector<int> d_atom;
    std::vector<int> d_mol;

    gcore::System *d_sys;

  };

  
// distance property class
  class DistanceProperty : public Property
  {    
  public:
    DistanceProperty(gcore::System &sys);
    virtual ~DistanceProperty();

    // overwritten to check arguments
    void parse(int count, std::string arguments[]);

    virtual float calc();
    virtual std::string average();  // averaging not implemented in class Property
    virtual std::string toString();
    virtual std::string toTitle();

    struct Exception: public gromos::Exception
    {
      Exception(const string &what): gromos::Exception("DistanceProperty", what) {}
    };

  protected:
    float d_average;
    float d_zrmsd;
    int d_count;
    
  };

// bond angle property class
  class AngleProperty : public Property
    {
    public:
      static const float pi = 3.14159265359;
      
      AngleProperty(gcore::System &sys);
      virtual ~AngleProperty();
      
      void parse(int count, std::string arguments[]);
      
      virtual float calc();
      virtual std::string average();
      virtual std::string toString();
      virtual std::string toTitle();
      
      struct Exception: public gromos::Exception
      {
	  Exception(const string &what): gromos::Exception("AngleProperty", what) {}
      };
      
    protected:
      float d_average;
      float d_zrmsd;    // <zero value> rmsd
      int d_count;
      
    };
      
// torsional angle property class
  class TorsionProperty : public Property
    {
    public:
      static const float pi = 3.14159265359;
      
      TorsionProperty(gcore::System &sys);
      virtual ~TorsionProperty();
      
      void parse(int count, std::string arguments[]);
      
      virtual float calc();
      virtual std::string average();
      virtual std::string toString();
      virtual std::string toTitle();
      
      struct Exception: public gromos::Exception
      {
	  Exception(const string &what): gromos::Exception("TorsionProperty", what) {}
      };
      
    protected:
      float d_average;
      float d_zrmsd;
      int   d_count;
    };    

  inline float Property::getValue()
    {
      return d_value;
    }
  
  inline float Property::getZValue()
    {
      return d_zvalue;
    }

}

#endif













