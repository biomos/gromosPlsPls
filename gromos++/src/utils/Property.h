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
  class Property;
  std::ostream &operator<<(std::ostream &os, Property &p);

  /**
   * Class Property 
   * Purpose: General base class for properties to be calculated
   * during an analysis run.
   *
   * Description:
   * This class defines (and implements) the general methods that a
   * property calculated during an analysis run should posses.
   * Specialized derived classes for (ie) bond lengths, angles or
   * dihedral angles exist.
   *
   * @class Property
   * @version Jul 31 15:00 2002
   * @author M. Christen
   * @ingroup utils
   * @sa utils::PropertyContainer
   * @sa utils::DistanceProperty
   * @sa utils::AngleProperty
   * @sa utils::TorsionProperty
   */

  class Property
  {
  public:
    /**
     * Maximum number of arguments.
     */
    static const int MAXARGUMENTS = 10;

    /**
     * Standard constructor. Not used for predefined properties.
     * (Because they need all a reference to the system to be able
     * to calculate their value)
     */
    // Property(); // for user defined properties that do not need a system
    /**
     * Standard constructor.
     * Takes a reference to the system, so that the properties can
     * calculate their value (ie by looking up the atom positions).
     */
    Property(gcore::System &sys);
    /**
     * Destructor.
     */
    virtual ~Property();

    // accessor
    std::string type() {return d_type;}

    /**
     * Return the value of the property.
     * In retrospect, i should have written the whole thing as a
     * template, so that anything could be returned here.
     * Maybe sometime somebody wants to change that.
     */
    double getValue();
    /**
     * Returns the ideal value of the property (if one has been specified).
     */
    double getZValue();
    /**
     * As most of the properties i can think of have something to do with
     * atoms and molecules, i define these vectors in the base class.
     * This is also done in order to be able to write one general
     * arguments parsing function in the base class.
     */
    // std::vector<int> atoms();
    /**
     * Vector of the molecule mol[i] corresponding to atom[i].
     * @sa utils::Property::atoms
     */
    // std::vector<int> mols();

    AtomSpecifier & atoms();

    // methods
    
    /**
     * Calculates the value of the property.
     * Override in derived classes.
     */
    virtual double calc();
    /**
     * After a call to utils::Property::calc() checkBounds() can check,
     * whether the value of the property lies within the specified range.
     * (This only works if calc() calculates one single value and stores
     * it into d_value)
     */
    virtual std::string checkBounds();
    /**
     * Write a title string specifying the property.
     */
    virtual std::string toTitle();
    /**
     * Returns the value in string representation.
     */
    virtual std::string toString();
    /**
     * Returns the average value over all calls to calc.
     */
    virtual std::string average();

    /**
     * Returns the type of the interaction from the
     * topology.
     */
    virtual int getTopologyType(gcore::System const &sys);

    /**
     * @struct Exception
     * Property exception
     */
    struct Exception: public gromos::Exception
    {
      /**
       * Constructor.
       */
      Exception(const std::string &what): 
	gromos::Exception("Property", what) {}
    };
    
  protected:
    /**
     * Parse the command line property specification.
     * This is the standard implementation. It knows about
     * molecules, atoms, zero value and boundaries.
     */
    virtual void parse(int count, std::string arguments[]);
    /**
     * Helper method to parse the atom part of the property specifier.
     * The argument is in AtomSpecifier format.
     */
    void parseAtoms(std::string atoms);
    /**
     * Helper method to parse the atoms belonging to one molecule.
     */
    // void _parseAtomsHelper(std::string substring, int &mol);

    /**
     * find the corresponding forcefield type of the property.
     * needs to be overwritten for the specific properties.
     */
    virtual int findTopologyType(gcore::MoleculeTopology const &mol_topo);
    
  protected:
    // member variables
    /**
     * Number of required arguments. Used in the standard parse function.
     * Set in the constructor.
     */
    int REQUIREDARGUMENTS; 
    /**
     * The property type (in string representation).
     */
    std::string d_type;

    /**
     * The bounds for 'boundary violation check'.
     * @sa utils::Property::checkBounds
     */
    double d_ubound, d_lbound;
    /**
     * The zero/equilibrium/standard value of the property.
     */
    double d_zvalue;
    /**
     * Stores the calculated value. This is used for subsequent toString
     * calls.
     * If in a user defined property, d_value is not used, those functions
     * must be overridden.
     */
    double d_value;
    
    /**
     * The atoms belonging to this property.
     */
    // std::vector<int> d_atom;
    /**
     * The molecule, the atoms belong to.
     */
    // std::vector<int> d_mol;
    AtomSpecifier d_atom;
    
    /**
     * Reference of the system.
     */
    gcore::System *d_sys;

  };

  
  /**
   * Class DistanceProperty
   * Purpose: Implements a distance property.
   *
   * Description:
   * This class implements a distance property. It is derived from the 
   * Property class.
   *
   * @class DistanceProperty
   * @version Wed Jul 31 2002
   * @author gromos++ development team
   * @sa utils::Property
   */

  class DistanceProperty : public Property
  {    
  public:
    /**
     * Constructor.
     */
    DistanceProperty(gcore::System &sys);
    /**
     * Destructor.
     */
    virtual ~DistanceProperty();
    /**
     * Parses the arguments. Is overridden to check the input, but basically
     * just calls Property::parse.
     * @sa utils::Property::parse
     */
    virtual void parse(int count, std::string arguments[]);
    /**
     * Calculates the distance.
     */
    virtual double calc();
    /**
     * Averages over the calculated values.
     * Not implemented ?! (in class Property?!)
     */
    virtual std::string average();
    /**
     * Get the average value of all calc() calls.
     */
    // virtual std::string toString();
    /**
     * Get a title string.
     */
    // virtual std::string toTitle();
    /**
     * @struct Exception
     * DistanceProperty exception.
     */
    struct Exception: public gromos::Exception
    {
      /**
       * Constructor.
       */
      Exception(const std::string &what): 
	gromos::Exception("DistanceProperty", what) {}
    };

  protected:

    /**
     * find the corresponding forcefield type of the property.
     * needs to be overwritten for the specific properties.
     */
    virtual int findTopologyType(gcore::MoleculeTopology const &mol_topo);

    /**
     * The average value.
     */
    double d_average;
    /**
     * The rmsd from the zero/initial/ideal value.
     */
    double d_zrmsd;
    /**
     * The number of times that calc() has been called.
     */
    int d_count;
    
  };

  /**
   * Class AngleProperty
   * Purpose: Implements an angle property.
   *
   * Description:
   * Implements an angle property. It is derived from the Property class.
   *
   * @class AngleProperty
   * @version Wed Jul 31 2002
   * @author gromos++ development team
   * @sa utils::Property
   */

  class AngleProperty : public Property
    {
    public:
      /**
       * Constructor.
       */
      AngleProperty(gcore::System &sys);
      /**
       * Destructor.
       */
      virtual ~AngleProperty();
      /**
       * Parse and check property specifier (given in arguments).
       * Calls Property::parse and checks the arguments.
       */
      virtual void parse(int count, std::string arguments[]);
      /**
       * Calculate the angle between the given atoms.
       */
      virtual double calc();
      /**
       * Calculate the average of all values calculated so far.
       */
      virtual std::string average();
      /**
       * Print the calculated value.
       */
      // virtual std::string toString();
      /**
       * Get a title string.
       */
      // virtual std::string toTitle();
      /**
       * @struct Exception
       * AngleProperty exception.
       */
      struct Exception: public gromos::Exception
      {
	/**
	 * Constructor.
	 */
	Exception(const std::string &what): 
	  gromos::Exception("AngleProperty", what) {}
      };
      
    protected:
      /**
       * find the corresponding forcefield type of the property.
       * needs to be overwritten for the specific properties.
       */
      virtual int findTopologyType(gcore::MoleculeTopology const &mol_topo);

      /**
       * The average value.
       */
      double d_average;
      /**
       * The rmsd from the ideal value.
       */
      double d_zrmsd;    // <zero value> rmsd
      /**
       * Number of times calc() has been called.
       */
      int d_count;
      
    };
      
  /**
   * Class TorsionProperty
   * Purpose: Implements a torsion property.
   *
   * Description:
   * This class implements a torsion property. It is derived from the 
   * Property class.
   *
   * @class TorsionProperty
   * @version Wed Jul 31 2002
   * @author gromos++ development team
   * @sa utils::Property
   */

  class TorsionProperty : public Property
    {
    public:
      /**
       * Constructor.
       */
      TorsionProperty(gcore::System &sys);
      /**
       * Destructor.
       */
      virtual ~TorsionProperty();
      /**
       * Parse the arguments. Calls Property::parse.
       */
      virtual void parse(int count, std::string arguments[]);
      /**
       * Calculate the torsional angle.
       */
      virtual double calc();
      /**
       * Get the average of the calculated values.
       */
      virtual std::string average();
      /**
       * Get the value as string.
       */
      // virtual std::string toString();
      /**
       * Get a title string.
       */
      // virtual std::string toTitle();
      /**
       * @struct Exception
       * TorsionProperty exception.
       */
      struct Exception: public gromos::Exception
      {
	/**
	 * Constructor.
	 */
	  Exception(const std::string &what): 
	    gromos::Exception("TorsionProperty", what) {}
      };
      
    protected:
      /**
       * find the corresponding forcefield type of the property.
       * needs to be overwritten for the specific properties.
       */
      virtual int findTopologyType(gcore::MoleculeTopology const &mol_topo);

      /**
       * The average value.
       */
      double d_average;
      /**
       * The rmsd from the ideal value.
       */
      double d_zrmsd;
      /**
       * Number of times calc() has been called.
       */
      int   d_count;
    };    

  /**
   * Class OrderProperty
   * Purpose: Implements an order property.
   *
   * Description:
   * Implements an order property
   * (angle between axis and vector specified
   * by two atoms).
   * It is derived from the Property class.
   *
   * @class OrderProperty
   * @version Jan 16 2004
   * @author gromos++ development team
   * @sa utils::Property
   */

  class OrderProperty : public Property
    {
    public:
      /**
       * Constructor.
       */
      OrderProperty(gcore::System &sys);
      /**
       * Destructor.
       */
      virtual ~OrderProperty();
      /**
       * Parse and check property specifier (given in arguments).
       * Calls Property::parse and checks the arguments.
       */
      virtual void parse(int count, std::string arguments[]);
      /**
       * Calculate the angle between the given atoms.
       */
      virtual double calc();
      /**
       * Calculate the average of all values calculated so far.
       */
      virtual std::string average();
      /**
       * Print the calculated value.
       */
      // virtual std::string toString();
      /**
       * Get a title string.
       */
      virtual std::string toTitle();
      /**
       * @struct Exception
       * AngleProperty exception.
       */
      struct Exception: public gromos::Exception
      {
	/**
	 * Constructor.
	 */
	Exception(const std::string &what): 
	  gromos::Exception("OrderProperty", what) {}
      };
      
    protected:
      /**
       * find the corresponding forcefield type of the property.
       * needs to be overwritten for the specific properties.
       */
      virtual int findTopologyType(gcore::MoleculeTopology const &mol_topo);

      /**
       * The average value.
       */
      double d_average;
      /**
       * The rmsd from the ideal value.
       */
      double d_zrmsd;    // <zero value> rmsd
      /**
       * Number of times calc() has been called.
       */
      int d_count;
      /**
       * axis with respect to which the angle is calculated.
       */
      gmath::Vec d_axis;
      
    };

  /**
   * Class OrderProperty
   * Purpose: Implements an order property.
   *
   * Description:
   * Implements an order property
   * (angle between axis and vector specified
   * by two atoms).
   * It is derived from the Property class.
   *
   * @class OrderProperty
   * @version Jan 16 2004
   * @author gromos++ development team
   * @sa utils::Property
   */

  class OrderParamProperty : public Property
    {
    public:
      /**
       * Constructor.
       */
      OrderParamProperty(gcore::System &sys);
      /**
       * Destructor.
       */
      virtual ~OrderParamProperty();
      /**
       * Parse and check property specifier (given in arguments).
       * Calls Property::parse and checks the arguments.
       */
      virtual void parse(int count, std::string arguments[]);
      /**
       * Calculate the angle between the given atoms.
       */
      virtual double calc();
      /**
       * Calculate the average of all values calculated so far.
       */
      virtual std::string average();
      /**
       * Print the calculated value.
       */
      // virtual std::string toString();
      /**
       * Get a title string.
       */
      virtual std::string toTitle();
      /**
       * @struct Exception
       * AngleProperty exception.
       */
      struct Exception: public gromos::Exception
      {
	/**
	 * Constructor.
	 */
	Exception(const std::string &what): 
	  gromos::Exception("OrderProperty", what) {}
      };
      
    protected:
      /**
       * find the corresponding forcefield type of the property.
       * needs to be overwritten for the specific properties.
       */
      virtual int findTopologyType(gcore::MoleculeTopology const &mol_topo);

      /**
       * The average value.
       */
      double d_average;
      /**
       * The rmsd from the ideal value.
       */
      double d_zrmsd;    // <zero value> rmsd
      /**
       * Number of times calc() has been called.
       */
      int d_count;
      /**
       * axis with respect to which the angle is calculated.
       */
      gmath::Vec d_axis;
      
    };

  inline double Property::getValue()
    {
      return d_value;
    }
  
  inline double Property::getZValue()
    {
      return d_zvalue;
    }

  /*
  inline std::vector<int> Property::atoms()
    {
      return d_atom;
    }
  
  inline std::vector<int> Property::mols()
    {
      return d_mol;
    }
  */
  inline utils::AtomSpecifier & Property::atoms()
  {
    return d_atom;
  }

  /**
   * Class PseudoRotationProperty
   * Purpose: Implements a property that can calculate a pseudo rotation.
   *
   * Description:
   * This class implements a pseudo rotation. It is derived from the 
   * Property class. The pseudo rotation (\Delta / 2 ) is defined according  
   * to:
   * Altona, C; Geise, HJ; Romers C; Tetrahedron 24 13-32 (1968) 
   * With the addition that if \theta_0 < 0, then 180 is added to the value
   * For a (DNA)-sugar furanose ring this means that if one specifies the 
   * atoms for this property as C1',C2',C3',C4',O4' you can determine the
   * puckering, see:
   * Altona, C; Sundaralingam, M; JACS 94 8205-8212 (1972)
   *
   * @class PseudoRotationProperty
   * @version Fri Apr 23 2004
   * @author gromos++ development team
   * @sa utils::Property
   */

  class PseudoRotationProperty : public Property
    {
    public:
      /**
       * Constructor.
       */
      PseudoRotationProperty(gcore::System &sys);
      /**
       * Destructor.
       */
      virtual ~PseudoRotationProperty();
      /**
       * Parse the arguments. Calls Property::parse.
       */
      virtual void parse(int count, std::string arguments[]);
      /**
       * Calculate the torsional angle.
       */
      virtual double calc();
      /**
       * Get the average of the calculated values.
       */
      virtual std::string average();
      /**
       * Get the value as string.
       */
      // virtual std::string toString();
      /**
       * Get a title string.
       */
      // virtual std::string toTitle();
      /**
       * @struct Exception
       * TorsionProperty exception.
       */
      struct Exception: public gromos::Exception
      {
	/**
	 * Constructor.
	 */
	  Exception(const std::string &what): 
	    gromos::Exception("PseudoRotationProperty", what) {}
      };
      
    protected:
      /**
       * find the corresponding forcefield type of the property.
       * needs to be overwritten for the specific properties.
       */
      virtual int findTopologyType(gcore::MoleculeTopology const &mol_topo);

      /**
       * The average value.
       */
      double d_average;
      /**
       * The rmsd from the ideal value.
       */
      double d_zrmsd;
      /**
       * Number of times calc() has been called.
       */
      int   d_count;
      /**
       * Function to calculate a torsion for four atoms
       */
      double _calcDihedral(int const a, int const b, int const c, int const d);
      /**
       * A constant that is needed every time
       * Should be sin(36) + sin(72);
       */
      double d_sin36sin72;
  };      

  /**
   * Class PuckerAmplitudeProperty
   * Purpose: Implements a property that can calculate a the amplitude of a 
   * pucker rotation.
   *
   * Description:
   * This class implements a pucker amplitude. It is derived from the 
   * PseudoRotationProperty class. The amplitude is defined according  
   * to:
   * Altona, C; Geise, HJ; Romers C; Tetrahedron 24 13-32 (1968) 
   * see also:
   * Altona, C; Sundaralingam, M; JACS 94 8205-8212 (1972)
   *
   * @class PseudoRotationProperty
   * @version Fri Apr 23 2004
   * @author gromos++ development team
   * @sa utils::PseudoRotationProperty
   */

  class PuckerAmplitudeProperty : public PseudoRotationProperty
    {
    public:
      /**
       * Constructor.
       */
      PuckerAmplitudeProperty(gcore::System &sys);
      /**
       * Destructor.
       */
      virtual ~PuckerAmplitudeProperty();
      /**
       * Parse the arguments. Calls Property::parse.
       */
      virtual double calc();
      /**
       * Get a title string.
       */
      // virtual std::string toTitle();
  };    
}

#endif

