/**
 * @file VectorSpecifier.h
 * VectorSpecifier methods
 */

// Class that contains a vector
// specified by cartesian coordinates,
// polar coordinates or
// 

#ifndef INCLUDED_UTILS_VECTORSPECIFIER
#define INCLUDED_UTILS_VECTORSPECIFIER

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

// minimal complete headers
#include "../utils/AtomSpecifier.h"

namespace gmath
{
  class Vec;
}

namespace bound
{
  class Boundary;
}

namespace utils
{
  /**
   * @class VectorSpecifier
   * @author C. Oostenbrink, M. Christen
   * @ingroup utils
   *
   * Class VectorSpecifier
   */
  class VectorSpecifier{

    AtomSpecifier d_atomspec;
    gmath::Vec d_vec;
    bound::Boundary * d_pbc;
    
  public: 
    // Constructors
    /** 
     * VectorSpecifier standard constructor
     */
    VectorSpecifier() : d_atomspec(), d_vec()
    {}

    /**
     * VectorSpecifier Constructor
     * @param sys The VectorSpecifier needs to know about the system. It 
     *            does not know about any atoms yet.
     */
    VectorSpecifier(gcore::System &sys, bound::Boundary * pbc)
      : d_atomspec(sys), d_vec(), d_pbc(pbc)
    {}
    
    /**
     * VectorSpecifier Constructor
     * @param sys The VectorSpecifier needs to know about the system.
     * @param s   A string of the correct format. Usually this is provided
     *            by the user, so it is assumed to start numbering at 1
     */
    VectorSpecifier(gcore::System &sys, bound::Boundary * pbc, std::string s);

    /**
     * copy constructor!
     */
    VectorSpecifier(VectorSpecifier const & vs);

    /**
     * VectorSpecifier Destructor
     */
    ~VectorSpecifier();

    /**
     * Method to set the system the atom specifier is referring to
     * @param sys the system
     */
    void setSystem(gcore::System &sys)
    {
      d_atomspec.setSystem(sys);
    }
    
    /**
     * set the boundary condition object.
     */
    void setBoundary(bound::Boundary *pbc)
    {
      d_pbc = pbc;
    }
    
    /**
     * Method to add parse a string to the VectorSpecifier.
     * @param s Is assumed to be user-specified, 
     * with numbering starting at 1
     */
    int setSpecifier(std::string s);

    /**
     * Member operator = copies one VectorSpecifier into the other
     */
    VectorSpecifier &operator=(const VectorSpecifier &vs);

    /**
     * Accessor, returns the vector
     */    
    gmath::Vec const & operator()();
    /**
     * set vector to zero, empty atom specifier
     */
    void clear();
    /**
     * Accessor, returns a pointer to the system on which the VectorSpecifier
     * is based
     */
    gcore::System & sys();
    /**
     * boundary condition object accessor
     */
    bound::Boundary * pbc();
    
    /**
     * Method, returns a vector of strings that
     * would reproduce the
     * VectorSpecifier
     */
    std::string toString();
    /**
     * @struct Exception
     * Throws an exception if something is wrong
     */
    struct Exception: public gromos::Exception{
      /**
       * @exception If called says VectorSpecifier, followed by the argument
       * @param what The string that is thrown
       */
      Exception(const std::string &what): 
	gromos::Exception("VectorSpecifier", what){}
    };
  protected:
    //Internal function
    /**
     * Parse the arguments string into the VectorSpecifier
     */
    void parse(std::string s);

    void parse_cartesian(std::string s);
    void parse_polar(std::string s);
    void parse_atom(std::string s);

    /**
     * find matching closing bracket for a given opening bracket
     * at position it in the string s
     */
    std::string::size_type find_matching_bracket(std::string s, 
						 char bra='(',
						 std::string::size_type it=0);
  };
  
}

#endif

