// utils_CheckTopo.h

// Class that runs some basic checks on a molecule topology
#ifndef INCLUDED_VECTOR
#include <vector>
#define INCLUDED_VECTOR
#endif
#ifndef INCLUDED_STRING
#include <string>
#define INCLUDED_STRING
#endif
#ifndef INCLUDED_MAP
#include <map>
#define INCLUDED_MAP
#endif
#ifndef INCLUDED_GMATH_EXPRESSION
#include "../gmath/Expression.h"
#define INCLUDED_GMATH_EXPRESSION
#endif
#ifndef INCLUDED_GIO_GINSTREAM
#include "../gio/Ginstream.h"
#define INCLUDED_GIO_GINSTREAM
#endif

namespace utils{
  /**
   * Class EnergyTraj
   * A class that contains all energy terms from (free) energy trajectories
   * 
   * The EnergyTraj class stores the energy terms which can be accessed via 
   * (user defined) names. It also allows the user to define how properties 
   * should be calculated from other values, defined through a 
   * gmath::Expression 
   *
   * @class EnergyTraj
   * @author B.C. Oostenbrink
   * @ingroup utils
   */
  class EnergyTraj
    {
      /**
       * A vector of doubles that will contain all energy trajectory data
       *
       * The (free) energy elements will be in a fixed order, bounded by the
       * constants that are defined below. First, come the energy elements 
       * ENER[], then the special energy elements ENERES[], the volume pressure
       * temperature elements VOLPRT[], the lambda value RLAM, the free energy 
       * derivatives FREN[]. After this come the energy group elements ENERLJ[]
       * ENERCL[] ENERRF[] and ENERRC[]. These come in the end because the
       * number of energy groups is non-constant.
       */
      std::vector<double> d_data;
      /**
       * A vector of expressions that define how a property should be computed
       *
       * The properties that need to be computed get internally a negative 
       * index.
       */
      std::vector<gmath::Expression> d_e;
      /**
       * A vector of dependencies for the expressions
       *
       * A property that needs to be calculated according to an expression
       * is always based on a number of known properties. The indices of 
       * the properties that are needed to calculate a property are stored
       * in this vector of integers.
       */
      std::vector<std::vector <int> > d_dep;
      /**
       * A vector that stores the result of the properties that need to be
       * calculated so that a second access does not recalculate.
       */
      std::vector<double> d_calc;
      /**
       * A vector of booleans that keeps track of which properties have
       * been calculated already and which one need to be recalculated.
       */
      std::vector<bool> d_recalc;
      /**
       * A type definition for the map
       */
      typedef std::map<std::string, int>::value_type MP;
      /** 
       * A map that links the name of a property to its index.
       * Since we do not know the number of energy groups in advance, these
       * energies are not put in the map, but are handled seperately in the 
       * index function
       */
      std::map<std::string, int> d_map;
      /**
       * a counter for the number of frames that have been read.
       */
      int d_en_frame;
      int d_fr_frame;
      /**
       * The number of energy groups that has been read in 
       */
      int d_negr;
      /** 
       * last index of the energy elements (ENER[1..22])
       */
      int num_ener;
      /**
       * last index of the special energies (ENERES[1..6])
       */
      int num_eneres;
      /**
       * last index of the VOLPRT block (VOLPRT[1..20])
       */
      int num_volprt;
      /**
       * last index of the RLAM (RLAM)
       */
      int num_rlam;
      /**
       * last index of the derivatives of the energies with respect
       * to lambda (FREN[1..22]
       */
      int num_fren;
      /**
       * number of ENER[] entries in the free energy file
       */
      int num_enerf;
      /**
       * A constant that will be returned if one tries to access
       * an unknown element
       */
      int unknownvariable;
      /**
       * A string to keep track of the first kind of block in the energy 
       * trajectory
       */
      std::string en_first;
      /**
       * A string to keep track of the first kind of block in the free
       * energy trajectory
       */
      std::string fr_first;
      
    public:
      /**
       * Constructor
       *
       * This constructor constructs a standard Energy Trajectory wich has
       * all the elements that are defined on page III-56 of the GROMOS96 
       * manual.
       */
      EnergyTraj();
      /**
       * Constructor
       *
       * This constructor allows the user to prepare the EnergyTraj for
       * energy trajectories with non-standard numbers of elements.
       * @param vector<int> numbers A vector of integers that contains the
       *                            number of energy elements in the following
       *                            order<br>
       *                            numbers[0]: number of ENER[] elements<br>
       *                            numbers[1]: number of ENERES[] elements<br>
       *                            numbers[2]: number of VOLPRT[] elements<br>
       *                            numbers[3]: number of lambda parameters<br>
       *                            numbers[4]: number of FREN[] elements<br>
       *                            numbers[5]: number of ENER[] elements
       *                                      in the free energy trajectory<br>
       */
      EnergyTraj(std::vector<int> numbers);
      /**
       * Accessor to get the i-th element of the data-set. The first 'num_fren'
       * elements refer to the data in the energy files
       */
      double operator[](int i);
      /**
       * Accessor to get the element in the data set that is referred to with
       * the name s
       *
       * The name of an element can be the name according to the (free) energy
       * files (page III-56 of the GROMOS manual) or a user defined name (see
       * function addKnown
       */
      double operator[](std::string s);
      /**
       * function to read in one frame from the energy trajectory.
       */
      void read_energy_frame(gio::Ginstream& is);
      /**
       * function to read in one frame from the free energy trajectory.
       */
      void read_free_energy_frame(gio::Ginstream& is);
      /**
       * function to teach the class about a new property
       *
       * A new property can be either a direct mapping of an existing known 
       * property, or it can be an expression, that can be calculated from 
       * previously known properties.
       * @param string s The name of the new property
       * @param string v An expression (or name) of the existing properties
       *                 that it refers to.
       */
      void addKnown(std::string s, std::string v);
      /**
       * A function that returns the index-number of the internal data array
       * belonging to the property with name s
       */
      int index(std::string s);
      /**
       * A function that writes out the definitions and connections of all 
       * known properties. All known properties are listed in alphabetical 
       * order, followed by their definition. This can be just a mapping to the
       * data array, or it can be a user defined property for which the
       * definition and the dependencies are also listed.
       */
      void write_map(std::ostream& os = std::cout);
    private:
      /**
       * A function that backtracks the first name that can be found,
       * which refers to the index number i
       */
      std::string back_index(int i);
      /**
       * A function to Tokenize a string into a vector of strings
       */
      void Tokenize(const std::string& str,
		    std::vector<std::string>& tokens,
		    const std::string& delimiters = " ");
      /**
       * A function to initialize the free energy trajectory. Set some
       * variables and learn about the standard names;
       */
      void init(std::vector<int> numbers);
    };
}



      
  
