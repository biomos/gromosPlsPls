#ifndef INCLUDED_ARGS_ARGUMENTS
#define INCLUDED_ARGS_ARGUMENTS

#ifndef INCLUDED_MAP
#include <map>
#define INCLUDED_MAP
#endif

#ifndef INCLUDED_SET
#include <set>
#define INCLUDED_SET
#endif

#ifndef INCLUDED_GROMOS_EXCEPTION
#include "../gromos/Exception.h"
#endif

#include <stdlib.h>

namespace args{
  class Arguments_i;
 /**
  * Class Argument_List
  * Purpose: add an easy addition interface for the arguments set
  *
  * @class Argument_List
  * @version $Date: Mon Sep 16 14:27:41 MEST 2008
  * @author  N. Schmid
  * @ingroup args
  * @sa args::Arguments
  */
  class Argument_List {
  public:
    std::set<std::string> known;
    /**
     * adds an argument (c string) to the argument list
     */
    Argument_List & operator<<(const char* arg) {
      known.insert(std::string(arg));
      return *this;
    }
  };

/**
 * Class Arguments
 * Purpose: Parse arguments from the command line or an input file.
 * 
 *
 * Description:
 * This class is used to parse arguments from the command line or an input file .
 *
 * 
 * @class Arguments
 * @version $Date: Mon Jul 15 14:17:41 MEST 2002
 * @author  R. Buergi
 * @author  M. Kastenholz
 * @ingroup args
 * @sa args::BoundaryParser
 * @sa args::GatherParser
 */

  class Arguments: public std::multimap<std::string,std::string>{
    Arguments_i *d_this;
    // not implemented
    Arguments();
    Arguments(const Arguments &);
    Arguments &operator=(const Arguments &);
  public:
    // global variables needed for @i08 etc. flags
    static bool inG96;
    static bool outG96;
/**
 * Arguments constructor.
 * Details.
 */
  Arguments(int argc, char **args, const Argument_List & known, const std::string &usage); 
/**
 * Arguments destructor.
 * Details.
 */
   ~Arguments();

/** 
 * Checks for whether the argument string had num_arg arguments
 * @param &str Takes a std::string as argument.
 * @param num_args Integer of the number of arguments.
 * @return check Integer to check for failure.
 */
  int check(const std::string &str, int num_args=0)const;

  /**
   * Returns the number of arguments that follow string
   * @param &str Takes a std::string as argument.
   * @return the number of arguments found for this argument. Returns -1 
   * if string was not found at all in the argument list.
   */
  int count(const std::string &str)const;
  
  // This has to be in to fix a bug in gcc Solaris 2.6 (?)
//   const const_iterator begin()const
 //   {return this->multimap<string,string>::begin();}

/** 
 * Member operator >>.
 * Details.
 */
  friend std::istream &operator>>(std::istream &is, Arguments &args);
  
/** 
 * Member operator [] used to access the members.
 * Details.
 */
  const std::string &operator[](const std::string &str)const;

/**
 * @struct Exception. 
 * Throws the Usage string if invoked.
 */
  struct Exception: public gromos::Exception{
/**
 * @exception const std::string &str throws Usage string if invoked.
 */
    Exception(const std::string &str): gromos::Exception("# Usage", str) {}
  };

};

}

#endif







