#ifndef INCLUDED_ARGS_ARGUMENTS
#define INCLUDED_ARGS_ARGUMENTS

#ifndef INCLUDED_MAP
#include <map>
#define INCLUDED_MAP
#endif

#ifndef INCLUDED_GROMOS_EXCEPTION
#include "../gromos/Exception.h"
#endif

namespace args{

  class Arguments_i;

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
 * @sa args::BoundaryParser
 * @sa args::GatherParser
 */

  class Arguments: public multimap<std::string,std::string>{
    Arguments_i *d_this;
    // not implemented
    Arguments();
    Arguments(const Arguments &);
    Arguments &operator=(const Arguments &);
  public:
/**
 * Arguments constructor.
 * Details.
 */
  Arguments(int argc, char **args, int nknown, char **known, const std::string &usage); 
/**
 * Arguments destructor.
 * Details.
 */
   ~Arguments();

/** 
 * Constructs the class and returns a member pointer to a Boundary.
 * Method parse parses the input from args["pbc"].
 * @param &str Takes a std::string as argument.
 * @param num_args Integer of the number of arguments.
 * @return check Integer to check for failure.
 * Details.
 */
  int check(const std::string &str, int num_args=0)const;
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
    Exception(const std::string &str): gromos::Exception("Usage", str) {}
  };


};

}

#endif







