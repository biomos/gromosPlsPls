// gromos_exception.h

#ifndef INCLUDED_GROMOS_EXCEPTION
#define INCLUDED_GROMOS_EXCEPTION

#include <exception>
#include <string>

namespace gromos
{
  /**
   * Class Exception
   * Gromos++ exception class.
   *
   * Recently the use of this class has been questioned. Now that I write 
   * the description I sort of seem to agree.
   *
   * @class Exception
   * @author R. Buergi
   * @ingroup gromos
   */
  class Exception: public std::exception{
    std::string d_mesg;
  public:
    /**
     * Exception constructor
     * @param from A description of which program / routine 
     *             throws the description
     * @param mess The error message
     */
    Exception(const std::string from,
	      const std::string mess) throw():
      d_mesg(from + ": " + mess) 
      {}
  virtual  const char *what() const throw(){
      return d_mesg.c_str();
    }
  virtual  ~Exception() throw(){}
  };
}



#endif
