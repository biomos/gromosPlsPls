// gromos_exception.h

#ifndef INCLUDED_GROMOS_EXCEPTION
#define INCLUDED_GROMOS_EXCEPTION

#ifndef INCLUDED_EXCEPTION
#include <exception>
#define INCLUDED_EXCEPTION
#endif

#ifndef INCLUDED_STRING
#include <string>
#define INCLUDED_STRING
#endif

using namespace std;

namespace gromos
{
  class Exception: public std::exception{
    string d_class;
    string d_what;
  public:
    Exception(const string &from,
	      const string &mess) throw():
      d_class(from),
      d_what(mess)
      {}
  virtual  const char *what() const throw(){
      return (d_class + ": " + d_what).c_str();
    }
  virtual  ~Exception() throw(){}
  };
}



#endif
