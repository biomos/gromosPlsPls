// StringTokenizer.h

#ifndef INCLUDED_GIO_STRINGTOKENIZER
#define INCLUCED_GIO_STRINGTOKENIZER

#ifndef INCLUDED_STRING
#include <string>
#define INCLUDED_STRING
#endif
#ifndef INCLUDED_VECTOR
#include <vector>
#define INCLUDED_VECTOR
#endif
#ifndef INCLUDED_GROMOS_EXCEPTION
#include "../gromos/Exception.h"
#endif


namespace gio{
  
  /**
   * Class StringTokenizer
   * Takes a string and returns tokens as vector of strings
   * separated by a delimiter.
   *
   * @class StringTokenizer
   * @ingroup gio
   * @author M.A. Kastenholz
   */
  class StringTokenizer{
    // not implemented
    StringTokenizer(const StringTokenizer&);
    StringTokenizer &operator=(const StringTokenizer&);
    std::string to_tokenize;
    std::string delimiter;

    public:
    // Constructors
    StringTokenizer(const std::string &str, const std::string& delim = " ");
    ~StringTokenizer() {};

    // Methods
    void setdelimiter(const std::string delim);
    void setstring(const std::string str);
    std::vector<std::string> tokenize();
    
    //Exceptions
    struct Exception: public gromos::Exception{
      Exception(const std::string& what_arg) : gromos::Exception("StringTokenizer", what_arg){}
    };
  };
}
#endif
