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

  class Arguments: public multimap<std::string,std::string>{
    Arguments_i *d_this;
    // not implemented
    Arguments();
    Arguments(const Arguments &);
    Arguments &operator=(const Arguments &);
  public:
    // Constructors
  Arguments(int argc, char **args, int nknown, char **known, const std::string &usage); 

   ~Arguments();

  // Methods
  int check(const std::string &str, int num_args=0)const;
  // This has to be in to fix a bug in gcc Solaris 2.6 (?)
//   const const_iterator begin()const
 //   {return this->multimap<string,string>::begin();}

  friend std::istream &operator>>(std::istream &is, Arguments &args);
  
  // Accessors
  const std::string &operator[](const std::string &str)const;

  // Exceptions
  struct Exception: public gromos::Exception{
    Exception(const std::string &str): gromos::Exception("Usage", str) {}
  };


};

}

#endif







