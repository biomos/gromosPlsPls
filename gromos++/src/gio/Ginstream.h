// gio_Ginstream.h

#ifndef INCLUDED_GIO_GINSTREAM
#define INCLUDED_GIO_GINSTREAM

#ifndef INCLUDED_FSTREAM
#include <fstream>
#define INCLUDED_FSTREAM
#endif
#ifndef INCLUDED_GROMOS_EXCEPTION
#include "../gromos/Exception.h"
#endif

namespace gmath{
  class Vec;
}

using gmath::Vec;

namespace gio{

class Ginstream: public std::ifstream {
  static const int MAXCHAR;
  std::string d_title;
  std::string d_name;
  void readTitle();

  // not implemented
  Ginstream(const Ginstream& g);
  Ginstream& operator=(const Ginstream& g);
public:
  // Constructors
  Ginstream();
  // Open file and read in TITLE block.
  Ginstream(const std::string &name, ios::openmode mode=ios_base::in);
 
  ~Ginstream();
  //open file and read title
  void open(const char *name, ios::openmode mode=ios_base::in);
  void open(const std::string &name, ios::openmode mode=ios_base::in);
  void close();

  // methods
  // skips comments and checks for EOF.
  int eof();
  
  Ginstream &getline(std::string &str, int max=MAXCHAR);
  Ginstream &operator>>(std::string &s) ;
  Ginstream &operator>>(int &c) ;
  Ginstream &operator>>(double &d) ;
  Ginstream &operator>>(float &d) ;
  Ginstream &operator>>(Vec &v);
  Ginstream &operator>>(char* c);

  // Checks if the next std::string is str.
  bool check(const std::string &str="END");

  /// Accessors
  const std::string &name() const;
  const std::string &title() const;
  

  // Exceptions
  struct Exception: public gromos::Exception{
    Exception(const std::string &what):gromos::Exception("Ginstream", what) {}
  };
};
}
    
#endif
