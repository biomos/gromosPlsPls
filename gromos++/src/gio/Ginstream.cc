/**
 * @file Ginstream.cc
 * basic input stream class definition.
 */

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include "../gromos/Exception.h"
#include "Ginstream.h"
#include <stdio.h>

gio::Ginstream::Ginstream(const std::string &s, std::ios::openmode mode)
:_is(NULL){
  open(s, mode);
}

gio::Ginstream::Ginstream(Ginstream::Ginstream &gin)
{
  _is = gin._is;
  _title = gin._title;
  _name = gin._name;
}

void gio::Ginstream::open(const std::string s, std::ios::openmode mode)
{
  // std::cerr << "trying to open a file" << std::endl;
  std::ifstream *is = new std::ifstream(s.c_str(), mode);
  if (is == NULL){
    throw gromos::Exception("Ginstream", "could not create a std::ifstream( " + s + ")");
  }
  // std::cerr << "created a new ifstream" << std::endl;
  if(!is->good()){
    // std::cerr << "errno: " << errno << std::endl;
    throw gromos::Exception("Ginstream", "Could not open file "+s);
  }  
  // std::cerr << "file opened succesfully" << std::endl;

  stream(*is);
  //std::cout << "and assigned it to gin" << std::endl;
  
  _name=s;
  //std::cout << "it worked" << std::endl;
  
}


void gio::Ginstream::close()
{
  // std::cerr << "closing _is" << std::endl;
  _is->close();
  _title="";
  _name="";
  delete _is;
}

void gio::Ginstream::readTitle() {

  std::vector<std::string> _b;

  getblock(_b);
  if (_b[0] != "TITLE")
    throw gromos::Exception("Ginstream", 
			    "TITLE block expected. Found: " + _b[0]);
  _title = gio::concatenate(_b.begin() + 1, _b.end() - 1, _title);
}
std::string gio::Ginstream::name()
{
  return _name;
}

std::string gio::Ginstream::title() 
{
  return _title;
}

inline std::ifstream& gio::Ginstream::getline(std::string& s, 
					     const char& sep,
					     const char& comm){
  unsigned short int ii;

  while (_is->good()) {
    std::getline(*_is, s, sep);
    ii = std::find(s.begin(), s.end(), comm) - s.begin();
    if(!s.size()) continue;         // empty line
    else if (ii == s.size()) break; // no comment
    else if (!ii) continue;         // comment on first position
    else s.erase(s.begin() + ii, s.end());
  }
  
  return *_is;
}

inline std::ifstream& gio::Ginstream::getblock(std::vector<std::string>& b, 
					      const std::string& sep)
{
  if (!b.size())
    b.push_back("");
  std::vector<std::string>::iterator dest = b.begin();

  while (1) {

    if (dest == b.end()) {
      b.push_back("");
      dest = b.end() - 1;
    }       
    getline(*dest);

    if(_is->eof()) {
      --dest;
      break;
    }
    
    if (*dest == sep)
      break;
   
    if (!_is->good()) 
      throw gromos::Exception("getblock", "error reading block."+*b.begin());
    
    
    ++dest;
  }
  
  ++dest;
  b.erase(dest, b.end());
  
  return *_is;
}


inline std::string& gio::concatenate(
		std::vector<std::string>::const_iterator begin,
		std::vector<std::string>::const_iterator end,
		std::string& s,
		const char& sep) {
  s.clear();
  while (begin != end) {
    s += *begin;
    s += sep;
    begin++;
  }
  
  return s;
}
