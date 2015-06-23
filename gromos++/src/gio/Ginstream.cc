/**
 * @file Ginstream.cc
 * basic input stream class definition.
 */

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <algorithm>    // std::remove_if
#include "../gromos/Exception.h"
#include "Ginstream.h"
#include <cstdio>
#include "gzstream.h"

template<class size_type>
inline std::basic_string<size_type>&
trim_right( std::basic_string<size_type>& str )
{
    return( str = str.substr( 0, str.find_last_not_of( ' ' ) + 1 ) );
}


template<class size_type>
inline std::basic_string<size_type>&
trim( std::basic_string<size_type>& str )
{
  if (str.find_first_not_of( ' ' ) == std::string::npos) return (str = "");
  return( trim_right( str ) );
}


gio::Ginstream::Ginstream(const std::string &s, std::ios::openmode mode)
  :_is(NULL)
{
  open(s, mode);
  _has_version = false;
}

gio::Ginstream::Ginstream(gio::Ginstream &gin)
{
  _is = gin._is;
  _title = gin._title;
  _name = gin._name;
  _has_version = gin._has_version;
  _version = gin._version;
}

void gio::Ginstream::open(const std::string s, std::ios::openmode mode)
{
  igzstream *gis = new igzstream(s.c_str(), mode);
  if (gis == NULL){
    throw gromos::Exception("Ginstream", "could not create a std::ifstream( " + s + ")");
  }
  if(!gis->good()){
    throw gromos::Exception("Ginstream", "Could not open file '" + s + "'");
  }  
  if(!gis->is_open()){
    throw gromos::Exception("Ginstream", "could not open file '" + s + "'");
  }
  
  stream(*gis, true);
  if (!_has_version) {
    // We tried to read ENEVERSION, but apparently, there was none.
    // Our stream isn't usable anymore - we started reading the 
    // next block - probably TIMESTEP.
    // Rewind operations (seekg() and the like) don't seem to work - 
    // so let's recreate it. (Yeah, not THAT elegant...)
    gis->close();
    delete gis;
    igzstream *gis = new igzstream(s.c_str(), mode);
    if (gis == NULL){
      throw gromos::Exception("Ginstream", "could not create a std::ifstream( " + s + ")");
    }
    if(!gis->good()){
      throw gromos::Exception("Ginstream", "Could not open file '" + s + "'");
    }  
    if(!gis->is_open()){
      throw gromos::Exception("Ginstream", "could not open file '" + s + "'");
    }
    
    stream(*gis, false);
  }
  _name=s;
}


void gio::Ginstream::close()
{
  // std::cerr << "closing _is" << std::endl;

  std::ifstream * ifs = dynamic_cast<std::ifstream *>(_is);
  if (ifs != NULL)
    ifs->close();
  else{
    igzstream * igzs = dynamic_cast<igzstream *>(_is);
    igzs->close();
  }

  // _is->close();
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

std::istream& gio::Ginstream::getline(std::string& s, 
					     const char& sep,
					     const char& comm){
  //unsigned short int ii;
  std::string::size_type ii;
  
  while (_is->good()) {
    std::getline(*_is, s, sep);
    //ii = std::find(s.begin(), s.end(), comm) - s.begin();
    ii=s.find(comm,0);
    
    

    if(!s.size()) continue;                 // empty line
    else if(ii == std::string::npos) break; // no comment
    else if (!ii) continue;                 // comment on first position
    else {
      s.erase(s.begin() + ii, s.end());
      if (!trim_right(s).size()) continue;  // line with comment only
      break;
    }
    
  }
  s = trim(s);

  return *_is;
}

std::istream& gio::Ginstream::getblock(std::vector<std::string>& b, 
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
    
    if (dest->find(sep) ==0)
      break;
   
    if (!_is->good()) 
      throw gromos::Exception("getblock", "error reading block."+*b.begin());
    
    
    ++dest;
  }
  
  ++dest;
  b.erase(dest, b.end());

  //std::cout << "B: " << b.size() << std::endl;
  // for (int i=0; i < (int) b.size(); ++i) b[i] = trim(b[i]);
  
  return *_is;
}

std::istream& gio::Ginstream::skipblock(const std::string& sep)
{
  std::string s;

  while (1) {

    getline(s);

    if(_is->eof()) {
      break;
    }
    
    if (s.find(sep) == 0)
      break;

    if (!_is->good()) 
      throw gromos::Exception("getblock", "error skipping block." + s);
  }
  
  return *_is;
}


std::string& gio::concatenate(
		std::vector<std::string>::const_iterator begin,
		std::vector<std::string>::const_iterator end,
		std::string& s,
		const char& sep) {
  //s.clear();
  s="";
  while (begin != end) {
    s += *begin;
    s += sep;
    begin++;
  }
  
  return s;
}

void gio::Ginstream::readVersion(){  
  std::vector<std::string> _b;

  getblock(_b);
  if (_b[0] != "ENEVERSION") {
    _has_version = false;
    return;
  }
  // else, save new version and set bool switch
  _version = gio::concatenate(_b.begin() + 1, _b.end() - 1, _version);
  // we're ignoring any whitespaces - less error-prone
  _version.erase( std::remove_if( _version.begin(), _version.end(), ::isspace ), _version.end() );

  _has_version = true;
}

bool gio::Ginstream::has_version(){
  return _has_version;
}

std::string gio::Ginstream::version() {
  if (_has_version) {
    return _version;
  } else {
    return "";
  }
}

