// StringTokenizer.cc
#include <string>
#include <vector>
#include "StringTokenizer.h"

using gio::StringTokenizer;

  StringTokenizer::StringTokenizer(const std::string &str, const std::string& delim) {
    to_tokenize = str;
    delimiter = delim;
  };
  


  void StringTokenizer::setdelimiter(const std::string delim) {
   delimiter = delim;
  }
  void StringTokenizer::setstring(const std::string str) {
   to_tokenize = str;
  }

  std::vector<std::string> StringTokenizer::tokenize() {
     std::vector< std::string> tokens;
   // Skip delimiters at beginning.
  std::string::size_type lastPos = to_tokenize.find_first_not_of(delimiter, 0);
  // Find first "non-delimiter".
  std::string::size_type pos     = to_tokenize.find_first_of(delimiter, lastPos);

  while (std::string::npos != pos || std::string::npos != lastPos)
    {
      // Found a token, add it to the vector.
      tokens.push_back(to_tokenize.substr(lastPos, pos - lastPos));
      // Skip delimiters.  Note the "not_of"
      lastPos = to_tokenize.find_first_not_of(delimiter, pos);
      // Find next "non-delimiter"
      pos = to_tokenize.find_first_of(delimiter, lastPos);
    }

   return tokens;
 
  }


