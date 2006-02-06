
#ifndef INCLUDED_UTILS_PARSE
#define INCLUDED_UTILS_PARSE

namespace utils
{

  /**
   * find matching bracket
   */
  std::string::size_type find_matching_bracket
  (
   std::string s,
   char bra='(',
   std::string::size_type it=0
   );

  /**
   * find character, taking care of brackets in string
   * ie: looking for ';' in va(prop1;prop2);prop3
   * will find second ';'
   */
  std::string::size_type find_par
  (
   std::string s,
   char c=';',
   std::string::size_type it=0,
   std::string bra = "([{<",
   std::string ket = ">}])"
   );

  /**
   * find character contained in c, taking care of brackets in string
   * ie: looking for ';' in va(prop1;prop2);prop3
   * will find second ';'
   */
  std::string::size_type find_par
  (
   std::string s,
   std::string c,
   std::string::size_type it=0,
   std::string bra = "([{<",
   std::string ket = ">}])"
   );

  /**
   * parse a range into an
   * index array
   */
  void parse_range(std::string s, std::vector<int> & range, int x=-1);
  
}


#endif
