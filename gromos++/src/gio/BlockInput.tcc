/**
 * @file blockinput.tcc
 * defines blockinput functions.
 */

inline 
std::istream& 
gio::getline(
	    std::istream& is, 
	    std::string& s, 
	    const char& sep,
	    const char& comm
	    )
{
  unsigned short int ii;

  while (is.good()) {
    std::getline(is, s, sep);
    ii = std::find(s.begin(), s.end(), comm) - s.begin();
    if (ii == s.size()) break; // no comment
    else if (!ii) continue;    // comment on first position
    else s.erase(s.begin() + ii, s.end());
  }
  
  return is;
}

inline 
std::istream& 
gio::getblock(
	     std::istream& is, 
	     std::vector<std::string>& b, 
	     const std::string& sep
  )
{

  if (!b.size())
    b.push_back("");
  std::vector<std::string>::iterator dest = b.begin();

  while (1) {

    if (dest == b.end()) {
      b.push_back("");
      dest = b.end() - 1;
    }       
    
    getline(is, *dest);

    if (*dest == sep || !is.good())
      break;

    ++dest;
  }

  ++dest;
  b.erase(dest, b.end());

  return is;
}

inline
std::string& 
gio::concatenate(
		std::vector<std::string>::const_iterator begin,
		std::vector<std::string>::const_iterator end,
		std::string& s,
		const char& sep
		)
{
  s.clear();
  while (begin != end) {
    s += *begin;
    s += sep;
    begin++;
  }
  
  return s;
}
