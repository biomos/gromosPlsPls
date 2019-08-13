#include "IntegerInputParser.h"

void utils::IntegerInputParser::parse(string const s, int maxnum)
{
  if(s=="ALL" || s=="all"){
    for(int i=0; i<maxnum; i++) insert(i+1);
    return;
  }
  std::string::size_type iterator;
  if((iterator=s.find(',')) != std::string::npos){
    parse(s.substr(0,iterator), maxnum);
    parse(s.substr(iterator+1, std::string::npos), maxnum);
  }
  else{
    istringstream is;
    int rangeBegin, rangeEnd;
    if((iterator=s.find('-')) != std::string::npos){
      is.str(s.substr(0,iterator));
      if(!(is >> rangeBegin))
	throw gromos::Exception("IntegerInputParser", 
				"Invalid begin of range "+ s);
      is.clear();
      is.str(s.substr(iterator+1, std::string::npos));
      if(!(is >> rangeEnd))
	throw gromos::Exception("IntegerInputParser", 
				"Invalid end of range " + s);
      for(int i=rangeBegin; i<= rangeEnd; ++i){
	if(i> maxnum)
	  throw gromos::Exception("IntegerInputParser",
				  "Requested number too high: "+s);
	insert(i);
      }
    }
    else{
      is.str(s);
      if(!(is >> rangeBegin))
	throw gromos::Exception("IntegerInputParser", 
				"Invalid number specified "+ s);
      if(rangeBegin > maxnum)
	throw gromos::Exception("IntegerInputParser",
				"Requested number too high: "+s);
      insert(rangeBegin);
    }
  }
}

void utils::IntegerInputParser::addSpecifier(string const s, int maxnum)
{
  parse(s, maxnum);
}
