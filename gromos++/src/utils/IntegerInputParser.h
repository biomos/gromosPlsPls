#include <set>
#include <string>

// A class that can take comma-separated identifiers, ranges, or a combination (e.g. 1,4,7-9)
// to parse user input in Integer format for various uses (e.g. selection of energy groups or similar)


namespace utils
{
  class IntegerInputParser: public std::set<int>
  {
  public:
   // Method to add the Integer numbers in a string (e.g. from user input)
   // maxnum is the highest number allowed (smallest is hardcoded to zero)
    void addSpecifier(std::string const s, int maxnum);
  protected:
    void parse(std::string const s, int maxnum);
  };
}
