#include <iostream>
#include <vector>
#include <string>

#include "BlockInput.h"

using namespace std;

int main()
{
  vector<string> block;

  try {
    while (gio::getblock(cin, block)) {
      for (unsigned int ii = 0; ii < block.size(); ii++) 
        cout << block[ii] << endl;
    }
  }
  catch (std::exception& e) {
    cerr << e.what() << endl;
    return 1;
  }

  return 0;
}
