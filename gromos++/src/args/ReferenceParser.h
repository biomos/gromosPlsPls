// args_ReferenceParser.h

#ifndef INCLUDED_ARGS_REFERENCEPARSER
#define INCLUDED_ARGS_REFERENCEPARSER

#include <vector.h>
#include "../fit/Reference.h"

namespace gcore{
  class System;
}

namespace args{
  class Arguments;
}

namespace fit{
  class Reference;
}

namespace args{
  class ReferenceParser{
    
    vector<int> myMols;
    gcore::System &mySys;
    const args::Arguments &myArgs;
    fit::Reference &myRef;
    bool myAdded;

  public:

    // constructor
    ReferenceParser::ReferenceParser(
      gcore::System &sys,
      const args::Arguments &args,
      fit::Reference &ref);

    void add_ref();
    void getMolecules();
    void classes();
    void atoms();
  };

}
#endif
