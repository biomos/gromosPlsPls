// args_BoundaryParser.h

#ifndef INCLUDED_ARGS_BOUNDARYPARSER
#define INCLUDED_ARGS_BOUNDARYPARSER

namespace gcore{
  class System;
}

namespace bound{
  class Boundary;
}

namespace args{
  class Arguments;

  class BoundaryParser{
    
    // not implemented
    BoundaryParser(const BoundaryParser&);
    BoundaryParser &operator=(const BoundaryParser&);

  public:
    BoundaryParser(){}
    ~BoundaryParser(){}
    static bound::Boundary *boundary(gcore::System &sys, const Arguments &args);
  };

}






#endif
