// utils_Rmsd.h

#ifndef INCLUDED_UTILS_RMSD
#define INCLUDED_UTILS_RMSD

namespace gcore{
  class System;
}

namespace fit{
  class Reference;
}

namespace utils{

  class Rmsd{
    const fit::Reference *d_ref;
    // not implemented
    Rmsd();
    Rmsd(const Rmsd &);
    Rmsd &operator=(const Rmsd&);
  public:
    Rmsd(const fit::Reference *);
    ~Rmsd(){}
    double rmsd(const gcore::System &)const;
  };
}

#endif
