// gcore_Exclusion.h

#ifndef INCLUDED_GCORE_EXCLUSION
#define INCLUDED_GCORE_EXCLUSION

namespace gcore{

class Exclusion_i;

class Exclusion{
  Exclusion_i *d_this;
 public:
  Exclusion();
  Exclusion(const Exclusion&);
  ~Exclusion();

  // Methods
  Exclusion &operator=(const Exclusion &);

  void erase(int i);
  // remove atom i from exclusion list

  void insert(int i);
  // add Atom i to exclusion list

  int size() const;
  // number of exclusions

  int atom(int i) const;
  // get atom number being number i in exclusion list
};

}

#endif
