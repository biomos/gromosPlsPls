#ifndef INCLUDED_ATOMPAIR
#define INCLUDED_ATOMPAIR

namespace gcore{

class AtomPair{
  int d_a[2];
  // not implemented
  AtomPair();
 public:
  AtomPair(int a, int b);
  AtomPair(const AtomPair &);
  int operator[](int i)const{return d_a[i];}
};

int operator<(const AtomPair &a, const AtomPair &b);
}
#endif
