// DihedralType.h
#ifndef INCLUDED_DIHEDRALTYPE
#define INCLUDED_DIHEDRALTYPE

namespace gcore{

class DihedralType
{
  double d_fc;
  double d_pd;
  int d_np;
 public:
  DihedralType(double fc=0, double pd=0, int np=0):
    d_fc(fc), d_pd(pd), d_np(np){}
  DihedralType(const DihedralType& b):d_fc(b.d_fc), d_pd(b.d_pd), d_np(b.d_np){}
  DihedralType &operator=(const DihedralType &b);
  double pd()const{return d_pd;}
  int np()const{return d_np;}
  double fc()const{return d_fc;}
};

}
#endif
