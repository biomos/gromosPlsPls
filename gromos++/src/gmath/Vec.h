// gmath_Vec

#ifndef INCLUDED_GMATH_VEC
#define INCLUDED_GMATH_VEC

#ifndef INCLUDED_CASSERT
#include <cassert>
#define INCLUDED_CASSERT
#endif

#ifndef INCLUDED_CMATH
#include <cmath>
#define INCLUDED_CMATH
#endif

namespace gmath{

class Vec{
  double d_v[3];

 public:
  Vec(double i=0, double j=0, double k=0);
  Vec(const Vec &v);
  ~Vec(){}
  // Methods
  Vec &operator=(const Vec &v);
  Vec operator-()const;
  Vec &operator+=(const Vec &v);
  Vec &operator-=(const Vec &v);
  Vec &operator*=(const Vec &v);
  Vec &operator*=(double d);
  Vec &operator/=(double d);
  bool operator==(const Vec &v)const;

  // Cross (outer) product of two gmath::Vectors
  Vec cross(const Vec &v)const;
  // Dot (inner) product of two gmath::Vectors
  double dot(const Vec &v)const;
  // returns the normalized gmath::Vector, i.e. the Vector divided by its absolute value
  Vec normalize()const;
  // The dot product of a gmath::Vector with itself (norm)
  double abs2()const;
  // The sqrt of the norm of the gmath::Vector (absolute value)
  double abs()const;


  // Accessor
  double operator[](int i)const;
  double &operator[](int i);
};


// Class Methods

 inline Vec::Vec(double i, double j, double k){
   d_v[0]=i;
   d_v[1]=j;
   d_v[2]=k;
 }
 
 inline Vec::Vec(const Vec &v){
   d_v[0]=v.d_v[0];
   d_v[1]=v.d_v[1];
   d_v[2]=v.d_v[2];
 }

 inline Vec Vec::operator-()const{
   return Vec(-d_v[0], -d_v[1], -d_v[2]);
 }
 
 inline Vec &Vec::operator=(const Vec &v){
   if(this != &v){
     d_v[0]=v.d_v[0];
     d_v[1]=v.d_v[1];
     d_v[2]=v.d_v[2];
   }   
   return *this;
 }

 inline Vec &Vec::operator+=(const Vec &v){
   d_v[0]+=v.d_v[0];
   d_v[1]+=v.d_v[1];
   d_v[2]+=v.d_v[2];
   return *this;
 }

 inline Vec &Vec::operator*=(const Vec &v){
   d_v[0]*=v.d_v[0];
   d_v[1]*=v.d_v[1];
   d_v[2]*=v.d_v[2];
   return *this;
 }

  
 inline Vec &Vec::operator-=(const Vec &v){
   d_v[0]-=v.d_v[0];
   d_v[1]-=v.d_v[1];
   d_v[2]-=v.d_v[2];
   return *this;
 }

 inline Vec &Vec::operator*=(double d){
   d_v[0]*=d;
   d_v[1]*=d;
   d_v[2]*=d;
   return *this;
 } 

 inline Vec &Vec::operator/=(double d){
   d_v[0]/=d;
   d_v[1]/=d;
   d_v[2]/=d;
   return *this;
 } 

 inline Vec Vec::cross(const Vec &v) const{
   return Vec(d_v[1]*v.d_v[2] - d_v[2]*v.d_v[1],
	      d_v[2]*v.d_v[0] - d_v[0]*v.d_v[2],
	      d_v[0]*v.d_v[1] - d_v[1]*v.d_v[0]);
 }



 inline double Vec::dot(const Vec &v) const{
   return d_v[0]*v.d_v[0] + d_v[1]*v.d_v[1] + d_v[2]*v.d_v[2];
 }

 inline double Vec::abs2()const{
   return this->dot(*this);
 }

 inline double Vec::abs()const{
   return std::sqrt(this->abs2());
}

 inline Vec Vec::normalize()const{
   double d = std::sqrt(this->abs2());
   return Vec(d_v[0]/d, d_v[1]/d, d_v[2]/d);
 }

 inline double Vec::operator[](int i)const{
   assert(i<3);
   return d_v[i];
 }

 inline double &Vec::operator[](int i){
   assert(i<3);
   return d_v[i];
 }

 inline Vec operator+(const Vec &a, const Vec &b){
   Vec v(a);
   v+=b;
   return v;
 }

 inline Vec operator-(const Vec &a, const Vec &b){
   Vec v(a);
   v-=b;
   return v;
 }

 inline Vec operator*(const Vec &a, const Vec &b){
   Vec v(a);
   v*=b;
   return v;
 }

 inline Vec operator*(double d, const Vec &b){
   Vec v(b);
   v*=d;
   return v;
 }

 inline Vec operator*(const Vec &a, double d){
   Vec v(a);
   v*=d;
   return v;
 }

 inline Vec operator/(const Vec &a, double d){
   Vec v(a);
   v/=d;
   return v;
 }

 
}
#endif





