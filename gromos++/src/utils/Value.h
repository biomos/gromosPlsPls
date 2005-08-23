/**
 * @file Value.h
 * Value
 */

/* 	$Id$	 */

#ifndef INCLUDED_UTILS_VALUE
#define INCLUDED_UTILS_VALUE

#ifndef INCLUDED_VECTOR
#include <vector>
#define INCLUDED_VECTOR
#endif

#ifndef INCLUDED_STRING
#include <string>
#define INCLUDED_STRING
#endif

#ifndef INCLUDED_UTILS_VECTORSPECIFIER
#include "VectorSpecifier.h"
#endif

#ifndef INCLUDED_GROMOS_EXCEPTION
#include "../gromos/Exception.h"
#endif

namespace utils
{
  enum value_enum { val_scalar = 0, val_vector = 1, val_vecspec = 2 };
  
  
  class Value
  {
  public:
    Value() 
      : d_value_type(val_scalar), d_scalar(0.0), d_vec(0.0) {}
    Value(double d)
      : d_value_type(val_scalar), d_scalar(d), d_vec(0.0) {}
    Value(gmath::Vec v)
      : d_value_type(val_vector), d_scalar(0.0), d_vec(v) {}
    Value(VectorSpecifier const & v)
      : d_value_type(val_vecspec), d_scalar(0.0), d_vecspec(v) {}
    Value(Value const &v)
      : d_value_type(v.d_value_type), d_scalar(v.d_scalar),
	d_vec(v.d_vec), d_vecspec(v.d_vecspec) {}

    gmath::Vec vec()const
    {
      switch(d_value_type){
	case val_scalar:
	  throw Exception("no vector value");
	case val_vector:
	  return d_vec;
	case val_vecspec:
	  return d_vecspec();
      }
      throw Exception("wrong value type");
    }
    double scalar()const
    {
      if (d_value_type != val_scalar)
	throw Exception("no scalar value");
      return d_scalar;
    }
    
    value_enum type()const
    {
      return d_value_type;
    }

    // maybe they should rather return the Value?
    double operator=(double d)
    {
      d_value_type = val_scalar;
      d_scalar = d;
      return d;
    }
    
    gmath::Vec operator=(gmath::Vec const & v)
    {
      d_value_type = val_vector;
      d_vec = v;
      return d_vec;
    }

    gmath::Vec operator=(VectorSpecifier const & v)
    {
      d_value_type = val_vecspec;
      d_vecspec = v;
      return d_vecspec();
    }
    
    std::string toString()const
    {
      std::ostringstream os;
      switch(d_value_type){
	case val_scalar:
	  os << d_scalar;
	  break;
	case val_vector:
	  os << "(" << d_vec[0] << ", "
	     << d_vec[1] << ", "
	     << d_vec[2] << ")";
	  break;
	case val_vecspec:
	  os << "(" << d_vecspec()[0] << ", "
	     << d_vecspec()[1] << ", "
	     << d_vecspec()[2] << ")";
	  break;
      }
      return os.str();
    }

    void parse(std::string s)
    {
      std::string::size_type bra = s.find('(');
      if (bra != std::string::npos){
	// assume it's a vector specifier
	d_value_type = val_vecspec;
	d_vecspec.setSpecifier(s);
      }
      else{
	std::istringstream is(s);
	if (!(is >> d_scalar))
	  throw Exception("could not parse value");
      }
    }

    struct Exception: public gromos::Exception{
      Exception(const std::string &what): 
	gromos::Exception("Value", what){}
    };
    
  protected:
    value_enum d_value_type;
    double d_scalar;
    gmath::Vec d_vec;
    VectorSpecifier d_vecspec;
  };
  
  inline std::ostream & operator<<(std::ostream & os, Value const & v)
  {
    os << v.toString();
    return os;
  }

  // some math
  inline Value operator+(Value const &v1, Value const &v2)
  {
    switch(v1.type()){
      case val_scalar:
	{
	  if (v2.type() != val_scalar)
	    throw Value::Exception("adding different types");
	  return Value(v1.scalar() + v2.scalar());
	}
      case val_vector:
      case val_vecspec:
	{
	  if (v2.type() == val_vector ||
	      v2.type() == val_vecspec)
	    return Value(v1.vec() + v2.vec());
	  else
	    throw Value::Exception("adding different types");
	}
    }
    throw Value::Exception("wrong value type");
  }
  inline Value operator-(Value const &v1, Value const &v2)
  {
    switch(v1.type()){
      case val_scalar:
	{
	  if (v2.type() != val_scalar)
	    throw Value::Exception("adding different types");
	  return Value(v1.scalar() - v2.scalar());
	}
      case val_vector:
      case val_vecspec:
	{
	  if (v2.type() == val_vector ||
	      v2.type() == val_vecspec)
	    return Value(v1.vec() - v2.vec());
	  else
	    throw Value::Exception("adding different types");
	}
    }
    throw Value::Exception("wrong value type");
  }
  inline Value operator*(Value const &v1, Value const &v2)
  {
    if (v1.type() == val_scalar){
      if (v2.type() == val_scalar)
	return Value(v1.scalar() * v2.scalar());
      else return Value(v1.scalar() * v2.vec());
    }
    else{
      if (v2.type() == val_scalar)
	return Value(v1.vec() * v2.scalar());
      else throw Value::Exception("use dot / cross to do vector products");
    }
  }
  inline Value operator/(Value const &v1, Value const &v2)
  {
    if (v2.type() != val_scalar)
      throw Value::Exception("division by vector");

    if (v1.type() == val_scalar)
      Value(v1.scalar() / v2.scalar());
    else
      return Value(v1.vec() / v2.scalar());

    throw Value::Exception("wrong value type");
  }
  inline Value abs(Value const &v)
  {
    if (v.type() == val_scalar)
      return Value(fabs(v.scalar()));
    return Value(v.vec().abs());
  }
  inline Value abs2(Value const &v)
  {
    if (v.type() == val_scalar)
      return Value(v.scalar() * v.scalar());
    return Value(v.vec().abs2());
  }
  inline Value dot(Value const &v1, Value const &v2)
  {
    if (v1.type() == val_scalar ||
	v2.type() == val_scalar)
      throw Value::Exception("dot product requires two vectors as arguments");

    return Value(v1.vec().dot(v2.vec()));
  }
  inline Value cross(Value const &v1, Value const &v2)
  {
    if (v1.type() == val_scalar ||
	v2.type() == val_scalar)
      throw Value::Exception("cross product requires two vectors as arguments");

    return Value(v1.vec().cross(v2.vec()));
  }
}

#endif

