// gmath_Matrix.h

#ifndef INCLUDED_GMATH_MATRIX
#define INCLUDED_GMATH_MATRIX

#ifndef INCLUDED_GROMOS_EXCEPTION
#include "../gromos/Exception.h"
#endif

#ifndef INCLUDED_STRING
#include <string>
#define INCLUDED_STRING
#endif

#ifndef INCLUDED_CASSERT
#include <cassert>
#define INCLUDED_CASSERT
#endif

namespace gmath{

  class Vec;

  class Matrix{
    double **d_val;
    int d_rows, d_columns;
    
    // not implemented
    Matrix();
  public:
    Matrix(int rows, int columns);
    Matrix(int rows, int columns, double value);
    Matrix(const Matrix &);
    Matrix(const Vec &v, const Vec &w);
    Matrix(const Vec &u, const Vec &v, const Vec &w);
    
    // "dyadic product" m_ij = v_i*w_j
    ~Matrix();


    // Methods
    Matrix &operator=(const Matrix &);
    void luDecomp(int *index, double *d);
      // LU Decomposition according to Num. Recipes p.46
      // indx[0..n-1] is an output vector that records the row permutation 
      // effected by the partial pivoting; d is output as   
      // 1 depending on whether the number of row interchanges 
      // was even or odd, respectively.
      // The input Matrix will be modified!
    void diagonaliseSymmetric(double *eigenValues);
      // diagonalise a symmetric Matrix and return eigenvalues.
    double det()const;

    // operators
    Matrix operator-()const;
    Matrix &operator+=(const Matrix &mat);
    Matrix &operator-=(const Matrix &mat);
    Matrix &operator*=(double d);

    // Accessors
    double operator()(int i, int j)const;
    double &operator()(int i, int j);
    int rows()const;
    int columns()const;

    // Exception
    struct Exception: public gromos::Exception{
      Exception(const string &what): gromos::Exception("Matrix",what){}
    };

  };

  // inline functions & free operators

  inline Matrix operator+(const Matrix &mat1, const Matrix &mat2){
    assert(mat1.rows()==mat2.rows() && mat1.columns()==mat2.columns());
    Matrix temp(mat1);
    temp+=mat2;
    return temp;
  }

  inline Matrix operator-(const Matrix &mat1, const Matrix &mat2){
    assert(mat1.rows()==mat2.rows() && mat1.columns()==mat2.columns());
    Matrix temp(mat1);
    temp-=mat2;
    return temp;
  }

  inline Matrix operator*(double d, const Matrix &m){
    Matrix temp(m);
    temp*=d;
    return temp;
  }

  Matrix operator*(const Matrix &m1, const Matrix &m2);

  Vec operator*(const Matrix &m, const Vec &v);

  inline Matrix &Matrix::operator+=(const Matrix &mat){
    assert(rows()==mat.rows() && columns()==mat.columns());
    for(int i=0; i<rows();++i)
      for(int j=0;j<columns();++j)
	d_val[i][j]+=mat.d_val[i][j];
    return *this;
  }

  inline Matrix &Matrix::operator-=(const Matrix &mat){
    assert(rows()==mat.rows() && columns()==mat.columns());
    for(int i=0; i<rows();++i)
      for(int j=0;j<columns();++j)
	d_val[i][j]-=mat.d_val[i][j];
    return *this;
  }

  inline Matrix &Matrix::operator*=(double d){
    for(int i=0; i<rows();++i)
      for(int j=0;j<columns();++j)
	d_val[i][j]*=d;
    return *this;
  }

  inline Matrix Matrix::operator-()const{
    Matrix temp(*this);
    for(int i=0; i<rows();++i)
      for(int j=0;j<columns();++j)
	temp.d_val[i][j]=-temp.d_val[i][j];
    return temp;
  }
    
  inline double Matrix::operator()(int i, int j)const{
    return d_val[i][j];
  }

  inline double &Matrix::operator()(int i, int j){
    return d_val[i][j];
  }

  inline int Matrix::rows()const{
    return d_rows;
  }
  
  inline int Matrix::columns()const{
    return d_columns;
  }

}


#endif
