// gmath_Matrix.cc

#include "Matrix.h"
#include "Vec.h"
#include <new>
#include <cassert>
#include "../gsl/matrix/gsl_matrix.h"
#include "../gsl/linalg/gsl_linalg.h"
#include "../gsl/header/gsl_math.h"
#include "../gsl/eigen/gsl_eigen.h"

using namespace std;

namespace gmath{

Matrix::Matrix(int rows, int columns){
  d_rows=rows;
  d_columns=columns;
  d_val=new (double* [d_rows]);
  for (int i=0;i<d_rows;++i){
    d_val[i]=new (double [d_columns]);
  }
}

Matrix::Matrix(int rows, int columns, double value){
  d_rows=rows;
  d_columns=columns;
  d_val=new (double* [d_rows]);
  for (int i=0;i<d_rows;++i){
    d_val[i]=new (double [d_columns]);
    for(int j=0;j<d_columns;++j){
      d_val[i][j]=value;
    }
  }
}

Matrix::Matrix(const Matrix &mat){
  d_rows=mat.d_rows;
  d_columns = mat.d_columns;
  d_val=new (double* [d_rows]);
  for (int i=0;i<d_rows;++i){
    d_val[i]=new (double [d_columns]);
    for(int j=0;j<d_columns;++j){
      d_val[i][j]=mat.d_val[i][j];
    }
  }
}

Matrix::Matrix(const Vec &v, const Vec &w){
  d_rows=3;
  d_columns=3;
  d_val=new (double* [d_rows]);
  for (int i=0;i<d_rows;++i){
    d_val[i]=new (double [d_columns]);
    for(int j=0;j<d_columns;++j){
      d_val[i][j]=v[i]*w[j];
    }
  }
}

Matrix::Matrix(const Vec &u, const Vec &v, const Vec &w){
  d_rows=3;
  d_columns=3;
  d_val=new (double* [d_rows]);
  for (int i=0;i<d_rows;++i){
    d_val[i]=new (double [d_columns]);
    d_val[i][0]=u[i];
    d_val[i][1]=v[i];
    d_val[i][2]=w[i];
  }
}


Matrix &Matrix::operator=(const Matrix &mat){
  if (this != &mat){
    this->~Matrix();
    new(this) Matrix(mat);
  }
  return *this;
}

Matrix::~Matrix(){
  for(int i=0;i<d_rows;++i)
    delete d_val[i];
  delete d_val;
}

Matrix Matrix::luDecomp(){
  assert(d_rows==d_columns);
  Matrix mat(*this);

  gsl_matrix * gsl_mat = gsl_matrix_alloc (mat.rows(), mat.columns());
  gsl_matrix_set_zero (gsl_mat);
  for (int i=0; i < mat.rows(); ++i){
   for (int j=0; j < mat.columns(); j++){
     gsl_matrix_set (gsl_mat, i , j , mat(i,j));
   }
  }  

  int s;
 gsl_permutation * p = gsl_permutation_alloc (mat.rows());
 gsl_linalg_LU_decomp (gsl_mat, p, &s);

 Matrix ret(mat.rows(), mat.columns());
 for (int i=0; i < mat.rows(); ++i){
   for (int j=0; j < mat.columns(); j++){
     ret(i,j)=gsl_matrix_get(gsl_mat, i, j);
   }
  }  

 return (ret);
}

double Matrix::det()const{
 assert(d_rows==d_columns);
 Matrix tmp(*this);
 
 gsl_matrix * gsl_mat = gsl_matrix_alloc (tmp.rows(), tmp.columns());
 gsl_matrix_set_zero (gsl_mat);
  for (int i=0; i < tmp.rows(); ++i){
   for (int j=0; j < tmp.columns(); j++){
     gsl_matrix_set (gsl_mat, i , j , tmp(i,j));
   }
  }  
 gsl_permutation * p = gsl_permutation_alloc (tmp.rows());
 int s; 
 gsl_linalg_LU_decomp (gsl_mat, p, &s);


 double d = gsl_linalg_LU_det(gsl_mat, s);
 return d;
}


Matrix Matrix::diagonaliseSymmetric(double *eigenValues){
  assert(d_rows==d_columns);
  Matrix mat(*this);
  
   gsl_matrix * gsl_mat = gsl_matrix_alloc (mat.rows(), mat.columns());
   gsl_matrix_set_zero (gsl_mat);
  for (int i=0; i < mat.rows(); ++i){
   for (int j=0; j < mat.columns(); j++){
     gsl_matrix_set (gsl_mat, i , j , mat(i,j));
   }
  }  

  gsl_vector *eval = gsl_vector_alloc (mat.rows());
  gsl_matrix *evec = gsl_matrix_alloc (mat.rows(), mat.columns());

  gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (mat.rows());

  gsl_eigen_symmv (gsl_mat, eval, evec, w);

  gsl_eigen_symmv_free(w);

  gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_VAL_DESC);

 Matrix ret(mat.rows(), mat.columns());
 for (int i=0; i < mat.rows(); ++i){
  eigenValues[i] = gsl_vector_get(eval, i);
   for (int j=0; j < mat.columns(); j++){
     ret(i,j)=gsl_matrix_get(evec, i, j);
   }
  }  

  return ret;
}

Vec operator*(const Matrix &m, const Vec &v){
  assert(m.rows()==3&&m.columns()==3);
  Vec temp;
  for (int i=0;i<3;++i)
    for(int j=0;j<3;++j)
      temp[i]+=m(i,j)*v[j];
  return temp;
}

  Matrix operator*(const Matrix &m1, const Matrix &m2){
    assert(m1.rows()==m2.columns()&&m1.columns()==m2.rows());
    Matrix temp(m1.rows(),m1.columns(),0);
    for(int i=0;i<m1.rows();++i)
      for(int j=0;j<m1.columns();++j)
	for(int k=0;k<m1.columns();++k)
	  temp(i,j)+=m1(i,k)*m2(k,j);
    return temp;
  }

}
