// gmath_Matrix.cc

#include "Matrix.h"
#include "Vec.h"
#include <new>
#include <cassert>

using namespace std;

namespace gmath{

// used for LU decomposition according to Numerical Recipes
static void ludcmp(double **a, int n, int *indx, double *d);
static const double TINY = 1.0e-80;
// Put a symmetric Matrix into tridiagonal form according to N.R.
static void tred2(double **a, int n, double d[], double e[]);
// Diagonalise a tridiagonal matrix
static void tqli(double d[], double e[], int n, double **z);
// sort eigenvectors and eigenvalues
static void sort(Matrix *mat, double* eigenvalues);
// swap columns
static void swapColumns(Matrix *mat, int i, int j);

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

void Matrix::luDecomp(int *index, double *d){
  assert(d_rows==d_columns);
  ludcmp(d_val,d_rows,index,d);
  // is it sensible to subsitute TINY back?
  //  for(int i=0;i<d_rows;++i)
  //    for(int j=0;j<d_columns;++j)
  //      if((*this)(i,j)==TINY)(*this)(i,j)=0;
}

double Matrix::det()const{
  assert(d_rows==d_columns);
  Matrix tmp(*this);
  int *index=new (int[d_rows]);
  double d;
  tmp.luDecomp(index,&d);
  for(int i=0;i<d_rows;++i) d*=tmp(i,i);
  delete index;
  return d;
}


void Matrix::diagonaliseSymmetric(double *eigenValues){
  assert(d_rows==d_columns);
  
  double *e;
  e = new double[d_rows];
  tred2(d_val,d_rows,eigenValues,e);
  tqli(eigenValues,e,d_rows,d_val);
  sort(this,eigenValues);
  
  delete[] e;

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

// ludcmp of Numerical Recipes P. 46
// Modified to count from 0.

static void ludcmp(double **a, int n, int *indx, double *d){
    int i,imax=0,j,k; 
    double big,dum,sum,temp; 
    double *vv; 
    vv=new (double [n]);
    *d=1.0; 
    for (i=0;i<n;i++) { 
      big=0.0; 
      for (j=0;j<n;j++) 
	if ((temp=fabs(a[i][j])) > big) big=temp; 
      if (big == 0.0) throw Matrix::Exception("ludcmp: Matrix singular.");
      vv[i]=1.0/big; 
    } 
    for (j=0;j<n;j++) { 
      for (i=0;i<j;i++) { 
	sum=a[i][j]; 
	for (k=0;k<i;k++) sum -= a[i][k]*a[k][j]; 
	a[i][j]=sum; 
      } 
      big=0.0; 
      for (i=j;i<n;i++) { 
	sum=a[i][j]; 
	for (k=0;k<j;k++)
	  sum -= a[i][k]*a[k][j]; 
	a[i][j]=sum; 
	if ( (dum=vv[i]*fabs(sum)) >= big) { 
	  big=dum; imax=i; 
	} 
      } 
      if (j != imax) { 
	for (k=0;k<n;k++) { 
	  dum=a[imax][k]; 
	  a[imax][k]=a[j][k]; 
	  a[j][k]=dum; 
	} 
	*d = -(*d); 
	vv[imax]=vv[j]; 
      } 
      indx[j]=imax; 
      if (a[j][j] == 0.0) a[j][j]=TINY; 
      if (j != n) { 
	dum=1.0/(a[j][j]); 
	for (i=j+1;i<n;i++) a[i][j] *= dum; 
      } 
    } 
    delete vv;
}

static void tred2(double **a, int n, double d[], double e[]){
  int l,k,j,i; 
  double scale,hh,h,g,f; 
  for (i=n-1;i>0;i--) { 
    l=i-1; 
    h=scale=0.0; 
    if(l > 0) { 
      for (k=0;k<=l;k++) scale += fabs(a[i][k]); 
      if (scale == 0.0) 
	e[i]=a[i][l]; 
      else { 
	for (k=0;k<=l;k++) { 
	  a[i][k] /= scale; 
	  h += a[i][k]*a[i][k]; // Form Ù in h. 
	} 
	f=a[i][l]; 
	g=(f >= 0.0 ? -sqrt(h) : sqrt(h)); 
	e[i]=scale*g; 
	h -= f*g; //Now h is equation (11.2.4). 
	a[i][l]=f-g; //Store u in the ith row of a. 
	f=0.0; 
	for (j=0;j<=l;j++) { 
	  /* Next statement can be omitted if eigenvectors not wanted */ 
	  a[j][i]=a[i][j]/h; // Store u=H in ith column of a. 
	  g=0.0; // Form an element of A   u in g. 
	  for (k=0;k<=j;k++) g += a[j][k]*a[i][k]; 
	  for (k=j+1;k<=l;k++) g += a[k][j]*a[i][k]; 
	  e[j]=g/h; // Form element of p in temporarily unused element of e.
	  f += e[j]*a[i][j]; 
	} 
	hh=f/(h+h); // Form K, equation (11.2.11). 
	for (j=0;j<=l;j++) { //Form q and store in e overwriting p. 
	  f=a[i][j]; 
	  e[j]=g=e[j]-hh*f; 
	  for (k=0;k<=j;k++) // Reduce a, equation (11.2.13). 
	    a[j][k] -= (f*e[k]+g*a[i][k]); 
	} 
      } 
    } 
    else e[i]=a[i][l]; d[i]=h; 
  } 
  /* Next statement can be omitted if eigenvectors not wanted */ 
  d[0]=0.0; e[0]=0.0; 
  /* Contents of this loop can be omitted if eigenvectors not wanted except for statement d[i]=a[i][i]; */ 
  for (i=0;i<n;i++) { // Begin accumulation of transformation ma- trices. 
    l=i-1; 
    if (d[i]) { // This block skipped when i=1. 
      for (j=0;j<=l;j++) { 
	g=0.0; 
	for (k=0;k<=l;k++) // Use u and u=H stored in a to form P   Q. 
	  g += a[i][k]*a[k][j]; 
	for (k=0;k<=l;k++) a[k][j] -= g*a[k][i]; 
      } 
    } 
    d[i]=a[i][i]; // This statement remains. 
    a[i][i]=1.0; // Reset row and column of a to identity matrix for next iteration. 
    for (j=0;j<=l;j++) a[j][i]=a[i][j]=0.0; 
  } 
}

inline double pythag(double a, double b)
{
  double absa,absb;
  absa=fabs(a);
  absb=fabs(b);
  if (absa > absb) return absa*sqrt(1.0+(absb/absa)*(absb/absa));
  else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+(absa/absb)*(absa/absb)));
}

static void tqli(double d[], double e[], int n, double **z){
  int m,l,iter,i,k; 
  double s,r,p,g,f,dd,c,b; 
  for (i=1;i<n;i++) e[i-1]=e[i]; //Convenient to renumber the el- ements of e. 
  e[n-1]=0.0; 
  for (l=0;l<n;l++) { 
    iter=0; 
    do { 
      for (m=l;m<n-1;m++) { // Look for a single small subdi- agonal element to split the matrix. 
	dd=fabs(d[m])+fabs(d[m+1]); 
	if ((double)(fabs(e[m])+dd) == dd) break; 
      } 
      if (m != l) { 
	if (iter++ == 30) throw Matrix::Exception("Too many iterations in tqli"); 
	g=(d[l+1]-d[l])/(2.0*e[l]); // Form shift. 
	r=pythag(g,1.0); 
	//	g=d[m]-d[l]+e[l]/(g+SIGN(r,g)); // This is dm   ks. 
        g=d[m]-d[l]+e[l]/(g+((g>=0) ? fabs(r) : -fabs(r)));
	s=c=1.0; 
	p=0.0;
	for (i=m-1;i>=l;i--) { // A plane rotation as in the origi- nal QL, 
	  // followed by Givens rotations to restore tridiag- onal form. 
	  f=s*e[i]; 
	  b=c*e[i]; 
	  e[i+1]=(r=pythag(f,g));
	  if (r == 0.0) { // Recover from under ow. 
	    d[i+1] -= p; 
	    e[m]=0.0; 
	    break; 
	  } 
	  s=f/r; 
	  c=g/r; 
	  g=d[i+1]-p; 
	  r=(d[i]-g)*s+2.0*c*b; 
	  d[i+1]=g+(p=s*r); 
	  g=c*r-b; 
	  /* Next loop can be omitted if eigenvectors not wanted*/ 
	  for (k=0;k<n;k++) { // Form eigenvectors. 
	    f=z[k][i+1]; 
	    z[k][i+1]=s*z[k][i]+c*f; 
	    z[k][i]=c*z[k][i]-s*f; 
	  } 
	} 
	if (r == 0.0 && i >= l) continue; 
	d[l] -= p; 
	e[l]=g; 
	e[m]=0.0; 
      } 
    } 
    while (m != l); 
  } 
}

static void sort(Matrix *mat, double* eigenvalues) {

  assert(mat->rows()==mat->columns());        // only quadratic!

  int min;
  double temp;

  for(int i=0;i<(mat->rows()-1);i++) {
    min = i;
    for(int j=i+1;j<mat->rows();j++) {
      if(eigenvalues[j] > eigenvalues[min])
        min = j;
    }
    // swap min,i;
    temp = eigenvalues[min];
    eigenvalues[min]=eigenvalues[i];
    eigenvalues[i] = temp;
    swapColumns(mat,i,min);         // sort eigenvectors
  }
}

static void swapColumns(Matrix *mat, int i,int j) {
  double temp;
  
  for(int k=0;k<mat->rows();k++) {
    temp = (*mat)(k,i);
    (*mat)(k,i) = (*mat)(k,j);
    (*mat)(k,j) = temp;
  }
}

}
