// fit_RotationalFit.cc
//includes explicit calls to gsl now

#include <cassert>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <iostream>
#include <vector>

#include "../gmath/Matrix.h"
#include "../gmath/Vec.h"

#include "FastRotationalFit.h"

using gmath::Matrix;
using gmath::Vec;
using namespace fit;
using namespace std;



double FastRotationalFit::rmsd(Matrix const &rot, 
			       vector<Vec> const &ref, 
			       vector<Vec> const &sys)const
{
  double rmsd2=0;
  double temp;
  if(d_rmsd_spec.size()){
    for(size_t i=0; i< ref.size(); ++i){
      for(size_t j=0; j< 3; ++j){
	if(d_rmsd_spec[i]){
	  temp=0;
	  for(int b=0;b<3;++b)
	    temp += rot(j,b) * sys[i][b];
	
	  rmsd2 += (ref[i][j] - temp) * (ref[i][j] - temp);
	}
      }
    }
    
    rmsd2 /= d_rmsd_num_atoms;
  }
  else{
    for(size_t i=0; i< ref.size(); ++i){
      for(size_t j=0; j< 3; ++j){
	temp=0;
	for(int b=0;b<3;++b)
	  temp += rot(j,b) * sys[i][b];
      
	rmsd2 += (ref[i][j] - temp) * (ref[i][j] - temp);
      }
    }
    
    rmsd2 /= ref.size();
  }
  
  return sqrt(rmsd2);
}

int FastRotationalFit::fit(vector<Vec> const &ref, 
			   vector<Vec> &sys)const{

  Matrix r(3,3,0);
  int error=fit(r,ref, sys);
  if(error)
    return error;
  size_t num = ref.size();
  
  for(size_t n=0; n < num; ++n){
    sys[n] = r * sys[n];
  }
  return 0;
}

  
int FastRotationalFit::fit(Matrix &rot,
			   vector<Vec> const &ref, 
			   vector<Vec> const &sys)const{
  
  size_t num = ref.size();
  
  Matrix U(3,3,0);
  if(d_fit_spec.size()){
    
    for(int i=0;i<3;++i)
      for(int j=0;j<3;++j)
	for(size_t n=0;n<num;++n)
	  if(d_fit_spec[n])
	    U(i,j)+= sys[n][i]* ref[n][j];
    U *= 1.0/d_fit_num_atoms;
  }
  else{
    for(int i=0;i<3;++i)
      for(int j=0;j<3;++j)
	for(size_t n=0;n<num;++n)
	  U(i,j)+= sys[n][i]* ref[n][j];
    U *= 1.0/num;
  }
  
  double det=U.fastdet3X3Matrix();
  
  int signU = ( det>0 ? 1 : -1);


  gsl_matrix * omega = gsl_matrix_alloc (6, 6);
  gsl_matrix_set_zero (omega);

  for(int i=0;i<3;++i){
    for(int j=0;j<3;++j){
     gsl_matrix_set (omega, i, j+3,U(i,j));
     gsl_matrix_set (omega, i+3, j,U(j,i));
    }
  }

  
  double *eigenvals = new double [6];
  
  gsl_vector *eval = gsl_vector_alloc (6);
  gsl_matrix *evec = gsl_matrix_alloc (6,6);

  gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (6);

  gsl_eigen_symmv (omega, eval, evec, w);

  gsl_eigen_symmv_free(w);

  gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_VAL_DESC);

   Matrix Omega(6,6,0);
   for (int i=0; i < 6; ++i){
    eigenvals[i] = gsl_vector_get(eval, i);
   for (int j=0; j < 6; ++j){
     Omega(i,j)=gsl_matrix_get(evec, i, j);
   }
   }

  gsl_matrix_free (omega);
  gsl_matrix_free (evec);
  gsl_vector_free (eval);

  //  if(det<0){
  //  delete[] eigenvals;
  //  return -1;
  //}
  if(det<0 && fabs(eigenvals[1]-eigenvals[2])<1.0e-5){
    delete[] eigenvals;
    return -2;
  }

  // Extract vectors from Omega.  
  Omega *= sqrt(2.0);
  Vec k1(Omega(0,0), Omega(1,0), Omega(2,0));
  Vec k2(Omega(0,1), Omega(1,1), Omega(2,1));
  Vec k3(Omega(0,2), Omega(1,2), Omega(2,2));
  Vec h1(Omega(3,0), Omega(4,0), Omega(5,0));
  Vec h2(Omega(3,1), Omega(4,1), Omega(5,1));
  Vec h3(Omega(3,2), Omega(4,2), Omega(5,2));

  double spat = h1.dot(h2.cross(h3));
  
  // turn 3rd vectors
  if(spat<0){
    h3=-h3;
    k3=-k3;
  }

  rot = Matrix(h1,k1) + Matrix(h2,k2) + signU*Matrix(h3,k3);

  delete[] eigenvals;

  return 0;
}
