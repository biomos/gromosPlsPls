// fit_RotationalFit.cc

#include "RotationalFit.h"
#include "Reference.h"
#include "PositionUtils.h"
#include "../gcore/System.h"
#include "../gcore/Molecule.h"
#include "../gmath/Matrix.h"
#include "../gmath/Vec.h"

using fit::RotationalFit;
using fit::RotationalFit_i;
using fit::Reference;
using gcore::System;
using gcore::Molecule;
using gmath::Matrix;
using gmath::Vec;

// constructs the rotation Matrix
static void rotationMatrix(gmath::Matrix *mat, const gcore::System &mol, const fit::Reference &w);

RotationalFit::RotationalFit(Reference *w){
  d_ref=w;
  PositionUtils::shiftToCog(&w->sys(),*w);
}

RotationalFit::~RotationalFit(){}

void RotationalFit::fit(gcore::System *sys)const{
  PositionUtils::shiftToCog(sys,*d_ref);
  Matrix rot(3,3);
  rotationMatrix(&rot,*sys,*d_ref);
  PositionUtils::rotate(sys,rot);
}

static void rotationMatrix(Matrix *mat, const System &sys, const Reference &r){

  const System &ref = r.sys();
  
  Matrix U(3,3,0);
  for(int i=0;i<3;++i)
    for(int j=0;j<3;++j)
      for(int m=0;m<ref.numMolecules();++m)
	for(int n=0;n<ref.mol(m).numAtoms();++n)
	  if(r.weight(m,n))
	    U(i,j)+=r.weight(m,n)*sys.mol(m).pos(n)[i]*ref.mol(m).pos(n)[j];
  double det=U.det();
  int signU = ( det>0 ? 1 : -1);
  
  Matrix Omega(6,6,0);
  for(int i=0;i<3;++i)
    for(int j=0;j<3;++j){
      Omega(i,j+3) = U(i,j);
      Omega(i+3,j) = U(j,i);
    }
  
  double *eigenvals = new double [6];
  Omega.diagonaliseSymmetric(eigenvals);
  if(det<0&&fabs(eigenvals[1]-eigenvals[2])<1.0e-5)
    throw RotationalFit::Exception("Rotation matrix degenerate!");

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

  *mat = Matrix(h1,k1) + Matrix(h2,k2) + signU*Matrix(h3,k3);

  delete eigenvals;

}
