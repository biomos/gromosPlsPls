// utils_Cluster.cc

#ifndef INCLUDED_GROMOS_EXCEPTION
#include "../gromos/Exception.h"
#endif

#include "RmsdMat.h"

/*
 * RmsdMat
*/
using utils::RmsdMat;

RmsdMat::RmsdMat(const unsigned int n) : width(n), matrix(n * (n - 1) / 2) {}

RmsdMat::~RmsdMat(){

  matrix.clear();
}

unsigned int RmsdMat::index(const unsigned int i, 
  const unsigned int j) const{

  int k, l;
  unsigned int r;
 
  if(i == j)
    throw gromos::Exception("RmsdMat",
      "Can't handle element. RmsdMat does not have a diagonal.");
  else if(i >= width || j >= width)
    throw gromos::Exception("RmsdMat",
      "Can't handle element. Index is larger than matrix width.");

  if(i < j){
    k = i;
    l = j;
  }
  else{
    k = j;
    l = i;
  }

  r = k * (width - 1) - (k * k + k) / 2 + l - 1;

  return r;
}

void RmsdMat::insert(const unsigned int i, 
  const unsigned int j, const float rmsd){

  matrix[index(i, j)] = rmsd;
}

float RmsdMat::retrieve(const unsigned int i, const unsigned int j) const{

  return matrix[index(i, j)];
}

unsigned int RmsdMat::size_n() const{

  return matrix.size();
}
