// TrajArray.cc

#include "../gmath/Vec.h"
#include "../gromos/Exception.h"
#include "TrajArray.h"

#include <algorithm>
#include <iostream>

using gcore::Box;
using gcore::Molecule;
using gcore::System;
using gmath::Vec;

// Constructor
namespace utils{

TrajArray::TrajArray(const System &sys) {

  nAtoms = 0;

  // count the atoms in the system to find frame size
  for (int molIndex = 0; molIndex < sys.numMolecules(); molIndex++) {
    nAtoms += sys.mol(molIndex).numAtoms();
  } 
}

// Destructor
TrajArray::~TrajArray() {
  for(unsigned int i = 0; i < trajectoryData.size(); i++)
    delete [] trajectoryData[i];
}

// Function to store a frame
void TrajArray::store(const gcore::System &sys, 
  const unsigned int frameIndex){

  int i;
  int nAtomsMol, molAtomIndex, molIndex;
  unsigned int nAtomsSystem = 0;
  double *framePointer;

  // count the atoms in the system to find frame size
  for (molIndex = 0; molIndex < sys.numMolecules(); molIndex++)
    nAtomsSystem += sys.mol(molIndex).numAtoms();

  if (frameIndex < 0 || nAtomsSystem != nAtoms )
    throw gromos::Exception("TrajArray", "Unable to store frame.\n");

  if(frameIndex + 1 >= trajectoryData.size()){
  // guard againts off-by-one
    trajectoryData.resize(frameIndex + 1);
  }
  else
    delete [] trajectoryData[frameIndex];

  trajectoryData[frameIndex] = new double[3 * nAtoms];
  framePointer = trajectoryData[frameIndex];

  // read all the coords from the sys and store them in the array
  for (molIndex = 0; molIndex < sys.numMolecules(); molIndex++) {
    nAtomsMol = sys.mol(molIndex).numAtoms();
    for (molAtomIndex = 0; molAtomIndex < nAtomsMol; molAtomIndex++) {
      for(i = 0; i < 3; i++){
        *framePointer = sys.mol(molIndex).pos(molAtomIndex)[i];
        framePointer++;
      }
    }
  }
}

// Function to extract a frame
void TrajArray::extract( gcore::System &sys,
  const unsigned int frameIndex ) const {

  int i, molIndex, nAtomsMol, molAtomIndex;
  unsigned int nAtomsSystem = 0;
  double *framePointer;

  // count the atoms in the system to find frame size
  for (molIndex = 0; molIndex < sys.numMolecules(); molIndex++)
    nAtomsSystem += sys.mol(molIndex).numAtoms();
   
  if (frameIndex < 0 || frameIndex >= trajectoryData.size())
    throw gromos::Exception("TrajArray", 
      "Can't extract. Invalid frame index.\n");
  if(nAtomsSystem != nAtoms)
    throw gromos::Exception("TrajArray", 
      "Can't extract. Invalid system.\n");
  if(!(framePointer = trajectoryData[frameIndex]))
    throw gromos::Exception("TrajArray", 
      "Can't extract. Empty frame.\n");

  for (molIndex = 0; molIndex < sys.numMolecules(); molIndex++) {
    nAtomsMol = sys.mol(molIndex).numAtoms();
    for (molAtomIndex = 0; molAtomIndex < nAtomsMol; molAtomIndex++) {
      for(i = 0; i < 3; i++){
        sys.mol(molIndex).pos(molAtomIndex)[i] = *framePointer;
        framePointer++;
      }
    }
  }
}

unsigned int TrajArray::numAtoms() {
    return nAtoms;
}
}
