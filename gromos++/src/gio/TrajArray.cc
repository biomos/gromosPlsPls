// TrajArray.cc

#include "../gmath/Vec.h"
#include "TrajArray.h"

#include <algorithm>

using gcore::Box;
using gcore::Molecule;
using gcore::System;
using gmath::Vec;

// Constructors
TrajArray::TrajArray(const System &sys, const unsigned int nfram) {

  // a frame contains framsize*3 coordinates and 3 box dimensions
  framesize=0;
  numframes=0;

  // count the atoms in the system to find frame size
  for (int mm=0; mm<sys.numMolecules(); mm++) {
    framesize+=sys.mol(mm).numAtoms();
  } 
  unsigned int totalsize=nfram*(framesize+1)*3;
  trajdata=new double[totalsize];
  for (unsigned int ii=0; ii<totalsize; ii++) {
    trajdata[ii]=0.0;
  }
  numframes=nfram;
}

TrajArray::TrajArray(const Molecule &mol, const unsigned int nfram) {

  // a frame contains framsize*3 coordinates and 3 box dimensions
  framesize=0;
  numframes=0;

  // count the atoms in the molecule to find frame size
  framesize=mol.numAtoms();

  unsigned int totalsize=nfram*(framesize+1)*3;
  trajdata=new double[totalsize];
  for (unsigned int ii=0; ii<totalsize; ii++) {
    trajdata[ii]=0.0;
  }
  numframes=nfram;
}


// Destructor
TrajArray::~TrajArray() {
  if (numframes > 0) {
    delete [] trajdata;
  }
}


// Function to store a frame
bool TrajArray::store( const gcore::System &sys,
                       const unsigned int frameidx ) {
  // count the atoms in the system to find frame size
  unsigned int sys_size=0;
  for (int mm=0; mm<sys.numMolecules(); mm++) {
    sys_size+=static_cast<unsigned int>(sys.mol(mm).numAtoms());
  } 
  if (  (frameidx < 0) ||
        (frameidx >= numframes) ||
        (sys_size != framesize)  ) {
    return false;
  }
  else {
    unsigned int iimin=frameidx*((framesize*3)+3);
    unsigned int iimax=iimin+(framesize*3);
    unsigned int jj=0;
    unsigned int natoms=0;
    for (int mm=0; mm<sys.numMolecules(); mm++) {
      natoms=static_cast<unsigned int>(sys.mol(mm).numAtoms());
      for (unsigned int ii=0; ii < natoms ; ii++) {
        jj=iimin+3*ii;
        trajdata[jj]  =sys.mol(mm).pos(ii)[0];
        trajdata[jj+1]=sys.mol(mm).pos(ii)[1];
        trajdata[jj+2]=sys.mol(mm).pos(ii)[2];
      }
      iimin+=natoms*3;
    }
    trajdata[iimax]  =sys.box()[0];
    trajdata[iimax+1]=sys.box()[1];
    trajdata[iimax+2]=sys.box()[2];
    return true;
  }
}

bool TrajArray::store( const gcore::Molecule &mol,
                       const gcore::Box &box,
                       const unsigned int frameidx ) {
  if (  (frameidx < 0) ||
        (frameidx >= numframes) ||
        (static_cast<unsigned int>(mol.numAtoms()) != framesize)  ) {
    return false;
  }
  else {
    unsigned int iimin=frameidx*((framesize*3)+3);
    unsigned int iimax=iimin+(framesize*3);
    unsigned int jj=0;
    for (unsigned int ii=0; ii < framesize ; ii++) {
      jj=iimin+3*ii;
      trajdata[jj]  =mol.pos(ii)[0];
      trajdata[jj+1]=mol.pos(ii)[1];
      trajdata[jj+2]=mol.pos(ii)[2];
    }
    trajdata[iimax]  =box[0];
    trajdata[iimax+1]=box[1];
    trajdata[iimax+2]=box[2];
    return true;
  }
}


// Function to extract a frame
bool TrajArray::extract( gcore::System &sys,
                         const unsigned int frameidx ) const {

  // count the atoms in the system to find frame size
  unsigned int sys_size=0;
  for (int mm=0; mm<sys.numMolecules(); mm++) {
    sys_size+=static_cast<unsigned int>(sys.mol(mm).numAtoms());
  } 
  if (  (frameidx < 0) ||
        (frameidx >= numframes) ||
        (sys_size != framesize)  ) {
    return false;
  }
  else {
    unsigned int iimin=frameidx*((framesize*3)+3);
    unsigned int iimax=iimin+(framesize*3);
    unsigned int jj=0;
    unsigned int natoms=0;
    for (int mm=0; mm<sys.numMolecules(); mm++) {
      natoms=static_cast<unsigned int>(sys.mol(mm).numAtoms());
      for (unsigned int ii=0; ii < natoms ; ii++) {
        jj=iimin+3*ii;
        sys.mol(mm).pos(ii)[0]=trajdata[jj];
        sys.mol(mm).pos(ii)[1]=trajdata[jj+1];
        sys.mol(mm).pos(ii)[2]=trajdata[jj+2];
      }
      iimin+=natoms*3;
    }
    sys.box()[0]=trajdata[iimax];
    sys.box()[1]=trajdata[iimax+1];
    sys.box()[2]=trajdata[iimax+2];
    return true;
  }

}

bool TrajArray::extract( gcore::Molecule &mol,
                         gcore::Box &box,
                         const unsigned int frameidx ) const {
  if (  (frameidx < 0) ||
        (frameidx >= numframes) ||
        (static_cast<unsigned int>(mol.numAtoms()) != framesize)  ) {
    return false;
  }
  else {
    unsigned int iimin=frameidx*((framesize*3)+3);
    unsigned int iimax=iimin+(framesize*3);
    unsigned int jj=0;
    for (unsigned int ii=0; ii < framesize ; ii++) {
      jj=iimin+3*ii;
      mol.pos(ii)[0]=trajdata[jj];
      mol.pos(ii)[1]=trajdata[jj+1];
      mol.pos(ii)[2]=trajdata[jj+2];
    }
    box[0]=trajdata[iimax];
    box[1]=trajdata[iimax+1];
    box[2]=trajdata[iimax+2];
    return true;
  }
}


inline unsigned int TrajArray::num_of_atoms() {
    return framesize;
}

