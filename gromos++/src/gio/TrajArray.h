// TrajArray.h
#ifndef TRAJARRAY_H
#define TRAJARRAY_H
#include "../gcore/Box.h"
#include "../gcore/Molecule.h"
#include "../gcore/System.h"

class TrajArray {

  public:
    // Constructors
    // nfram is number of frames to hold in array
    TrajArray(const gcore::System &sys, const unsigned int nfram);
    TrajArray(const gcore::Molecule &mol, const unsigned int nfram);

    // Destructor
    ~TrajArray();

    // Function to store a frame
    bool store( const gcore::System &sys,
                const unsigned int frameidx );

    bool store( const gcore::Molecule &mol,
                const gcore::Box &box,
                const unsigned int frameidx );

    // Function to extract a frame
    bool extract( gcore::System &sys,
                  const unsigned int frameidx ) const;

    bool extract( gcore::Molecule &mol,
                  gcore::Box &box,
                  const unsigned int frameidx ) const;

    // Accessor for framesize
    unsigned int num_of_atoms();

  private:
    // the array of data
    double* trajdata;
    // the framesize
    unsigned int framesize;
    // the number of frames in the array
    unsigned int numframes;
};
#endif                                                                 
