#include <cassert>
#include <iostream>
#include <stdio.h>
#include <string>
#include <numeric>
#include <functional>
#include "../args/Arguments.h"
#include "../gio/Ginstream.h"
#include "../gcore/System.h"
#include "../gcore/Molecule.h"
#include "../gcore/MoleculeTopology.h"
#include "../gcore/Solvent.h"
#include "../gcore/SolventTopology.h"
#include "../gcore/AtomTopology.h"
#include "../gmath/Vec.h"

#include "AtomSpecifier.h"
#include "Hbond.h"
#include "Neighbours.h"

using namespace std;

using utils::Hbond;



void Hbond::adddistance(double i)
{

  d_dist += i;
  
} //end adddistance

void Hbond::addangle(double i)
{
  
  d_angle += i;
  
} //end addangle

void Hbond::addnum()
{
  ++d_num;
}

void Hbond::setIndices(int donor, int boundto, int acceptor)
{
  d_don_ind = donor;
  d_b_ind = boundto;
  d_ac_ind = acceptor;
  
} 




void Hbond::calcmean()
{
  d_mean_dist  = d_dist/d_num;
  d_mean_angle = d_angle/d_num;
  
} //end Hbondcalc::calcmean()

void Hbond::clear()
{
  d_dist = 0;
  d_angle = 0;
  d_num = 0;
}


