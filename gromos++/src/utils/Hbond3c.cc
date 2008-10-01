#include <cassert>
#include <iostream>
#include <cstdio>
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
#include "Hbond3c.h"
#include "Neighbours.h"

using namespace std;

using utils::Hbond3c;


void Hbond3c::adddistances(double i, double j)
{
  d_dist_don_A1 += i;
  d_dist_don_A2 += j;
} //end adddistance

void Hbond3c::addangles(double i, double j)
{
  d_angle_b_don_A1 += i;
  d_angle_b_don_A2 += j;
} //end addangle

void Hbond3c::addanglesum(double i)
{
  d_angle_sum += i;
}

void Hbond3c::adddihedral(double i)
{
  d_dihedral += i;
}

void Hbond3c::addnum()
{
  ++d_num;
}

void Hbond3c::setIndices(int donor, int boundto, int acceptor1, int acceptor2)
{
  d_don_ind = donor;
  d_b_ind = boundto;
  d_ac1_ind = acceptor1;
  d_ac2_ind = acceptor2;
} 


void Hbond3c::calcmean()
{
  d_mean_dist_don_A1    = d_dist_don_A1 / d_num;
  d_mean_dist_don_A2    = d_dist_don_A2 / d_num;
  d_mean_angle_b_don_A1 = d_angle_b_don_A1 / d_num;
  d_mean_angle_b_don_A2 = d_angle_b_don_A2 / d_num;
  d_mean_angle_sum      = d_angle_sum / d_num;
  d_mean_dihedral       = d_dihedral / d_num;
  
} //end Hbondcalc::calcmean()

void Hbond3c::clear()
{
  d_dist_don_A1 = 0;
  d_dist_don_A2 = 0;
  d_angle_b_don_A1 = 0;
  d_angle_b_don_A2 = 0;
  d_angle_sum = 0;
  d_dihedral = 0;
  d_num = 0;
}


