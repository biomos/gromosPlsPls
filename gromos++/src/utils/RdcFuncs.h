/*
 * This file is part of GROMOS.
 * 
 * Copyright (c) 2011, 2012, 2016, 2018, 2021, 2023 Biomos b.v.
 * See <https://www.gromos.net> for details.
 * 
 * GROMOS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

/* 
 * File:   RdcFuncs.h
 * Author: jra, lnw
 *
 * Created on November 24, 2009, 4:09 PM, improved mid-2015
 */

#ifndef RDCFUNCS_H
#define RDCFUNCS_H

#include <iomanip>
#include <vector>
#include <string>
#include <ostream>

#include "../gromos/Exception.h"
#include "../utils/VirtualAtom.h"

namespace utils {

// struct for storing data to describe one specific RDC
class rdcparam {
  rdcparam();

  public:
  VirtualAtom atom1;
  VirtualAtom atom2;
  // parameters
  double w;    // weight factor
  double exp;  // experimental rdc
  double dD0;  // allowed deviation (half of flat-bottom width)
  double gi;   // gyromagnetic ratio of atom i
  double gj;   // gyromagnetic ratio of atom j
  double r0;  // standard distance read from RDCRESSPEC
  double dmax; // maximum possible rdc for atoms ij (and ik)
  double dmax_r3; // maximum possible rdc for atoms ij (and ik) * rij^3
  double dij; // the actual distance between atoms i and j
  gmath::Vec rij;
  std::string atomnum1, atomnum2, atomname1, atomname2;


  rdcparam(VirtualAtom atom1, VirtualAtom atom2, double r0, double gi, double gj, double exp, double dD0, double w);

  rdcparam(const rdcparam &rdcp) : atom1(rdcp.atom1), atom2(rdcp.atom2),
              r0(rdcp.r0), gi(rdcp.gi), gj(rdcp.gj), w(rdcp.w), exp(rdcp.exp), 
              dD0(rdcp.dD0), dmax(rdcp.dmax), dmax_r3(rdcp.dmax_r3),
              atomnum1(rdcp.atomnum1), atomnum2(rdcp.atomnum2),
              atomname1(rdcp.atomname1), atomname2(rdcp.atomname2), dij(rdcp.dij), rij(rdcp.rij) {}

  rdcparam &operator=(const rdcparam &rdcp) {
    atom1 = rdcp.atom1;
    atom2 = rdcp.atom2;
    w = rdcp.w;
    exp = rdcp.exp;
    gi = rdcp.gi;
    gj = rdcp.gj;
    dD0 = rdcp.dD0;
    r0 = rdcp.r0;
    dmax = rdcp.dmax;
    dmax_r3 = rdcp.dmax_r3;
    atomnum1 = rdcp.atomnum1;
    atomnum2 = rdcp.atomnum2;
    atomname1 = rdcp.atomname1;
    atomname2 = rdcp.atomname2;
    dij = rdcp.dij;
    rij = rdcp.rij;
    return *this;
  }

  friend std::ostream& operator<<(std::ostream& s, const rdcparam& rdc){
    s << "{" << rdc.atom1.toString() << ", " << rdc.atom2.toString() << ", " 
      << std::setprecision(5) << rdc.w << ", " << std::scientific << rdc.exp;
    s.unsetf(std::ios_base::floatfield);
    s << ", " << rdc.gi << ", " << rdc.gj <<  ", " << rdc.r0 <<  ", " << rdc.dD0 <<  ", "
      << std::scientific << rdc.dmax <<  ", " << rdc.dij << "}";
    return s;
  }
};

typedef std::vector<rdcparam> rdcdata_t;


// function to read in RDC data
rdcdata_t read_rdc(const std::vector<std::string> &buffer, gcore::System &sys);

// function to read in RDC grousp
std::vector<std::vector<unsigned int> > read_groups(const std::vector<std::string> &buffer, const unsigned n_rdc);

// compute the RDC vector
void calc_rij_dmax(rdcdata_t &fit_data, bool norm_r0);

// compute the coefficients of the matrix describing bond vector fluctuations for fitting
void calc_coef_fit(const gcore::System &sys, const rdcdata_t &fit_data, double coef_mat[]);

// compute the coefficients of the matrix describing bond vector fluctuations for back-calculation
void calc_coef_bc(const gcore::System &sys, const rdcdata_t &fit_data, double coef_mat_j[], double coef_mat_k[]);

// get inter spin vector, normalised (\vec r / |\vec r|)
gmath::Vec get_inter_spin_vector_normalised(const gcore::System &sys, const rdcdata_t &fit_data, int index);

//  returns {xx, yy, xy, xz, yz}
std::vector<double> lls_fit(const gcore::System &sys, const rdcdata_t &fit_data, double coef_mat[]);
std::vector<double> svd_fit(const gcore::System &sys, const rdcdata_t &fit_data, double coef_mat[]);

void diagonalize_tensor(const std::vector<double> &t5, std::vector<double> &eval, std::vector<std::vector<double> >  &evec);


// calculate Q according to Cornilescu
double calc_Q(const std::vector<double> &calc, const std::vector<double> &expt);
// calculate the R value
double calc_R(const std::vector<double> &calc, const std::vector<double> &expt);
// calculate the RMSD
double calc_RMSD(const std::vector<double> &calc, const std::vector<double> &expt, const rdcdata_t &rdc_data);

}

#endif /* _RDCFUNCS_H */

