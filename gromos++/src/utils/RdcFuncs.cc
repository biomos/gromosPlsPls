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
 * File:   RdcFuncs.cc
 * Author: jra, lnw
 *
 * Created 2009, improved mid-2015
 */
#include "RdcFuncs.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <math.h>
#include <sstream>
#include <string>
#include <vector>

// //time measurement
// #include <sys/time.h>

#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_vector_double.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

#include "../gcore/AtomTopology.h"     // added for DEBUG
#include "../gcore/Molecule.h"
#include "../gcore/MoleculeTopology.h" // added for DEBUG
#include "../gcore/System.h"
#include "../gmath/Physics.h"
#include "../gmath/Vec.h"
#include "../gromos/Exception.h"
#include "../utils/VirtualAtom.h"
#include "../utils/AtomSpecifier.h"
#include "../utils/debug.h"

using namespace std;
using namespace gcore;

//using gcore::MoleculeTopology; // added for DEBUG
//using gcore::AtomTopology;     // added for DEBUG


#define container_output(container) \
  template <typename T> ostream& operator<<(ostream& s, const container<T>& v) \
  { \
  s << "{"; \
  for(typename container<T>::const_iterator x(v.begin());x!=v.end();){ \
    s << *x; \
    if(++x!=v.end()) s << ","; \
  } \
  s << "}"; \
  return s; \
}
container_output(vector);



utils::rdcparam::rdcparam(VirtualAtom atom1, VirtualAtom atom2, double rij, double gi, double gj, double exp, double dD0, double w):
          atom1(atom1), atom2(atom2),
          rij(rij), gi(gi), gj(gj), w(w), exp(exp), 
          dD0(dD0) {

  // compute and store Dmax
  const double enumerator = -1.0 * gmath::physConst.get_mu0() * gmath::physConst.get_h() * gi * gj;
  const double denominator = pow(2.0 * M_PI * rij, 3);
  dmax = enumerator / denominator;
  if (atom1.type() == 0) {
    atomname1=atom1.conf().name(0);
    atomnum1=atom1.conf().toString(0);
  } else {
    atomname1="VA";
    atomnum1=atom1.toString();
  }
  if (atom2.type() == 0) {
    atomname2=atom2.conf().name(0);
    atomnum2=atom2.conf().toString(0);
  } else {
    atomname2="VA";
    atomnum2=atom2.toString();
  }

}


// // not currently used
// void utils::rdcparam::recalculate_bond_lengths(const System &sys){
//   // calculate the inter-nuclear distances from the reference coordinates
//   rij = (sys.mol(mol).pos(i) - sys.mol(mol).pos(j)).abs();
//   // and if it's an NH sidechain, calculate rik too
//   if (type == 5 || type == 6) {
//     rik = (sys.mol(mol).pos(i) - sys.mol(mol).pos(k)).abs();
//   }
// }

// function to read RDC data
utils::rdcdata_t utils::read_rdc(const vector<string> &buffer, System &sys, bool fit){
  rdcdata_t rdcp;
  double dish, disc;
  istringstream is(buffer[1]);
  is >> dish >> disc;
  for (unsigned int jj = 2; jj < buffer.size() - 1; jj++) {

    // first get atom numbers
    is.clear();
    is.str(buffer[jj]);
    int i, j, k, l, t;
    std::string atomnum1, atomnum2, atomname1, atomname2;
    is >> i >> j >> k >> l >> t;
    DEBUG(15, "i,j,k,l:" << i << ", " << j << ", " << k << ", " << l)
    if (i < 1) {
      throw gromos::Exception("read_rdc", "Disallowed atom number in RDC file\n" + buffer[jj]);
    }

    // adjust to gromos numbering
    VirtualAtom va1(sys, VirtualAtom::virtual_type(t), { i-1, j-1, k-1, l-1 }, dish, disc);
    is >> i >> j >> k >> l >> t;
    DEBUG(15, "i,j,k,l:" << i << ", " << j << ", " << k << ", " << l)
    if (i < 1) {
      throw gromos::Exception("read_rdc", "Disallowed atom number in RDC file\n" + buffer[jj]);
    }
    // adjust to gromos numbering
    VirtualAtom va2(sys, VirtualAtom::virtual_type(t), { i-1, j-1, k-1, l-1 }, dish, disc);

    double R0, G1, G2, D0, dD0, wrdc;

    is >> R0 >> G1 >> G2 >> D0 >> dD0 >> wrdc;

    // now get the other parameters
    DEBUG(15, "R0,G1,G2,D0,dD0,WRDC:" << R0 << ", " << G1 << ", " << G2 << ", " << D0 << ", " << dD0 << ", " << wrdc)

    if (is.fail()) throw gromos::Exception("read_rdc", "Bad line in RDC file\n" + buffer[jj]);

    const double atomic_mass_unit = 1.6605655e-27;
    const double elementary_charge = 1.6021892e-19;
    const double g_cf = (atomic_mass_unit / elementary_charge) * 1.0e6;
    // convert into gromos units
    D0 *= 1e-12; // s^-1 --> ps^-1
    dD0 *= 1e-12; // s^-1 --> ps^-1
    G1 *= g_cf;
    G2 *= g_cf;


    // store temporary RDC parameters
    rdcp.push_back(rdcparam(va1, va2, R0, G1, G2, D0, dD0, wrdc));
    DEBUG(15, "rdc before appending: " << rdcparam(va1, va2, R0, G1, G2, D0, dD0, wrdc))
  }

  return rdcp;
}


vector<vector<unsigned int> > utils::read_groups(const vector<string> &buffer, const unsigned n_rdc) {
  vector<vector<unsigned int> > rdc_groups;
  unsigned int int_buf = 0;
  
  for (unsigned int line = 1; line < buffer.size() - 1; line++) {
    istringstream is(buffer[line]);

    rdc_groups.push_back(vector<unsigned int>());
    while(is >> int_buf){ // rdc in groups are written starting at 1
      rdc_groups.back().push_back(int_buf);
    }
  }

  // test sanity
  for (unsigned i=0; i< rdc_groups.size(); i++){
    for (unsigned j=0; j< rdc_groups[i].size(); j++){
      if (rdc_groups[i][j] > n_rdc) {
        throw gromos::Exception("fit_rdc", "There are higher RDC group entries, than there are RDCs");
      }
    }
  }
  vector<unsigned int> group_entries(n_rdc, 0);
  for (unsigned i=0; i< rdc_groups.size(); i++){
    for (unsigned j=0; j< rdc_groups[i].size(); j++){
      group_entries[rdc_groups[i][j]-1]++;
    }
  }
  if (find(group_entries.begin(), group_entries.end(), 0) != group_entries.end()){
    cout << "There are RDCs which are in no group and which will be ignored.  That is possible, but may not be intentional." << endl;
  }
  DEBUG(7, "group entries: " << group_entries)

  return rdc_groups;
}


// compute the coefficients of the matrix describing bond vector fluctuations
void utils::calc_coef_fit(const System &sys, const rdcdata_t &fit_data, double coef_mat[]) {
  const int n_rdc = fit_data.size();
  const int n_ah = 5;
  for (int k=0; k<n_rdc; k++){
    Vec mu = get_inter_spin_vector_normalised(sys, fit_data, k);
    coef_mat[k*n_ah + 0] = mu[0] * mu[0] - mu[2] * mu[2];
    coef_mat[k*n_ah + 1] = mu[1] * mu[1] - mu[2] * mu[2];
    coef_mat[k*n_ah + 2] = 2.0 * mu[0] * mu[1];
    coef_mat[k*n_ah + 3] = 2.0 * mu[0] * mu[2];
    coef_mat[k*n_ah + 4] = 2.0 * mu[1] * mu[2];
    for(int h=0; h<5; ++h){
      DEBUG(12, "ck (" << k << "," << h << "): " << scientific << coef_mat[k*n_ah+h])
    }
  }
}


// compute the coefficients of the matrix describing bond vector fluctuations
// for back-calculated RDCs
/*
void utils::calc_coef_bc(const System &sys, const rdcdata_t &bc_data,
                         double coef_mat_j[], double coef_mat_k[]) {
  const int n_rdc = bc_data.size();
  const int n_ah = 5;
  for(int d=0; d<n_rdc; d++) {
    Vec mu_j = get_inter_spin_vector_normalised(sys, bc_data, d);
    // put into gsl matrix for ij
    coef_mat_j[d*n_ah + 0] = mu_j[0] * mu_j[0] - mu_j[2] * mu_j[2];
    coef_mat_j[d*n_ah + 1] = mu_j[1] * mu_j[1] - mu_j[2] * mu_j[2];
    coef_mat_j[d*n_ah + 2] = 2.0 * mu_j[0] * mu_j[1];
    coef_mat_j[d*n_ah + 3] = 2.0 * mu_j[0] * mu_j[2];
    coef_mat_j[d*n_ah + 4] = 2.0 * mu_j[1] * mu_j[2];

    // if side-chain NH, do it all again for ik
    if (bc_data[d].type == 5 || bc_data[d].type == 6) {
      unsigned int i = bc_data[d].i;
      unsigned int k = bc_data[d].k;
      unsigned int m = bc_data[d].mol;

      // get diff in coords for ik
      double mu_x_k = sys.mol(m).pos(i)[0] - sys.mol(m).pos(k)[0];
      double mu_y_k = sys.mol(m).pos(i)[1] - sys.mol(m).pos(k)[1];
      double mu_z_k = sys.mol(m).pos(i)[2] - sys.mol(m).pos(k)[2];
      const double mu_r_k = bc_data[d].rik;

      // scale
      mu_x_k /= mu_r_k;
      mu_y_k /= mu_r_k;
      mu_z_k /= mu_r_k;

      // put into gsl matrix for ik
      coef_mat_k[d*n_ah + 0] = mu_x_k * mu_x_k - mu_z_k * mu_z_k;
      coef_mat_k[d*n_ah + 1] = mu_y_k * mu_y_k - mu_z_k * mu_z_k;
      coef_mat_k[d*n_ah + 2] = 2.0 * mu_x_k * mu_y_k;
      coef_mat_k[d*n_ah + 3] = 2.0 * mu_x_k * mu_z_k;
      coef_mat_k[d*n_ah + 4] = 2.0 * mu_y_k * mu_z_k;

    } // end if side-chain NH
  }   // end RDC loop
}
*/


inline gmath::Vec utils::get_inter_spin_vector_normalised(const System &sys, const rdcdata_t &fit_data, int index) {
  // get atom numbers
  VirtualAtom va_i = fit_data[index].atom1;
  VirtualAtom va_j = fit_data[index].atom2;
  double mu_r = fit_data[index].rij;
  double mu_x, mu_y, mu_z;

  // RDCs for atoms that are explicitly present in GROMOS
 
  // get diff in coords
  mu_x = va_i.pos()[0] - va_j.pos()[0];
  mu_y = va_i.pos()[1] - va_j.pos()[1];
  mu_z = va_i.pos()[2] - va_j.pos()[2];

  DEBUG(12, "mu " << v2s(Vec(mu_x, mu_y, mu_z)))
  // normalise
  mu_x /= mu_r;
  mu_y /= mu_r;
  mu_z /= mu_r;
  DEBUG(12, "mu normalized" << v2s(Vec(mu_x, mu_y, mu_z)))

  return Vec(mu_x, mu_y, mu_z);
}


// calculate the Q value according to Cornilescu, Alexandrescu, Bax, Zweckstetter...
double utils::calc_Q(const vector<double> &calc, const vector<double> &expt) {
  double sdev2 = 0, sexpt2 = 0;
  for (unsigned int i = 0; i < calc.size(); i++) {
    sdev2 += pow(expt[i] - calc[i], 2);
    sexpt2 += expt[i] * expt[i];
  }
  return sqrt(sdev2 / sexpt2);
}


// calculate the R value (crystallographic)
double utils::calc_R(const vector<double> &calc, const vector<double> &expt) {
  double enumerator = 0, denominator = 0;
  for (unsigned int i = 0; i < calc.size(); i++) {
    enumerator += abs(abs(expt[i]) - abs(calc[i]));
    denominator += abs(expt[i]);
  }
  return enumerator / denominator;
}

// calculate the RMSD
double utils::calc_RMSD(const vector<double> &calc, const vector<double> &expt) {
  double sum = 0;
  for (unsigned int i = 0; i < calc.size(); i++) {
    sum += pow(expt[i] - calc[i], 2);
  }
  return sqrt(sum);
}

vector<double> utils::lls_fit(const System &sys, const rdcdata_t &fit_data, double coef_mat[]) {
//   struct timeval tv;
//   gettimeofday(&tv, NULL);
//   suseconds_t t0 = tv.tv_usec;

  const int n_rdc = fit_data.size();
  const int n_ah = 5;
  vector<double> tensor(n_ah, 0.0);

  gsl_vector *x = gsl_vector_alloc(n_ah);
  gsl_permutation *p = gsl_permutation_alloc(n_ah);
  int signum; // not used

  double matrix_array[n_ah*n_ah];
  for(int h=0; h<n_ah*n_ah; ++h){matrix_array[h]=0.0;} // init
  double result_array[n_ah];
  for(int h=0; h<n_ah; ++h){result_array[h]=0.0;} // init

//  gettimeofday(&tv, NULL);
//  suseconds_t t1 = tv.tv_usec;

  for(int k=0; k<n_rdc; k++){
    for(int h_prime=0; h_prime<n_ah; ++h_prime){
      for(int h=0; h<n_ah; ++h){
        matrix_array[h_prime*n_ah + h] += fit_data[k].w * coef_mat[k*n_ah + h] * coef_mat[k*n_ah + h_prime];
      }
      result_array[h_prime] += fit_data[k].w * (fit_data[k].exp / fit_data[k].dmax) * coef_mat[k*n_ah + h_prime];
    }
  }

  gsl_matrix_view m = gsl_matrix_view_array (matrix_array, n_ah, n_ah);
  gsl_vector_view b = gsl_vector_view_array (result_array, n_ah);

//  gettimeofday(&tv, NULL);
//  suseconds_t t2 = tv.tv_usec;

  // solve the SLE using gsl
  gsl_linalg_LU_decomp(&m.matrix, p, &signum);

// i might switch on the error handling later
//  // switch off and save error handler
//  gsl_error_handler_t *old_handler = gsl_set_error_handler_off();
//  int status = gsl_linalg_LU_solve(m, p, b, x);
  gsl_linalg_LU_solve(&m.matrix, p, &b.vector, x);
//  gsl_set_error_handler(old_handler);

//  gettimeofday(&tv, NULL);
//  suseconds_t t3 = tv.tv_usec;

//  if (status == GSL_SUCCESS) { // copy back a values
    for (int h = 0; h < n_ah; ++h) {
      tensor[h] = gsl_vector_get(x, h); // [nm]
//      cout << gsl_vector_get(x, h) << endl;
    }
//   } else if (status == GSL_EDOM) {
//     // ignore singular matrices which could very rarely occur
//     // do not update a_h
//   } else {
//     cout << "something unusual happend." << endl;
//   }

  gsl_permutation_free(p);
  gsl_vector_free(x);

//  gettimeofday(&tv, NULL);
//  suseconds_t t4 = tv.tv_usec;
//  cout << "time: " << t1-t0 << " us (allocate)" << endl;
//  cout << "time: " << t2-t1 << " us (prepare)" << endl;
//  cout << "time: " << t3-t2 << " us (do stuff)" << endl;
//  cout << "time: " << t4-t3 << " us (write out, deallocate)" << endl;
//  cout << "time: " << t4-t0 << " us (total)" << endl;

  return tensor;
}

// this function modifies coef_mat!
vector<double> utils::svd_fit(const System &sys, const rdcdata_t &fit_data, double coef_mat_array[]) {
//  struct timeval tv;
//  gettimeofday(&tv, NULL);
//  suseconds_t t0 = tv.tv_usec;

  const int n_rdc = fit_data.size();
  const int n_ah = 5;
  vector<double> tensor(n_ah, 0.0);

  // set up gsl for svd fitting
  gsl_vector *exp_reduced = gsl_vector_alloc(n_rdc);
  gsl_vector *S = gsl_vector_alloc(n_ah);
  gsl_vector *Stmp = gsl_vector_alloc(n_ah);
  gsl_vector *work = gsl_vector_alloc(n_ah);
  gsl_matrix *V = gsl_matrix_alloc(n_ah, n_ah);

//  gettimeofday(&tv, NULL);
//  suseconds_t t1 = tv.tv_usec;

  gsl_matrix_view coef_mat = gsl_matrix_view_array(coef_mat_array, n_rdc, n_ah);

  for(int k=0; k<n_rdc; k++){
    gsl_vector_set(exp_reduced, k, fit_data[k].exp / fit_data[k].dmax);
  }

//  gettimeofday(&tv, NULL);
//  suseconds_t t2 = tv.tv_usec;

  // calculate alignment tensor (solve SLE by SVD)
  // \mtx coef_mat * \vec S = \vec exp_reduced
  gsl_linalg_SV_decomp(&coef_mat.matrix, V, Stmp, work);
  gsl_linalg_SV_solve(&coef_mat.matrix, V, Stmp, exp_reduced, S);

//  gettimeofday(&tv, NULL);
//  suseconds_t t3 = tv.tv_usec;

  for(int i=0; i<n_ah; i++){
    tensor[i] = gsl_vector_get(S, i);
//    cout << gsl_vector_get(S, i) << endl;
  }

  gsl_vector_free(work);
  gsl_vector_free(exp_reduced);
  gsl_vector_free(Stmp);
  gsl_vector_free(S);
  gsl_matrix_free(V);

//  gettimeofday(&tv, NULL);
//  suseconds_t t4 = tv.tv_usec;
//  cout << "time: " << t1-t0 << " us (allocate)" << endl;
//  cout << "time: " << t2-t1 << " us (prepare)" << endl;
//  cout << "time: " << t3-t2 << " us (do stuff)" << endl;
//  cout << "time: " << t4-t3 << " us (write out, deallocate)" << endl;
//  cout << "time: " << t4-t0 << " us (total)" << endl;

  return tensor;
}

void utils::diagonalize_tensor(const std::vector<double> &t5, std::vector<double> &eval, std::vector<std::vector<double> >  &evec)
{
  eval.resize(3);
  evec.resize(3,vector<double>(3));

  // extract components
  const double Sxx = t5[0], Syy = t5[1], Szz = - Sxx - Syy, Sxy = t5[2], Sxz = t5[3], Syz = t5[4];

  // allocate
  gsl_matrix *evec_gsl = gsl_matrix_alloc(3, 3);
  gsl_vector *eval_gsl = gsl_vector_alloc(3);
  gsl_matrix *Smat = gsl_matrix_alloc(3, 3);

  gsl_matrix_set(Smat, 0, 0, Sxx);
  gsl_matrix_set(Smat, 0, 1, Sxy);
  gsl_matrix_set(Smat, 0, 2, Sxz);
  gsl_matrix_set(Smat, 1, 0, Sxy); // == Syx
  gsl_matrix_set(Smat, 1, 1, Syy);
  gsl_matrix_set(Smat, 1, 2, Syz);
  gsl_matrix_set(Smat, 2, 0, Sxz); // == Szx
  gsl_matrix_set(Smat, 2, 1, Syz); // == Szy
  gsl_matrix_set(Smat, 2, 2, Szz);
  // allocate space (O4n, where n is the size given (for nxn matrix))
  gsl_eigen_symmv_workspace *ework = gsl_eigen_symmv_alloc(3);
  // compute the eigenvalues and eigenvectors
  gsl_eigen_symmv(Smat, eval_gsl, evec_gsl, ework);
  // sort the absolute values of the eigenvalues in ascending order, so abs(z) is the largest and abs(x) the smallest (bax-2001)
  gsl_eigen_symmv_sort(eval_gsl, evec_gsl, GSL_EIGEN_SORT_ABS_ASC);

  // get the eigenvalues
  for(int i=0; i<3; i++){
    eval[i] = gsl_vector_get(eval_gsl, i);
    for(int j=0; j<3; j++){
      evec[i][j] = gsl_matrix_get(evec_gsl, i, j);
    } 
  }

  // free space
  gsl_eigen_symmv_free(ework);
  gsl_matrix_free(Smat);
  gsl_vector_free(eval_gsl);
  gsl_matrix_free(evec_gsl);
}

