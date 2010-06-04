// utils_RdcFuncs.cc

#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <cassert>
#include <sstream>
#include <vector>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_eigen.h>

#include "RdcFuncs.h"
#include "../args/Arguments.h"
#include "../fit/Reference.h"
#include "../fit/PositionUtils.h"
#include "../gcore/System.h"
#include "../gcore/Molecule.h"
#include "../gcore/MoleculeTopology.h" // added for DEBUG
#include "../gcore/AtomTopology.h" // added for DEBUG
#include "../gmath/Matrix.h"
#include "../gmath/Physics.h"
#include "../gmath/Vec.h"
#include "../gromos/Exception.h"
#include "../gio/Ginstream.h"

using namespace std;
using namespace fit;
using namespace gio;
using namespace gcore;
using namespace gmath;

using gcore::System;
//using fit::Reference;
using args::Arguments;
using utils::RdcFuncs;
using gcore::MoleculeTopology; // added for DEBUG
using gcore::AtomTopology; // added for DEBUG

// Constructor
RdcFuncs::RdcFuncs(System &sys, Arguments &args) {}

// Destructor
RdcFuncs::~RdcFuncs() {}

// function to read RDC data
void RdcFuncs::read_rdc(vector<string> buffer, const System &sys,
        vector<RDCData::rdcparam> &rdcp, bool calc_rij, bool fit) {

  // read into rdcparam class
  for (unsigned int jj = 1; jj < buffer.size() - 1; jj++) {

    // first get atom numbers
    istringstream is(buffer[jj]);
    unsigned int i, j, k;
    is >> i >> j >> k;
    // adjust to gromos numbering
    i--;
    j--;
    if (i < 0 || j < 0) {
      throw gromos::Exception("RdcFuncs::read_rdc",
              "Disallowed atom number in rdc file\n" + buffer[jj]);
    }

    // define local set of RDC parameters   
    RDCData::rdcparam rp;

    // get molecule number, offset atom numbers
    unsigned int m = 0, offset = 0;
    for (; i >= sys.mol(m).numAtoms() + offset; m++)
      offset += sys.mol(m).numAtoms();
    rp.mol = m;
    rp.i = i - offset;
    rp.j = j - offset;
    if (rp.i < 0 || rp.j < 0) {
      throw gromos::Exception("RdcFuncs::read_rdc",
              "Offsetting for molecule number produces negative atom number\n");
    }

    // now get the other parameters
    is >> rp.w >> rp.exp >> rp.gi >> rp.gj >> rp.rij >> rp.rik >> rp.type;

    // we can't fit to HH RDCs
    if (rp.type < 0 && fit) {
      throw gromos::Exception("RdcFuncs::read_rdc",
              "You cannot fit to HH RDCs because the sign is not known\n" + buffer[jj]);
    }

    // if we have side-chain NH RDCs, adjust and check atom k
    else if (rp.type > 0) {

      // we can't fit to NH RDCs
      if (fit) {
        throw gromos::Exception("RdcFuncs::read_rdc",
                "You cannot fit to side-chain NH RDCs because they are a sum\n" + buffer[jj]);
      } else {

        k--;
        if (k < 0) {
          throw gromos::Exception("RdcFuncs::read_rdc",
                  "Disallowed atom number in rdc file\n" + buffer[jj]);
        }
        rp.k = k - offset;
        if (rp.k < 0) {
          throw gromos::Exception("RdcFuncs::read_rdc",
                  "Offsetting for molecule number produces negative atom number\n");
        }
        // check that type(atom j) = type(atom k) = H
        if (!sys.mol(m).topology().atom(rp.j).isH() || !sys.mol(m).topology().atom(rp.k).isH() ||
                (sys.mol(m).topology().atom(rp.j).iac() != sys.mol(m).topology().atom(rp.k).iac())) {
          throw gromos::Exception("RdcFuncs::read_rdc",
                  "Atoms j and k must be of the same type (hydrogen) \n");
        }
      }
    }

    // convert into gromos units
    const double rdc_cf = gmath::physConst.get_pico();
    const double atomic_mass_unit = gmath::physConst.get_atomic_mass_unit();
    const double elementary_charge = gmath::physConst.get_elementary_charge();
    const double g_cf = (atomic_mass_unit / elementary_charge) * 1.0e7;
    rp.exp *= rdc_cf;
    rp.gi *= g_cf;
    rp.gj *= g_cf;

    if (!calc_rij) {
      if (rp.rij <= 0.0) {
        throw gromos::Exception("RdcFuncs::read_rdc",
                "If using rij from rdc file, it must be positive and non-zero\n");
      }
      else if (rp.type > 0 && rp.rik <= 0.0) {
        throw gromos::Exception("RdcFuncs::read_rdc",
                "If using rik from rdc file, it must be positive and non-zero\n");
      }
    }

    /*
    // DEBUG: convert into SI units
    double g_cf = 1.0e7;
    rp.gi *= g_cf;
    rp.gj *= g_cf;
     */

    if (is.fail())
      throw gromos::Exception("RdcFuncs::read_rdc", "Bad line in rdc-file\n" + buffer[jj]);

    // store temporary RDC parameters
    rdcp.push_back(rp);
  }

}

// function to compute the prefactor, Dmax, with 8 * pi^3 * rij^3 as denominator
void RdcFuncs::calc_dmax_8pi3rij3(const System &sys, vector<RDCData::rdcparam> &rdcp,
        bool calc_rij) {

  const double pi = gmath::physConst.get_pi();
  const double pi3 = pi * pi * pi;
  const double mu0 = gmath::physConst.get_mu0();
  const double h = gmath::physConst.get_h();
  const unsigned int nrdc = rdcp.size();
  for (unsigned int d = 0; d < nrdc; d++) {

    if (calc_rij) {

      // calculate the inter-nuclear distances from the reference coordinates
      double rij;
      // get the distances rij
      rij = (sys.mol(rdcp[d].mol).pos(rdcp[d].i) - sys.mol(rdcp[d].mol).pos(rdcp[d].j)).abs();
      // alter the stored rij (read from file)
      rdcp[d].rij = rij;
      
      // and if it's an NH sidechain, calculate rik too
      if (rdcp[d].type > 0) {
        double rik;
        // get the distances rij
        rik = (sys.mol(rdcp[d].mol).pos(rdcp[d].i) - sys.mol(rdcp[d].mol).pos(rdcp[d].k)).abs();
        // alter the stored rij (read from file)
        rdcp[d].rik = rij;
      }
    }

    // compute and store Dmax (should be the same for ij and jk as we only deal with NH side-chains
    rdcp[d].dmax = (mu0 * h * rdcp[d].gi * rdcp[d].gj) /
            (8.0 * pi3 * rdcp[d].rij * rdcp[d].rij * rdcp[d].rij);
  }
}

// function to compute the prefactor, Dmax, with 16 * pi^3 * rij^5 as denominator
void RdcFuncs::calc_dmax_16pi3rij3(const System &sys, vector<RDCData::rdcparam> &rdcp) {

  const unsigned int nrdc = rdcp.size();
  const double pi = gmath::physConst.get_pi();
  const double pi3 = pi * pi * pi;
  const double mu0 = gmath::physConst.get_mu0();
  const double h = gmath::physConst.get_h();
  for (unsigned int d = 0; d < nrdc; d++) {
    // get the distance rij
    double rij = (sys.mol(rdcp[d].mol).pos(rdcp[d].i) - sys.mol(rdcp[d].mol).pos(rdcp[d].j)).abs();
    // store rij
    rdcp[d].rij = rij;
    // compute and store Dmax
    rdcp[d].dmax = (mu0 * h * rdcp[d].gi * rdcp[d].gj) /
            (16.0 * pi3 * rij * rij * rij);
  }
}

// function to read weights for individual frames from file
void RdcFuncs::read_weights(vector<string> buffer, vector<RDCWeights::weights> &weight_data) {

  // read into weights vector
  for (unsigned int d = 1; d < buffer.size() - 1; d++) {

    // define local weight struct
    RDCWeights::weights w;
    // get values from file
    istringstream is(buffer[d]);
    is >> w.frame >> w.weight;
    // add to weight_data
    weight_data.push_back(w);

  }
}

// compute the coefficients of the matrix describing bond vector fluctuations for fit RDCs
void RdcFuncs::calc_coef_fit(const System &sys, vector<RDCData::rdcparam> &fit_data,
        gsl_matrix *coef_mat, unsigned int nrdc, double w)
{

  for (unsigned int d = 0; d < nrdc; d++) {

    // get atom numbers
    unsigned int i = fit_data[d].i;
    unsigned int j = fit_data[d].j;
    unsigned int m = fit_data[d].mol;

    // get diff in coords
    double mu_x = sys.mol(m).pos(i)[0] - sys.mol(m).pos(j)[0];
    double mu_y = sys.mol(m).pos(i)[1] - sys.mol(m).pos(j)[1];
    double mu_z = sys.mol(m).pos(i)[2] - sys.mol(m).pos(j)[2];

    double mu_r = sqrt(mu_x * mu_x + mu_y * mu_y + mu_z * mu_z);
    mu_x /= mu_r;
    mu_y /= mu_r;
    mu_z /= mu_r;

    // put into gsl matrix
    gsl_matrix_set(coef_mat, d, 0,
            gsl_matrix_get(coef_mat, d, 0) + w * (mu_x * mu_x - mu_z * mu_z));
    gsl_matrix_set(coef_mat, d, 1,
            gsl_matrix_get(coef_mat, d, 1) + w * (mu_y * mu_y - mu_z * mu_z));
    gsl_matrix_set(coef_mat, d, 2,
            gsl_matrix_get(coef_mat, d, 2) + w * (2.0 * mu_x * mu_y));
    gsl_matrix_set(coef_mat, d, 3,
            gsl_matrix_get(coef_mat, d, 3) + w * (2.0 * mu_x * mu_z));
    gsl_matrix_set(coef_mat, d, 4,
            gsl_matrix_get(coef_mat, d, 4) + w * (2.0 * mu_y * mu_z));

  }
}

// compute the coefficients of the matrix describing bond vector fluctuations for back-calculated RDCs
void RdcFuncs::calc_coef_bc(const System &sys, vector<RDCData::rdcparam> &bc_data,
        gsl_matrix *coef_mat_j, gsl_matrix *coef_mat_k, unsigned int nrdc, double w) {

  for (unsigned int d = 0; d < nrdc; d++) {

    // get atom numbers
    unsigned int i = bc_data[d].i;
    unsigned int j = bc_data[d].j;
    unsigned int m = bc_data[d].mol;

    // get diff in coords for ij
    double mu_x_j = sys.mol(m).pos(i)[0] - sys.mol(m).pos(j)[0];
    double mu_y_j = sys.mol(m).pos(i)[1] - sys.mol(m).pos(j)[1];
    double mu_z_j = sys.mol(m).pos(i)[2] - sys.mol(m).pos(j)[2];

    double mu_r_j = sqrt(mu_x_j * mu_x_j + mu_y_j * mu_y_j + mu_z_j * mu_z_j);
    mu_x_j /= mu_r_j;
    mu_y_j /= mu_r_j;
    mu_z_j /= mu_r_j;

    // put into gsl matrix for ij
    gsl_matrix_set(coef_mat_j, d, 0,
            gsl_matrix_get(coef_mat_j, d, 0) + w * (mu_x_j * mu_x_j - mu_z_j * mu_z_j));
    gsl_matrix_set(coef_mat_j, d, 1,
            gsl_matrix_get(coef_mat_j, d, 1) + w * (mu_y_j * mu_y_j - mu_z_j * mu_z_j));
    gsl_matrix_set(coef_mat_j, d, 2,
            gsl_matrix_get(coef_mat_j, d, 2) + w * (2.0 * mu_x_j * mu_y_j));
    gsl_matrix_set(coef_mat_j, d, 3,
            gsl_matrix_get(coef_mat_j, d, 3) + w * (2.0 * mu_x_j * mu_z_j));
    gsl_matrix_set(coef_mat_j, d, 4,
            gsl_matrix_get(coef_mat_j, d, 4) + w * (2.0 * mu_y_j * mu_z_j));

    // if side-chain NH, do it all again for ik
    if (bc_data[d].type > 0) {

      unsigned int k = bc_data[d].k;

      // get diff in coords for ik
      double mu_x_k = sys.mol(m).pos(i)[0] - sys.mol(m).pos(k)[0];
      double mu_y_k = sys.mol(m).pos(i)[1] - sys.mol(m).pos(k)[1];
      double mu_z_k = sys.mol(m).pos(i)[2] - sys.mol(m).pos(k)[2];

      double mu_r_k = sqrt(mu_x_k * mu_x_k + mu_y_k * mu_y_k + mu_z_k * mu_z_k);
      mu_x_k /= mu_r_k;
      mu_y_k /= mu_r_k;
      mu_z_k /= mu_r_k;

      // put into gsl matrix for ik
      gsl_matrix_set(coef_mat_k, d, 0,
              gsl_matrix_get(coef_mat_k, d, 0) + w * (mu_x_k * mu_x_k - mu_z_k * mu_z_k));
      gsl_matrix_set(coef_mat_k, d, 1,
              gsl_matrix_get(coef_mat_k, d, 1) + w * (mu_y_k * mu_y_k - mu_z_k * mu_z_k));
      gsl_matrix_set(coef_mat_k, d, 2,
              gsl_matrix_get(coef_mat_k, d, 2) + w * (2.0 * mu_x_k * mu_y_k));
      gsl_matrix_set(coef_mat_k, d, 3,
              gsl_matrix_get(coef_mat_k, d, 3) + w * (2.0 * mu_x_k * mu_z_k));
      gsl_matrix_set(coef_mat_k, d, 4,
              gsl_matrix_get(coef_mat_k, d, 4) + w * (2.0 * mu_y_k * mu_z_k));

    }
  }
}

// fill a gsl vector with normalised RDCs (divided by Dmax)
void RdcFuncs::fill_rdcvec_norm(const vector<RDCData::rdcparam> &R, gsl_vector *v) {

  unsigned int nrdc = R.size();
  for (unsigned int i = 0; i < nrdc; i++) {
    if (R[i].dmax != 0) {
      gsl_vector_set(v, i, R[i].exp / R[i].dmax);
    } else {
      throw gromos::Exception("RdcFuncs::fill_rdcvec_norm",
              "Dmax not calculated yet\n");
    }
  }
}

// compute Euler angles according to x-convention (z,x',z'')
void RdcFuncs::euler_x(gmath::Matrix &A, double &alpha1, double &alpha2, double &beta1,
        double &beta2, double &gamma1, double &gamma2) {

  // -- "x-convention": Euler angles for rotation about z,x',z''
  const double pi = gmath::physConst.get_pi();
  // beta (about x')
  const double a33 = A(2, 2);
  if (abs(a33) != 1.0) {
    beta1 = acos(a33);
    beta2 = -beta1;
    // gamma (about z'')
    const double a13 = A(0, 2);
    const double a23 = A(1, 2);
    const double sign1 = sin(beta1);
    const double sign2 = sin(beta2);
    gamma1 = atan2(a13 / sign1, a23 / sign1);
    gamma2 = atan2(a13 / sign2, a23 / sign2);
    // alpha (about z)
    const double a31 = -A(2, 0);
    const double a32 = A(2, 1);
    alpha1 = atan2(a31 / sign1, a32 / sign1);
    alpha2 = atan2(a31 / sign2, a32 / sign2);
  } else {
    // special case: a33 = +/- 1 -> sin(theta) = 0
    alpha1 = alpha2 = gamma2 = 0; // arbitrary
    const double a11 = A(0, 0);
    const double a12 = A(0, 1);
    if (a33 == -1) {
      beta1 = pi;
      beta2 = -pi;
      gamma1 = alpha1 - atan2(a12, a11);
    } else {
      beta1 = 0.0;
      beta2 = pi * 2.0;
      gamma1 = atan2(a12, a11) - alpha1;
    }
  }
  const double radian2degree = gmath::physConst.get_radian2degree();
  alpha1 *= radian2degree;
  alpha2 *= radian2degree;
  beta1 *= radian2degree;
  beta2 *= radian2degree;
  gamma1 *= radian2degree;
  gamma2 *= radian2degree;
}

// compute Euler angles according to y-convention (z,y',z'')
void RdcFuncs::euler_y(gmath::Matrix &A, double &alpha1, double &alpha2, double &beta1,
        double &beta2, double &gamma1, double &gamma2) {

  // -- "y-convention": Euler angles for rotation about z,y',z''
  const double pi = gmath::physConst.get_pi();
  // beta (about y')
  const double a33 = A(2, 2);
  if (abs(a33) != 1.0) {
    beta1 = acos(a33);
    beta2 = -beta1;
    // alpha (about z)
    const double a32 = A(2, 1);
    const double a31 = A(2, 0);
    const double sign1 = sin(beta1);
    const double sign2 = sin(beta2);
    alpha1 = atan2(a32 / sign1, a31 / sign1);
    alpha2 = atan2(a32 / sign2, a31 / sign2);
    // gamma (about z'')
    const double a23 = -A(1, 2);
    const double a13 = A(0, 2);
    gamma1 = atan2(a23 / sign1, a13 / sign1);
    gamma2 = atan2(a23 / sign2, a13 / sign2);
  } else {
    // special case: a33 = +/- 1 -> sin(theta) = 0
    alpha1 = alpha2 = gamma2 = 0; // arbitrary
    const double a11 = A(0, 0);
    const double a12 = A(0, 1);
    if (a33 == -1) {
      beta1 = pi;
      beta2 = -pi;
      gamma1 = atan2(a12, -a11) + alpha1;
    } else {
      beta1 = 0.0;
      beta2 = pi * 2.0;
      gamma1 = atan2(a12, a11) - alpha1;
    }
  }
  const double radian2degree = gmath::physConst.get_radian2degree();
  alpha1 *= radian2degree;
  alpha2 *= radian2degree;
  beta1 *= radian2degree;
  beta2 *= radian2degree;
  gamma1 *= radian2degree;
  gamma2 *= radian2degree;
}

// compute Euler angles according to pitch-roll-yaw convention (x,y,z)
void RdcFuncs::euler_pry(gmath::Matrix &A, double &psi1, double &psi2, double &theta1,
        double &theta2, double &phi1, double &phi2) {

  // -- "pitch-roll-yaw-convention": Euler angles for rotation about x,y,z
  const double pi = gmath::physConst.get_pi();
  // theta (about y axis)
  const double a31 = A(2, 0);
  if (abs(a31) != 1.0) {
    theta1 = -asin(a31);
    theta2 = pi - theta1;
    // psi (about x-axis)
    const double a32 = A(2, 1);
    const double a33 = A(2, 2);
    const double sign1 = cos(theta1);
    const double sign2 = cos(theta2);
    psi1 = atan2(a32 / sign1, a33 / sign1);
    psi2 = atan2(a32 / sign2, a33 / sign2);
    // phi (about z-axis)
    const double a21 = A(1, 0);
    const double a11 = A(0, 0);
    phi1 = atan2(a21 / sign1, a11 / sign1);
    phi2 = atan2(a21 / sign2, a11 / sign2);
  } else {
    // special case: a31 = +/- 1 -> cos(theta) = 0
    phi1 = phi2 = psi2 = theta2 = 0; // arbitrary
    const double a12 = A(0, 1);
    const double a13 = A(0, 2);
    if (a31 == -1) {
      theta1 = pi * 0.5;
      psi1 = phi1 + atan2(a12, a13);
    } else {
      theta1 = pi * -0.5;
      psi1 = -phi1 + atan2(-a12, -a13);
    }
  }
  // convert to degrees
  const double radian2degree = gmath::physConst.get_radian2degree();
  psi1 *= radian2degree;
  psi2 *= radian2degree;
  theta1 *= radian2degree;
  theta2 *= radian2degree;
  phi1 *= radian2degree;
  phi2 *= radian2degree;
}

// sum the back-calculated RDCs
void RdcFuncs::sum_rdc(unsigned int nrdc, gsl_vector *bc_j, gsl_vector *bc_k, gsl_vector *bc_sum) {

  for (unsigned int i = 0; i < nrdc; i++) {
    gsl_vector_set(bc_sum, i, gsl_vector_get(bc_j, i) + gsl_vector_get(bc_k, i));
  }
}

// unnormalise RDCs (i.e. multiply by Dmax)
void RdcFuncs::unnorm_rdc(const vector<RDCData::rdcparam> &R, gsl_vector *v, gsl_vector *w) {

  unsigned int nrdc = R.size();
  double tmp;
  for (unsigned int i = 0; i < nrdc; i++) {
    if (R[i].dmax != 0.0) {
      tmp = gsl_vector_get(v, i) * R[i].dmax;
      gsl_vector_set(w, i, tmp);
    } else {
      throw gromos::Exception("RdcFuncs::unnorm_rdc",
              "Dmaxij not calculated yet\n");
    }
  }
}

// calculate the Q value (goodness of fit) Rob Best's way
double RdcFuncs::calc_Q1(gsl_vector *calc, gsl_vector *expt)
{
	double sdev2 = 0;
	double scalc2 = 0;
	unsigned int ndat = calc->size;
	for (unsigned int i=0; i<ndat; i++) {
		double calc_i = gsl_vector_get(calc,i);
		double expt_i = gsl_vector_get(expt,i);
		double tmpf = calc_i-expt_i;
		sdev2 += tmpf*tmpf;
		scalc2 += calc_i*calc_i;
	}
	return sqrt(sdev2/scalc2);
}

// calculate the Q value according to Cornilescu, Alexandrescu, Bax, Zweckstetter...
double RdcFuncs::calc_Q2(gsl_vector *calc, gsl_vector *expt)
{
	double sdev2 = 0;
	double sexpt2 = 0;
	const unsigned int ndat = calc->size;
        const double n = ndat;
	for (unsigned int i=0; i<ndat; i++) {
		double calc_i = gsl_vector_get(calc,i);
		double expt_i = gsl_vector_get(expt,i);
		double tmpf = expt_i-calc_i;
		sdev2 += tmpf*tmpf;
		sexpt2 += expt_i*expt_i;
	}
        sdev2 /= n;
        sexpt2 /= n;
        double sdev = sqrt(sdev2);
        double sexpt = sqrt(sexpt2);
	return sdev/sexpt;
}

// function to compute the angle between the magnetic field direction
// and the internuclear vector, then the rdc
void RdcFuncs::calc_rdc_H(const gcore::System &sys,
        std::vector<utils::RDCData::rdcparam> &rdcp,
        double theta, double phi, gsl_vector *rdc_tmp) {

  unsigned int nrdc = rdcp.size();
  for (unsigned int d = 0; d < nrdc; d++) {

    // get the vector rij
    Vec rij = (sys.mol(rdcp[d].mol).pos(rdcp[d].i) - sys.mol(rdcp[d].mol).pos(rdcp[d].j));

    // get length (note that this is stored as rij in rdcparam (for ref))
    const double length = rij.abs();

    // convert H(theta,phi) to H(x,y,z) (length H = 1.0)
    Vec H = gmath::Vec(sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta));

    // compute angle between H and rij, normalise by length (length H = 1.0)
    const double cos_theta = rij.dot(H) / length;

    // compute RDC (without prefactor, as experimental one has been scaled to remove prefactor)
    // note denominator of 2 is in prefactor
    double rdc = 3. * cos_theta * cos_theta - 1;

    // DEBUG: "un"scale
    rdc *= rdcp[d].dmax;

    // store temporarily
    gsl_vector_set(rdc_tmp, d, rdc);

  }
}

// function to multiply each RDC by Y^l_m and store in the matrix for solving
void RdcFuncs::rdcYlm(unsigned int nrdc, gsl_vector *rdc_tmp,
        unsigned int this_sh, double ylm, gsl_matrix *rdc_ylm) {

  for (unsigned int d = 0; d < nrdc; d++) {

    // multiply by spherical harmonic
    const double rdc_ylm_tmp = gsl_vector_get(rdc_tmp, d) * ylm;

    // add to matrix rdc_ylm[nrdc, nsh]
    gsl_matrix_set(rdc_ylm, d, this_sh, gsl_matrix_get(rdc_ylm, d, this_sh) + rdc_ylm_tmp);

  }

}

// function to rotate two structures onto each other (without performing translation)
void RdcFuncs::rot_fit(gcore::System &sys, const fit::Reference &ref) {

  Matrix rot(3, 3);
  rotationMatrix(&rot, sys, ref);
  PositionUtils::rotate(&sys, rot);
}

// function to compute the rotation matrix
void RdcFuncs::rotationMatrix(gmath::Matrix *mat, const gcore::System &sys,
        const fit::Reference &r) {

  const System &ref = r.sys();

  Matrix U(3,3,0);
  for(int i=0;i<3;++i)
    for(int j=0;j<3;++j)
      for(int m=0;m<ref.numMolecules();++m)
  	for(int n=0;n<ref.mol(m).numAtoms();++n)
  	  if(r.weight(m,n))
	    U(i,j)+=r.weight(m,n)*sys.mol(m).pos(n)[i]*ref.mol(m).pos(n)[j];

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

  if(det<0 && fabs(eigenvals[1] - eigenvals[2]) < 1.0e-5){

    std::cerr << "determinant = " << det << "\n"
	      << "eigenval[0] = " << eigenvals[0] << "\n"
	      << "eigenval[1] = " << eigenvals[1] << "\n"
	      << "eigenval[2] = " << eigenvals[2] << "\n" << std::endl;

    throw RdcFuncs::Exception("Rotation matrix degenerate!");
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

  *mat = Matrix(h1,k1) + Matrix(h2,k2) + signU*Matrix(h3,k3);

  delete[] eigenvals;

}
