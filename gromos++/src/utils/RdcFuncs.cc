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
        vector<RDCData::rdcparam> &rdcp) {

  // read into rdcparam class
  for (unsigned int jj = 1; jj < buffer.size() - 1; jj++) {

    // first get atom numbers
    istringstream is(buffer[jj]);
    int i, j;
    is >> i >> j;

    // adjust i and j to gromos numbering
    i--;
    j--;

    // define local set of RDC parameters   
    RDCData::rdcparam rp;

    // get molecule number, offset atoms numbers
    int m = 0, offset = 0;
    for (; i >= sys.mol(m).numAtoms() + offset; m++)
      offset += sys.mol(m).numAtoms();
    rp.mol = m;
    rp.i = i - offset;
    rp.j = j - offset;

    // now get the other parameters
    is >> rp.w >> rp.exp >> rp.gi >> rp.gj;

    // convert into gromos units
    double rdc_cf = pico;
    double g_cf = (atomic_mass_unit / elementary_charge) * 1.0e7;
    rp.exp *= rdc_cf;
    rp.gi *= g_cf;
    rp.gj *= g_cf;

    /*
    // DEBUG: convert into SI units
    double g_cf = 1.0e7;
    rp.gi *= g_cf;
    rp.gj *= g_cf;
     */

    if (is.fail())
      throw gromos::Exception("RdcFuncs::read_rdc", "Bad line in rdc-file\n" + buffer[jj]);

    rdcp.push_back(rp);
  }

}

// function to compute the prefactor, Dmax, with 8 * pi^3 * rij^3 as denominator

void RdcFuncs::calc_dmax_8pi3rij3(const System &sys, vector<RDCData::rdcparam> &rdcp,
        bool calc_rij) {

  int nrdc = rdcp.size();
  double rij;
  for (int d = 0; d < nrdc; d++) {

    if (!calc_rij) {

      // assign internuclear distances (from Rob Best) -- now scaled into gromos units (nm)
      int iaci = sys.mol(rdcp[d].mol).topology().atom(rdcp[d].i).iac();
      int iacj = sys.mol(rdcp[d].mol).topology().atom(rdcp[d].j).iac();

      // NOTE: stored IACs are one less than those in topology!!

      // N:H or H:N
      if (((iaci == 7) || (iaci == 5) && (iacj == 20)) || ((iaci == 20) && ((iacj == 7) || (iacj == 5)))) {
        //rij = 1.04e-10;
        rij = 0.104;
        // C:H(n) or H(n):C
      } else if (((iaci == 11) && (iacj == 20)) || ((iaci == 20) && (iacj == 11))) {
        //rij = 2.04e-10;
        rij = 0.204;
        // N:C or C:N
      } else if ((((iaci == 7) || (iaci == 5)) && (iacj == 11)) || ((iaci == 11) && (iacj == 7) || (iacj == 5))) {
        //rij = 1.33e-10;
        rij = 0.133;
        // CA:C or C:CA
      } else if (((iaci == 13) || (iaci == 14)) && (iacj == 11) || (iaci == 11) && ((iacj == 13) || (iacj == 14))) {
        //rij = 1.53e-10;
        rij = 0.153;
      } else {
        throw gromos::Exception("RdcFuncs::read_rdc", "Can't find IAC for RDC " + d + 1);
      }

    }
      // DEBUG: compute and store Dmax using SI mu0 and hbar
      //double SImu0 = 4.*pi*1e-7;
      //double SIhbar = 1.05e-34;
      // note now 4.0*pi^2 due to hbar
      //rdcp[d].dmax = (SImu0 * SIhbar * rdcp[d].gi * rdcp[d].gj) / (4.0 * pi * pi * rij * rij * rij);

    else {

      // get the distance rij
      rij = (sys.mol(rdcp[d].mol).pos(rdcp[d].i) - sys.mol(rdcp[d].mol).pos(rdcp[d].j)).abs();

    }
    // store rij
    rdcp[d].rij = rij;

    // compute and store Dmax
    rdcp[d].dmax = (mu0 * h * rdcp[d].gi * rdcp[d].gj) /
            (8.0 * pi * pi * pi * rij * rij * rij);

  }
}

// function to compute the prefactor, Dmax, with 16 * pi^3 * rij^5 as denominator
void RdcFuncs::calc_dmax_16pi3rij3(const System &sys, vector<RDCData::rdcparam> &rdcp) {

  int nrdc = rdcp.size();
  for (int d = 0; d < nrdc; d++) {

    // get the distance rij
    double rij = (sys.mol(rdcp[d].mol).pos(rdcp[d].i) - sys.mol(rdcp[d].mol).pos(rdcp[d].j)).abs();

    // store rij
    rdcp[d].rij = rij;

    // compute and store Dmax
    rdcp[d].dmax = (mu0 * h * rdcp[d].gi * rdcp[d].gj) /
            (16.0 * pi * pi * pi * rij * rij * rij);

  }
}

// function to compute the prefactor, Dmax, for an arbitrary NH vector
double RdcFuncs::calc_dmax_NH() {

  // hard-wired NH distance
  double rNH = 0.104;
  // hard-wired gyromagnetic ratios (note gN is positive)
  double gN = 2.712 * (atomic_mass_unit / elementary_charge) * 1.0e7;
  double gH = 26.7520 * (atomic_mass_unit / elementary_charge) * 1.0e7;
  // compute dmaxNH
  double dmax = ( mu0 * h * gN * gH ) / (8.0 * pi * pi * pi * rNH * rNH * rNH);
  return dmax;
  
}

// function to read weights for individual frames from file
void RdcFuncs::read_weights(vector<string> buffer, vector<RDCWeights::weights> &weight_data) {

  // read into weights vector
  for (int d = 1; d < buffer.size() - 1; d++) {

    // define local weight struct
    RDCWeights::weights w;
    // get values from file
    istringstream is(buffer[d]);
    is >> w.frame >> w.weight;
    // add to weight_data
    weight_data.push_back(w);

  }
}

// compute the coefficients of the matrix describing bond vector fluctuations
void RdcFuncs::calc_coef(const System &sys, vector<RDCData::rdcparam> &fit_data,
        gsl_matrix *coef_mat, int nrdc, double w)
{
  int i, j, m;
  double mu_x, mu_y, mu_z, mu_r;

  for (int d = 0; d < nrdc; d++) {

    // get atom numbers
    i = fit_data[d].i;
    j = fit_data[d].j;
    m = fit_data[d].mol;

    // get diff in coords
    mu_x = sys.mol(m).pos(i)[0] - sys.mol(m).pos(j)[0];
    mu_y = sys.mol(m).pos(i)[1] - sys.mol(m).pos(j)[1];
    mu_z = sys.mol(m).pos(i)[2] - sys.mol(m).pos(j)[2];

    mu_r = sqrt(mu_x * mu_x + mu_y * mu_y + mu_z * mu_z);
    mu_x /= mu_r;
    mu_y /= mu_r;
    mu_z /= mu_r;

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

// fill a gsl vector with normalised RDCs (divided by Dmax)
void RdcFuncs::fill_rdcvec_norm(const vector<RDCData::rdcparam> &R, gsl_vector *v) {

  int nrdc = R.size();
  for (int i = 0; i < nrdc; i++) {
    if (R[i].dmax != 0) {
      gsl_vector_set(v, i, R[i].exp / R[i].dmax);
    } else {
      throw gromos::Exception("RdcFuncs::fill_rdcvec_norm",
              "Dmax not calculated yet\n");
    }
  }
}

// unnormalise RDCs (i.e. multiply by Dmax)
void RdcFuncs::unnorm_rdc(const vector<RDCData::rdcparam> &R, gsl_vector *v, gsl_vector *w) {

  int nrdc = R.size();
  double tmp;
  for (int i = 0; i < nrdc; i++) {
    if (R[i].dmax != 0) {
      tmp = gsl_vector_get(v, i);
      gsl_vector_set(w, i, tmp * R[i].dmax);
    } else {
      throw gromos::Exception("RdcFuncs::unnorm_rdc",
              "Dmax not calculated yet\n");
    }
  }
}

// calculate the Q value (goodness of fit)
double RdcFuncs::calc_Q(gsl_vector *calc, gsl_vector *expt)
{
	int ndat;
	double sdev2,scalc2,tmpf,expt_i,calc_i;
	sdev2 = 0;
	scalc2 = 0;
	ndat = calc->size;
	for (int i=0; i<ndat; i++) {
		calc_i = gsl_vector_get(calc,i);
		expt_i = gsl_vector_get(expt,i);
		tmpf = calc_i-expt_i;
		sdev2 += tmpf*tmpf;
		scalc2 += calc_i*calc_i;
	}
	return sqrt(sdev2/scalc2);
}

// function to compute the angle between the magnetic field direction
// and the internuclear vector, then the rdc
void RdcFuncs::calc_rdc_H(const gcore::System &sys,
        std::vector<utils::RDCData::rdcparam> &rdcp,
        double theta, double phi, gsl_vector *rdc_tmp) {

  int nrdc = rdcp.size();
  for (int d = 0; d < nrdc; d++) {

    // get the vector rij
    Vec rij = (sys.mol(rdcp[d].mol).pos(rdcp[d].i) - sys.mol(rdcp[d].mol).pos(rdcp[d].j));

    // get length (note that this is stored as rij in rdcparam (for ref))
    double length = rij.abs();

    // convert H(theta,phi) to H(x,y,z) (length H = 1.0)
    Vec H = (sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta));

    // compute angle between H and rij, normalise by length (length H = 1.0)
    double cos_theta = rij.dot(H) / length;

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
void RdcFuncs::rdcYlm(int nrdc, gsl_vector *rdc_tmp,
        unsigned int this_sh, double ylm, gsl_matrix *rdc_ylm) {

  for (int d = 0; d < nrdc; d++) {

    // multiply by spherical harmonic
    double rdc_ylm_tmp = gsl_vector_get(rdc_tmp, d) * ylm;

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
