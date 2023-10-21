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

// pb_FDPoissonBoltzmann.cc

#include <new>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <cassert>
#include <set>
#include "../fit/PositionUtils.h"
#include "../utils/AtomSpecifier.h"
#include "../gmath/Physics.h"
#include "../gcore/System.h"

#include "FDPoissonBoltzmann.h"
#include "PB_Parameters.h"

#include "FDPoissonBoltzmann_ICCG_PBC.h"
#include "FDPoissonBoltzmann_ICCG_NPBC.h"

using pb::FDPoissonBoltzmann;
using pb::PB_Parameters;



FDPoissonBoltzmann::FDPoissonBoltzmann(utils::AtomSpecifier atoms,utils::AtomSpecifier atoms_to_charge,
				       int gridpointsX, int gridpointsY, int gridpointsZ, double gridspace,
				       bool pbc, double epssolvent, ofstream &os):ppp(epssolvent, os){
  this->atoms=atoms;
  this->atoms_to_charge=atoms_to_charge;
  this->pbc=pbc;
  this->epssolvent=epssolvent;
  this->epssolute=ppp.getEpssolute();

  // set the grid point variables
  GPX=gridpointsX;
  GPY=gridpointsY;
  GPZ=gridpointsZ;
  gridspacing=gridspace;
  os << "# gridspacing " << gridspacing << endl;
  //  if (pbc == false) {
  //  increaserid(os);
  //}
  GPXGPY=GPX*GPY;
  GPXGPYGPZ=GPX*GPY*GPZ;

  // read the radii
  for (unsigned int i=0; i<atoms.size(); i++){
    this->radii.push_back(atoms.radius(i));
  }


  os << "# GPX " << GPX << endl;
  os << "# GPY " << GPY << endl;
  os << "# GPZ " << GPZ << endl;
  os << "# epssolvent " << epssolvent << endl;
  for (unsigned int i=0; i<atoms_to_charge.size(); i++){
    os << "# atom " << i+1 << " charge " << atoms_to_charge.charge(i) << endl;
    os << "# atom " << i+1 << " radius " << atoms_to_charge.radius(i) << endl;
  }
  
  
  phigrid.resize(GPXGPYGPZ, 0.0);
  rhogrid.resize(GPXGPYGPZ, 0.0);
  epsIgrid.resize(GPXGPYGPZ, 0.0);
  epsJgrid.resize(GPXGPYGPZ, 0.0);
  epsKgrid.resize(GPXGPYGPZ, 0.0);
  epsCgrid.resize(GPXGPYGPZ, 0.0);
}

void FDPoissonBoltzmann::setupGrid(bool newphi, ofstream &os, double gridstartX, double gridstartY, double gridstartZ, double gridcenterX, double gridcenterY, double gridcenterZ){
  
  //the boolean determines whether to use a previous solution --
  //if not, we use the previous phigrid as well as the calculated grid
  //dimensions (anything else does not make much sense)
  //the grid dimension check at the end is still done, though
  if (newphi){
    // zero all the grids
    for (int i=0;  i<GPXGPYGPZ; i++){
      
      phigrid[i]=0.0;
      rhogrid[i]=0.0;
    }
  }
  else
    {
      for (int i=0;  i<GPXGPYGPZ; i++){
	rhogrid[i]=0.0;
      }
    }// end of if newphi

  
  
  // gridcenter
  //double gridcenter[3];
  //it should be the cog of the atoms
  //UUU
  if (gridcenterX == 0 && gridcenterY == 0 && gridcenterZ == 0) {
  //   gmath::Vec gridcenter = fit::PositionUtils::cog(atoms.sys(),atoms);
    gmath::Vec coormin = fit::PositionUtils::getmincoordinates(atoms.sys(), false);
    gmath::Vec coormax = fit::PositionUtils::getmaxcoordinates(atoms.sys(), false);
    os << "# Calculating new gridcenter" << endl;
    os << "# mincoordinates of all atoms: " << coormin[0] << " " << coormin[1] << " " << coormin[2] << endl;
    os << "# maxcoordinates of all atoms: " << coormax[0] << " " << coormax[1] << " " << coormax[2] << endl;  
  
    //gridcenter[0]=fit::PositionUtils::cog(*atoms.sys(), atoms)[0];
    //gridcenter[1]=fit::PositionUtils::cog(*atoms.sys(), atoms)[1];
    //gridcenter[2]=fit::PositionUtils::cog(*atoms.sys(), atoms)[2];
    for (int i=0;  i<3; i++)
      gridcenter[i] = (coormin[i] + coormax[i]) * 0.5;

    // os << "# cog x,y,z = " << fit::PositionUtils::cog(*atoms.sys(), atoms)[0] << " "<< \
	 // fit::PositionUtils::cog(*atoms.sys(), atoms)[1] << " " << fit::PositionUtils::cog(*atoms.sys(), atoms)[2] << endl;
  } else {
    gridcenter[0] = gridcenterX;
    gridcenter[1] = gridcenterY;
    gridcenter[2] = gridcenterZ;
  }
  
  os << "# GPX " << GPX << endl;
  os << "# GPY " << GPY << endl;
  os << "# GPZ " << GPZ << endl;
  
  if (gridstartX == 0 && gridstartY == 0 && gridstartZ == 0) {
    // build the grid around the cog;
    // set the boundary by adding gridspacing/2 to the last grid point
    os << "# Calculating new gridstart" << endl;
    gridstart[0]=gridcenter[0]-0.5*(gridspacing*GPX)-gridspacing/2.0;
    os << "# gridstart[0] " << gridstart[0] << endl;
    gridstart[1]=gridcenter[1]-0.5*(gridspacing*GPY)-gridspacing/2.0;
    os << "# gridstart[1] " << gridstart[1] << endl;
    gridstart[2]=gridcenter[2]-0.5*(gridspacing*GPZ)-gridspacing/2.0;
    os << "# gridstart[2] " << gridstart[2] << endl;
  } else {
    gridstart[0] = gridstartX;
    gridstart[1] = gridstartY;
    gridstart[2] = gridstartZ;
  }

  os << "# gridcenter x,y,z = " << gridcenter[0] << " " << gridcenter[1] << " " << gridcenter[2] << endl;
  os << "# gridstart x,y,z = " << gridstart[0] << " " << gridstart[1] << " " << gridstart[2] << endl;
  os << "# gridspacing " << gridspacing << endl;
  
  // if necessary, shift some atoms next to the end of the box inwards (may sometimes be a problem with ions)

  if (pbc == true) {
    //    os << "Now checking if one or more atoms are close to the border of the box. I might increase the number of gridpoints by 1!" << endl;
    //    increasebox(os);
    os << endl << "# Now doing a GRIDCHECK for all atoms with potential shrinking of radii when atoms are near the border of the box" << endl;
    atomshift(os);
  }
  
  // check whether everything is on the grid
  os << endl << "# Now doing a GRIDCHECK for all atoms with try and catch procedure for the case that anywthing went wrong" << endl;
  gridcheck(os);
  
  
  
  // charge the grid
  chargeGridtrilinear(os);
 
  
  //assign Debye Hueckel boundary charge density at the edge of the box
  //only if NPBC!
  if (!pbc) {
    DebyeHueckel(gridspacing, gridstart, ppp.getKappa(), rhogrid, os);
  }
  
  
  
  
  
  
  
  //  os << "# do radiusboundaries" << endl;
  // make the permittivity maps
  radiusboundaryEPSI(epsIgrid);
  radiusboundaryEPSJ(epsJgrid);
  radiusboundaryEPSK(epsKgrid);
  
  // create the diagonal elements
  //  double epsCGrid[GPXGPYGPZ];
  
  // os << "# do setboundarySolvent" << endl;
  // handle the points on the boundary
  
  setboundarySolvent();
  
  //os << "# ddd" << endl;
  //double epss4PI=epssolvent/(4*ppp.getPI());
  
  
  int count=0;
  
  
  for (int k=0; k <  GPZ; ++k) {
    for (int j=0; j <  GPY; ++j) {
      for (int i=0; i <  GPX; ++i) {
	if (i != 0)     epsCgrid[count] += epsIgrid[count - 1];
	if (i != GPX-1) epsCgrid[count] += epsIgrid[count];
	if (j != 0)     epsCgrid[count] += epsJgrid[count - (GPX)];
	if (j != GPY-1) epsCgrid[count] += epsJgrid[count];
	if (k != 0)     epsCgrid[count] += epsKgrid[count - ((GPX) * (GPY))];
	if (k != GPZ-1) epsCgrid[count] += epsKgrid[count];
	
	count++;
      }
    }
  }  
}


bool FDPoissonBoltzmann::solveforpotential_pbc(int maxits, double acceptance, FDPoissonBoltzmann_ICCG_PBC iccg, ofstream &os){
  // solve the linearized PB equation for the electrostatic potential
  // phigrid using the incomplete cholesky conjugate gradient algorithm
  
  
  // In principle we want to solve this: A x = b
  // where A is the symmentric coefficient matrix, b the constant
  // vector of the source charges (rhogrid) and x the to be calculated
  // solution.
  // The following algorithm as reported and implemented in the UHBD
  // package has the advantage that it does not require to store the
  // full matrix A as such. Obviously, the elements of A are already
  // somehow contained in the vectors/arrays we have set up above,
  // i.e. the permittivity grids. Therefore:
  // epsCgrid --> the diagonal of A
  // epsIgrid --> the 1st subdiagonal of A
  // epsJgrid --> the Imth subdiagonal of A
  // epsKgrid --> the Im * Jmth subdiagonal of A
  
  
  
  bool converged = false;
  
  std::vector<double> zvec;
  std::vector<double> pvec;
  std::vector<double> ldiag;
  
  
  
  double znorm;
  
  
  
  zvec.resize(GPXGPYGPZ, 0.0);
  pvec.resize(GPXGPYGPZ, 0.0);
  ldiag.resize(GPXGPYGPZ, 0.0);
  
  
  
  //this makes the diagonal elemets of the preconditioner matrix
  iccg.initpciccg(ldiag, epsCgrid, epsIgrid, epsJgrid, epsKgrid);
  
  znorm = 0;
  for (int i=0; i < GPXGPYGPZ; ++i) {
    pvec[i] = phigrid[i];
    //   znorm += Math.abs(rhogrid[i]);
    znorm+=fabs(rhogrid[i]);
  }
  
  
  os << "# @ ZNORM "  << znorm << endl;
  
  
  //build the actual preconditioner matrix (zvec)
  
  iccg.gqact(zvec,pvec,epsCgrid,epsIgrid,epsJgrid,epsKgrid);
  
  double aalpha;
  double bbeta = 0;
  double pdotz;
  
  double rdotz1 = 0;
  double rdotz2;
  double anorm = 0;
  
  int iter = 1;
  
  anorm = 0;
  
  //Initially do a first pass. This indeed means some
  //code duplication, but kills two inner loop 'if' statements.
  //The latter was really bad in terms of performance...
  for (int i=0; i < GPXGPYGPZ; ++i) {
    rhogrid[i] -= zvec[i];
    // anorm += Math.abs(rhogrid[i]);
    anorm += fabs(rhogrid[i]);
    
  }
  
  
  //check for convergence
  if (anorm <= acceptance * znorm) {
    // CONVERGENCE REACHED - EXIT ROUTINE
    os << "# CONVERGED AFTER "  << iter << " iterations..." << endl;
    os << "# Exit: solveforpotential()" << endl ;
    converged = true;
    return converged;
  }
  
  else if(iter == maxits) {
    
    os << "# NO CONVERGENCE AFTER " << iter << " iterations..." << endl ;
    os << "# RESULTS MAY NOT BE CORRECT!" << endl;
    os << "# Iterations "  << iter << endl;
    os << "# anorm " << anorm << endl;
    os << "# acceptance * znorm " << (acceptance*znorm) << endl;
    os << "# Exit: solveforpotential()" << endl;
    
    return false;
  }
  
  
  
  
  iccg.pciccg(zvec,ldiag,rhogrid, epsIgrid, epsJgrid, epsKgrid);
  
  rdotz2 = 0;
  for (int i=0; i < GPXGPYGPZ; ++i) {
    rdotz2 += rhogrid[i] * zvec[i];
    
  }
  
  
  
  
  for (int i=0; i < GPXGPYGPZ; ++i) pvec[i] = bbeta * pvec[i] + zvec[i];
  
  rdotz1 = rdotz2;
  
  //C***** GET Z = A * P
  iccg.gqact(zvec, pvec,epsCgrid, epsIgrid,epsJgrid, epsKgrid);
  
  
  //start main iteration loop
  
  while (iter != maxits) {
    ++iter;
    
    pdotz = 0;
    for (int i=0; i < GPXGPYGPZ; ++i) pdotz += pvec[i] * zvec[i];
    aalpha = rdotz1/pdotz;
    
    
    anorm = 0.0;
    
    for (int i=0; i < GPXGPYGPZ; ++i) {
      phigrid[i] += aalpha * pvec[i];
      rhogrid[i] -= aalpha * zvec[i];
      //anorm += Math.abs(rhogrid[i]);
      anorm += fabs(rhogrid[i]);
    }
    
    if (iter%100 == 0)
      os << "# anorm " << anorm << endl;
    
    //check for convergence
    
    if (anorm <= acceptance * znorm) {
      // CONVERGENCE REACHED - EXIT ROUTINE
      os << "# iter anorm acceptance*znorm " << iter << " " << anorm << " " << acceptance*znorm << endl;
      os << "# CONVERGED AFTER " << iter << " iterations..." << endl;
      os << "# Exit: solveforpotential()" << endl;
      converged = true;
      break ; // mainloop;
    }
    
    else if(iter == maxits) {
      
      os << "# NO CONVERGENCE AFTER " << iter << " iterations..." << endl;
      os << "# RESULTS MAY NOT BE CORRECT!" << endl;
      os << "# Iterations "  << iter << endl;
      os << "# anorm " << anorm << endl;
      os << "# acceptance * znorm " << (acceptance*znorm) << endl;
      os << "# Exit: solveforpotential()" << endl;
      break ; // mainloop;
    }
    
    
    
    
    iccg.pciccg(zvec,ldiag,rhogrid, epsIgrid, epsJgrid, epsKgrid);
    
    rdotz2 = 0;
    //omp parallel for reduction(+: double rdotz2)
    for (int i=0; i < GPXGPYGPZ; ++i) rdotz2 += rhogrid[i] * zvec[i];
    
    bbeta = rdotz2/rdotz1;
    
    
    for (int i=0; i < GPXGPYGPZ; ++i) pvec[i] = bbeta * pvec[i] + zvec[i];
    
    rdotz1 = rdotz2;
    
    //C***** GET Z = A * P
    iccg.gqact(zvec, pvec,epsCgrid, epsIgrid,epsJgrid, epsKgrid);
    
  }//while end
  
  
  
  return converged;
  
}





bool FDPoissonBoltzmann::solveforpotential_npbc(int maxits, double acceptance, FDPoissonBoltzmann_ICCG_NPBC iccg, ofstream &os){
  // phigrid using the incomplete cholesky conjugate gradient algorithm
  
  
  // In principle we want to solve this: A x = b
  // where A is the symmentric coefficient matrix, b the constant
  // vector of the source charges (rhogrid) and x the to be calculated
  // solution.
  // The following algorithm as reported and implemented in the UHBD
  // package has the advantage that it does not require to store the
  // full matrix A as such. Obviously, the elements of A are already
  // somehow contained in the vectors/arrays we have set up above,
  // i.e. the permittivity grids. Therefore:
  // epsCgrid --> the diagonal of A
  // epsIgrid --> the 1st subdiagonal of A
  // epsJgrid --> the Imth subdiagonal of A
  // epsKgrid --> the Im * Jmth subdiagonal of A
  
  
  
  bool converged = false;
  
  std::vector<double> zvec;
  std::vector<double> pvec;
  std::vector<double> ldiag;
  
  
  
  double znorm;
  
  
  
  zvec.resize(GPXGPYGPZ, 0.0);
  pvec.resize(GPXGPYGPZ, 0.0);
  ldiag.resize(GPXGPYGPZ, 0.0);
  
  
  
  //this makes the diagonal elemets of the preconditioner matrix
  iccg.initpciccg(ldiag, epsCgrid, epsIgrid, epsJgrid, epsKgrid);
  //for (int i=0;i<GPXGPYGPZ;i++){
  //  os << "# i "  << i << " ldiag[i] " << ldiag[i] << endl;
  //}
  
  znorm = 0;
  for (int i=0; i < GPXGPYGPZ; ++i) {
    pvec[i] = phigrid[i];
    //   znorm += Math.abs(rhogrid[i]);
    znorm+=fabs(rhogrid[i]);
  }
  
  
  os << "# @ ZNORM "  << znorm << endl;
  
  //build the actual preconditioner matrix (zvec)
  
  iccg.gqact(zvec,pvec,epsCgrid,epsIgrid,epsJgrid,epsKgrid);
  
  double aalpha;
  double bbeta = 0;
  double pdotz;
  
  double rdotz1 = 0;
  double rdotz2;
  double anorm = 0;
  
  int iter = 1;
  
  anorm = 0;
  
  //Initially do a first pass. This indeed means some
  //code duplication, but kills two inner loop 'if' statements.
  //The latter was really bad in terms of performance...
  for (int i=0; i < GPXGPYGPZ; ++i) {
    rhogrid[i] -= zvec[i];
    // anorm += Math.abs(rhogrid[i]);
    anorm += fabs(rhogrid[i]);
    
  }
  
  
  //check for convergence
  if (anorm <= acceptance * znorm) {
    // CONVERGENCE REACHED - EXIT ROUTINE
    os << "# CONVERGED AFTER "  << iter << " iterations..." << endl;
    os << "# Exit: solveforpotential()" << endl ;
    converged = true;
    return converged;
  }
  
  else if(iter == maxits) {
    
    os << "# NO CONVERGENCE AFTER " << iter << " iterations..." << endl ;
    os << "# RESULTS MAY NOT BE CORRECT!" << endl;
    os << "# Iterations "  << iter << endl;
    os << "# anorm " << anorm << endl;
    os << "# acceptance * znorm " << (acceptance*znorm) << endl;
    os << "# Exit: solveforpotential()" << endl;
    
    return false;
  }
  
  
  
  
  iccg.pciccg(zvec,ldiag,rhogrid, epsIgrid, epsJgrid, epsKgrid);
  
  rdotz2 = 0;
  for (int i=0; i < GPXGPYGPZ; ++i) {
    rdotz2 += rhogrid[i] * zvec[i];
    
  }
  
  
  
  
  for (int i=0; i < GPXGPYGPZ; ++i) pvec[i] = bbeta * pvec[i] + zvec[i];
  
  rdotz1 = rdotz2;
  
  //C***** GET Z = A * P
  iccg.gqact(zvec, pvec,epsCgrid, epsIgrid,epsJgrid, epsKgrid);
  
  
  //start main iteration loop
  
  while (iter != maxits) {
    ++iter;
    
    pdotz = 0;
    for (int i=0; i < GPXGPYGPZ; ++i) pdotz += pvec[i] * zvec[i];
    aalpha = rdotz1/pdotz;
    
    
    anorm = 0.0;
    
    for (int i=0; i < GPXGPYGPZ; ++i) {
      phigrid[i] += aalpha * pvec[i];
      rhogrid[i] -= aalpha * zvec[i];
      //anorm += Math.abs(rhogrid[i]);
      anorm += fabs(rhogrid[i]);
    }
    
    os << "# iter anorm acceptance*znorm " << iter << " " << anorm << " " << acceptance*znorm << endl;
    
    //check for convergence
    
    if (anorm <= acceptance * znorm) {
      // CONVERGENCE REACHED - EXIT ROUTINE
      os << "# iter anorm acceptance*znorm " << iter << " " << anorm << " " << acceptance*znorm << endl;
      os << "# CONVERGED AFTER " << iter << " iterations..." << endl;
      os << "# Exit: solveforpotential()" << endl;
      converged = true;
      break ; // mainloop;
    }
    
    else if(iter == maxits) {
      
      os << "# NO CONVERGENCE AFTER " << iter << " iterations..." << endl;
      os << "# RESULTS MAY NOT BE CORRECT!" << endl;
      os << "# Iterations "  << iter << endl;
      os << "# anorm " << anorm << endl;
      os << "# acceptance * znorm " << (acceptance*znorm) << endl;
      os << "# Exit: solveforpotential()" << endl;
      break ; // mainloop;
    }
    
    
    
    
    iccg.pciccg(zvec,ldiag,rhogrid, epsIgrid, epsJgrid, epsKgrid);
    
    rdotz2 = 0;
    //omp parallel for reduction(+: double rdotz2)
    for (int i=0; i < GPXGPYGPZ; ++i) rdotz2 += rhogrid[i] * zvec[i];
    
    bbeta = rdotz2/rdotz1;
    
    
    for (int i=0; i < GPXGPYGPZ; ++i) pvec[i] = bbeta * pvec[i] + zvec[i];
    
    rdotz1 = rdotz2;
    
    //C***** GET Z = A * P
    iccg.gqact(zvec, pvec,epsCgrid, epsIgrid,epsJgrid, epsKgrid);
    
  }//while end
  
  
  
  return converged;
  
}




double FDPoissonBoltzmann::dGelec(ofstream &os, vector<double> *potentials){
  
  
  //double potential = 0;
  double potential_rest = 0;
  
  //int I, J, K;
  
  
  double phitot = 0.0;
  
  
  for (int i=0; i < GPXGPYGPZ;++i) {
    phitot += phigrid[i];
  }
  
  
  os << "# Average potential [kJ/mol] " <<  phitot/GPXGPYGPZ*ppp.getFPEPSI() << endl;
  
  //in the periodic case, subtract the average potential...
  if ( pbc ) {
    
    //multiply by FPEPSI
    for (int i=0; i < GPXGPYGPZ;++i) phigrid[i] *= ppp.getFPEPSI();
    
    double phiaver = 0, phiaver2 = 0;
    
    
    for (int i=0; i < GPXGPYGPZ;++i) phiaver += phigrid[i];
    phiaver /= GPXGPYGPZ;

    for (int i=0; i < GPXGPYGPZ;++i) {
      phigrid[i] -= phiaver;
      phiaver2 += phigrid[i];
    }

    //divide by FPEPSI
    for (int i=0; i < GPXGPYGPZ;++i)
      phigrid[i] /=ppp.getFPEPSI();

    os << "# AVERAGE POTENTIAL SUBTRACTED [kJ/mol]: " << phiaver << " -> " <<  phiaver2/GPXGPYGPZ << endl;

  }

  //                potential = getdG();
  //                potential = potential * ppp.getFPEPSI();
  //                os << "# DG RETURNED [kJ/mol]: " <<  potential << endl;
  potential_rest = getdG_restricted(os, potentials);
  potential_rest = potential_rest * ppp.getFPEPSI();
  os << "# DG RETURNED [kJ/mol]: " <<  potential_rest << endl;
  //  return potential;
  return potential_rest;



}





double FDPoissonBoltzmann::getdG(){

  int size;
  size=atoms.size();

  double ijk[3] ;
  double fraction[3] ;
  double pottmp = 0.0, potential = 0.0;

  int I = 0, J = 0, K = 0;

  for (int i=0; i < size; ++i) { // loop over atoms




    //convert to grid units
    for (int j=0; j < 3; ++j) { // loop over x, y, z coors of that atom
      ijk[j] = ((atoms.pos(i))[j] - gridstart[j])/gridspacing;
      //get the fractional remainder
      fraction[j] = ijk[j] - (int) ijk[j];
    }

    //subtract 1, because we are counting from 0
    I = (int) ijk[0] - 1;
    J = (int) ijk[1] - 1;
    K = (int) ijk[2] - 1;



    pottmp =        (1-fraction[0]) * (1-fraction[1]) * (1-fraction[2]) * phigrid[index(I,J,K)]
      + (fraction[0]) * (1-fraction[1]) * (1-fraction[2]) * phigrid[index(I+1,J,K)]
      + (1-fraction[0]) * (fraction[1]) * (1-fraction[2]) * phigrid[index(I,J+1,K)]
      + (fraction[0]) * (fraction[1]) * (1-fraction[2]) * phigrid[index(I+1,J+1,K)]
      + (1-fraction[0]) * (1-fraction[1]) * (fraction[2]) * phigrid[index(I,J,K+1)]
      + (fraction[0]) * (1-fraction[1]) * (fraction[2]) * phigrid[index(I+1,J,K+1)]
      + (1-fraction[0]) * (fraction[1]) * (fraction[2]) * phigrid[index(I,J+1,K+1)]
      + (fraction[0]) * (fraction[1]) * (fraction[2]) * phigrid[index(I+1,J+1,K+1)];


    potential += 0.5 * pottmp * atoms.charge(i);

                      

  } // end of loop over atoms

  return potential;
}




double FDPoissonBoltzmann::getdG_restricted(ofstream &os, vector<double> *potentials){

  int size;
  size=atoms_to_charge.size();

  double ijk[3] ;
  double fraction[3] ;
  double pottmp = 0.0, potential = 0.0;

  int I = 0, J = 0, K = 0;

  for (int i=0; i < size; ++i) { // loop over atoms

              


    //convert to grid units
    for (int j=0; j < 3; ++j) { // loop over x, y, z coors of that atom
      ijk[j] = ((atoms_to_charge.pos(i))[j] - gridstart[j])/gridspacing;
      //get the fractional remainder
      fraction[j] = ijk[j] - (int) ijk[j];
    }

    //subtract 1, because we are counting from 0
    I = (int) ijk[0] - 1;
    J = (int) ijk[1] - 1;
    K = (int) ijk[2] - 1;



    pottmp =        (1-fraction[0]) * (1-fraction[1]) * (1-fraction[2]) * phigrid[index(I,J,K)]
      + (fraction[0]) * (1-fraction[1]) * (1-fraction[2]) * phigrid[index(I+1,J,K)]
      + (1-fraction[0]) * (fraction[1]) * (1-fraction[2]) * phigrid[index(I,J+1,K)]
      + (fraction[0]) * (fraction[1]) * (1-fraction[2]) * phigrid[index(I+1,J+1,K)]
      + (1-fraction[0]) * (1-fraction[1]) * (fraction[2]) * phigrid[index(I,J,K+1)]
      + (fraction[0]) * (1-fraction[1]) * (fraction[2]) * phigrid[index(I+1,J,K+1)]
      + (1-fraction[0]) * (fraction[1]) * (fraction[2]) * phigrid[index(I,J+1,K+1)]
      + (fraction[0]) * (fraction[1]) * (fraction[2]) * phigrid[index(I+1,J+1,K+1)];

    // print out the electrostatic potential at that atom site
    double pottmpfpep = pottmp  * ppp.getFPEPSI();
    if (potentials != NULL ) { // only fill with potentials if vector was passed to the function
      potentials->push_back(pottmpfpep);
    }
    os << "# atom  " << i+1 << " charge " << atoms_to_charge.charge(i)  << " ele. pot. " << pottmpfpep << endl;

    potential += 0.5 * pottmp * atoms_to_charge.charge(i);



  } // end of loop over atoms

  return potential;
}

void FDPoissonBoltzmann::increasebox(ofstream &os){
  int size = atoms.size();
  double gridposition_X_lower = gridstart[0];
  double gridposition_X_upper = gridstart[0] + GPX * gridspacing;
  double gridposition_Y_lower = gridstart[1];
  double gridposition_Y_upper = gridstart[1] + GPZ * gridspacing;
  double gridposition_Z_lower = gridstart[2];
  double gridposition_Z_upper = gridstart[2] + GPZ * gridspacing;
  double protruding_X = 0;
  double protruding_Y = 0;
  double protruding_Z = 0;
  double temp_lower = 0;
  double temp_upper = 0;

  for (int i=0; i < size; ++i) {
    // X AXIS
    temp_lower = gridposition_X_lower - atoms.pos(i)[0];
    temp_upper = atoms.pos(i)[0] - gridposition_X_upper;
    if ( temp_lower > protruding_X ) {
      protruding_X = temp_lower;
    }
    if ( temp_upper > protruding_X ) {
      protruding_X = temp_upper;
    }

    // Y AXIS
    temp_lower = gridposition_Y_lower - atoms.pos(i)[1];
    temp_upper = atoms.pos(i)[1] - gridposition_Y_upper;
    if ( temp_lower > protruding_Y ) {
      protruding_Y = temp_lower;
    }
    if ( temp_upper > protruding_Y ) {
      protruding_Y = temp_upper;
    }

        // Z AXIS
    temp_lower = gridposition_Z_lower - atoms.pos(i)[2];
    temp_upper = atoms.pos(i)[2] - gridposition_Z_upper;
    if ( temp_lower > protruding_Z ) {
      protruding_Z = temp_lower;
    }
    if ( temp_upper > protruding_Z ) {
      protruding_Z = temp_upper;
    }
  }
  if (protruding_X < gridspacing) {
    GPX += 1;
    os << "Increased the number of gridpoints in the X direction by 1 because one or more atoms were too close to the border" << endl;
  } else {
    throw gromos::Exception("FDPoissonBolzmann","Atom extending grid. Exiting ...");
  }
  if (protruding_Y < gridspacing) {
    GPY += 1;
    os << "Increased the number of gridpoints in the Y direction by 1 because one or more atoms were too close to the border" << endl;
  } else {
    throw gromos::Exception("FDPoissonBolzmann","Atom extending grid. Exiting ...");
  }
    if (protruding_Z < gridspacing) {
    GPX += 1;
    os << "Increased the number of gridpoints in the Z direction by 1 because one or more atoms were too close to the border" << endl;
  } else {
    throw gromos::Exception("FDPoissonBolzmann","Atom extending grid. Exiting ...");
  }
}
  
    

void FDPoissonBoltzmann::atomshift(ofstream &os){
  int size = atoms.size();

  for (int i=0; i < size; ++i) {
    os << "# GRIDCHECK for atom " << i+1 << " radius " << fixed << std::setprecision(7) << radii[i]
       << " charge " << fixed << std::setprecision(3) << setw(6) << atoms.charge(i)
       << " pos " << atoms.pos(i)[0] << " " << atoms.pos(i)[1] << " " << atoms.pos(i)[2] << endl;
    
      bool gridfail = 0;
      char coordfail = '0';
      double atomposition = 0;
      double gridposition = 0;
      double protruding = 0;
      
      //X-axis upper end atom exceeding
      gridposition = gridstart[0] + GPX * gridspacing;
      if ( atoms.pos(i)[0] > gridposition ) {
	coordfail='X';
	protruding = atoms.pos(i)[0] - gridposition;
	atoms.pos(i)[0] = gridstart[0]+protruding;
	os << "# Atom extending grid. Atom " << i+1 << " extends the " << coordfail << " axis by " << protruding << " nm" << endl;
	os << "# You may consider the @increasegrid option of dGslv_pbsolv" << endl;
	throw gromos::Exception("FDPoissonBolzmann","Atom extending grid. See the log data (@verbose flag) and search for the line starting with <# Atom extending grid>. Exiting ...");
	// os << "# Needed to shift atom " << i+1 << " on the " << coordfail << " axis to the other side of the box " << endl;
	// os << "# New coordinates atom " << i+1  << " " << atoms.pos(i)[0] << " " << atoms.pos(i)[1] << " " << atoms.pos(i)[2] << endl;
      };

      //X-axis upper end vdw sphere exceeding
      gridposition = gridstart[0] + GPX * gridspacing;
      if ( atoms.pos(i)[0] + radii[i] > gridposition) {
	coordfail='X';
	protruding = atoms.pos(i)[0] + radii[i] - gridposition + gridspacing/100;
	radii[i] = radii[i] - protruding;
	os << "# Needed to shrink the radius of atom " << i+1 << " on the " << coordfail << " axis since it exceeded the grid " << endl;
	os << "# New radius atom " << i+1  << " " << radii[i] << endl;
      };

      //X-axis lower end atom exceeding
      gridposition = gridstart[0];
      if ( atoms.pos(i)[0] < gridposition ) {
	coordfail='X';
	protruding = atoms.pos(i)[0] - gridposition;
	atoms.pos(i)[0] = gridstart[0] + GPX * gridspacing + protruding;
	os << "# Atom extending grid. Atom " << i+1 << " extends the " << coordfail << " axis by " << protruding << " nm" << endl;
	os << "# You may consider the @increasegrid option of dGslv_pbsolv" << endl;
	throw gromos::Exception("FDPoissonBolzmann","Atom extending grid. See the log data (@verbose flag) and search for the line starting with <# Atom extending grid>. Exiting ...");
	// os << "# Needed to shift atom " << i+1 << " on the " << coordfail << " axis to the other side of the box " << endl;
	// os << "# New coordinates atom " << i+1  << " " << atoms.pos(i)[0] << " " << atoms.pos(i)[1] << " " << atoms.pos(i)[2] << endl;
      };

      //X-axis lower end vdw sphere exceeding
      gridposition = gridstart[0];
      if ( atoms.pos(i)[0] - radii[i] < gridposition ) {
	coordfail='X';
	protruding = atoms.pos(i)[0] - radii[i] - gridposition - gridspacing/100;
	radii[i] = radii[i] + protruding;
	os << "# Needed to shrink the radius of atom " << i+1 << " on the " << coordfail << " axis since it exceeded the grid " << endl;
	os << "# New radius atom " << i+1  << " " << radii[i] << endl;
      };

      //X-axis upper end vdw sphere exceeding (another one, for the case it was shifted from the lower to the upper end, but radius still exceeds)
      gridposition = gridstart[0] + GPX * gridspacing;
      if ( atoms.pos(i)[0] + radii[i] > gridposition) {
	coordfail='X';
	protruding = atoms.pos(i)[0] + radii[i] - gridposition + gridspacing/100;
	radii[i] = radii[i] - protruding;
	os << "# Needed to shrink the radius of atom " << i+1 << " on the " << coordfail << " axis since it exceeded the grid " << endl;
	os << "# New radius atom " << i+1  << " " << radii[i] << endl;
      };




      //Y-axis upper end atom exceeding
      gridposition = gridstart[1] + GPY * gridspacing;
      if ( atoms.pos(i)[1] > gridposition) {
	coordfail='Y';
	protruding = atoms.pos(i)[1] - gridposition;
	atoms.pos(i)[1] = gridstart[1]+protruding;
	os << "# Atom extending grid. Atom " << i+1 << " extends the " << coordfail << " axis by " << protruding << " nm" << endl;
	os << "# You may consider the @increasegrid option of dGslv_pbsolv" << endl;
	throw gromos::Exception("FDPoissonBolzmann","Atom extending grid. See the log data (@verbose flag) and search for the line starting with <# Atom extending grid>. Exiting ...");
	// os << "# Needed to shift atom " << i+1 << " on the " << coordfail << " axis to the other side of the box " << endl;
	// os << "# New coordinates atom " << i+1  << " " << atoms.pos(i)[0] << " " << atoms.pos(i)[1] << " " << atoms.pos(i)[2] << endl;
      };

      //Y-axis upper end vdw sphere exceeding
      gridposition = gridstart[1] + GPY * gridspacing;
      if ( atoms.pos(i)[1] + radii[i] > gridposition) {
	coordfail='Y';
	protruding = atoms.pos(i)[1] + radii[i] - gridposition + gridspacing/100;
	radii[i] = radii[i] - protruding;
	os << "# Needed to shrink the radius of atom " << i+1 << " on the " << coordfail << " axis since it exceeded the grid " << endl;
	os << "# New radius atom " << i+1  << " " << radii[i] << endl;
      };

      //Y-axis lower end atom exceeding
      gridposition = gridstart[1];
      if ( atoms.pos(i)[1] < gridposition ) {
	coordfail='Y';
	protruding = atoms.pos(i)[1] - gridposition;
	atoms.pos(i)[1] = gridstart[1] + GPY * gridspacing + protruding;
	os << "# Atom extending grid. Atom" << i+1 << " extends the " << coordfail << " axis by " << protruding << " nm" << endl;
	os << "# You may consider the @increasegrid option of dGslv_pbsolv" << endl;
	throw gromos::Exception("FDPoissonBolzmann","Atom extending grid. See the log data (@verbose flag) and search for the line starting with <# Atom extending grid>. Exiting ...");
	// os << "# Needed to shift atom " << i+1 << " on the " << coordfail << " axis to the other side of the box " << endl;
	// os << "# New coordinates atom " << i+1  << " " << atoms.pos(i)[0] << " " << atoms.pos(i)[1] << " " << atoms.pos(i)[2] << endl;
      };

      //Y-axis lower end vdw sphere exceeding
      gridposition = gridstart[1];
      if ( atoms.pos(i)[1] - radii[i] < gridposition ) {
	coordfail='Y';
	protruding = atoms.pos(i)[1] - radii[i] - gridposition - gridspacing/100;
	radii[i] = radii[i] + protruding;
	os << "# Needed to shrink the radius of atom " << i+1 << " on the " << coordfail << " axis since it exceeded the grid " << endl;
	os << "# New radius atom " << i+1  << " " << radii[i] << endl;
      };

      //Y-axis upper end vdw sphere exceeding (another one, for the case it was shifted from the lower to the upper end, but radius still exceeds)
      gridposition = gridstart[1] + GPY * gridspacing;
      if ( atoms.pos(i)[1] + radii[i] > gridposition ) {
	coordfail='Y';
	protruding = atoms.pos(i)[1] + radii[i] - gridposition + gridspacing/100;
	radii[i] = radii[i] - protruding;
	os << "# Needed to shrink the radius of atom " << i+1 << " on the " << coordfail << " axis since it exceeded the grid " << endl;
	os << "# New radius atom " << i+1  << " " << radii[i] << endl;
      };
      

      //Z-axis upper end atom exceeding
      gridposition = gridstart[2] + GPZ * gridspacing;
      if ( atoms.pos(i)[2] > gridposition ) {
	coordfail='Z';
	protruding = atoms.pos(i)[2] - gridposition;
	atoms.pos(i)[2] = gridstart[2]+protruding;
	os << "# Atom extending grid. Atom " << i+1 << " extends the " << coordfail << " axis by " << protruding << " nm" << endl;
	os << "# You may consider the @increasegrid option of dGslv_pbsolv" << endl;
	throw gromos::Exception("FDPoissonBolzmann","Atom extending grid. See the log data (@verbose flag) and search for the line starting with <# Atom extending grid>. Exiting ...");
	// os << "# Needed to shift atom " << i+1 << " on the " << coordfail << " axis to the other side of the box " << endl;
	// os << "# New coordinates atom " << i+1  << " " << atoms.pos(i)[0] << " " << atoms.pos(i)[1] << " " << atoms.pos(i)[2] << endl;
      };

      //Z-axis upper end vdw sphere exceeding
      gridposition = gridstart[2] + GPZ * gridspacing;
      if ( atoms.pos(i)[2] + radii[i] > gridposition) {
	coordfail='Z';
	protruding = atoms.pos(i)[2] + radii[i] - gridposition + gridspacing/100;
	radii[i] = radii[i] - protruding;
	os << "# Needed to shrink the radius of atom " << i+1 << " on the " << coordfail << " axis since it exceeded the grid " << endl;
	os << "# New radius atom " << i+1  << " " << radii[i] << endl;
      };

      //Z-axis lower end atom exceeding
      gridposition = gridstart[2];
      if ( atoms.pos(i)[2] < gridposition) {
	coordfail='Z';
	protruding = atoms.pos(i)[2] - gridposition;
	atoms.pos(i)[2] = gridstart[2] + GPZ * gridspacing + protruding;
	os << "# Atom extending grid. Atom " << i+1 << " extends the " << coordfail << " axis by " << protruding << " nm" << endl;
	os << "# You may consider the @increasegrid option of dGslv_pbsolv" << endl;
	throw gromos::Exception("FDPoissonBolzmann","Atom extending grid. See the log data (@verbose flag) and search for the line starting with <# Atom extending grid>. Exiting ...");
	// os << "# Needed to shift atom " << i+1 << " on the " << coordfail << " axis to the other side of the box " << endl;
	// os << "# New coordinates atom " << i+1  << " " << atoms.pos(i)[0] << " " << atoms.pos(i)[1] << " " << atoms.pos(i)[2] << endl;
      };

      //Z-axis lower end vdw sphere exceeding
      gridposition = gridstart[2];
      if ( atoms.pos(i)[2] - radii[i] < gridposition ) {
	coordfail='Z';
	protruding = atoms.pos(i)[2] - radii[i] - gridposition - gridspacing/100;
	radii[i] = radii[i] + protruding;
	os << "# Needed to shrink the radius of atom " << i+1 << " on the " << coordfail << " axis since it exceeded the grid " << endl;
	os << "# New radius atom " << i+1  << " " << radii[i] << endl;
      };

      //Z-axis upper end vdw sphere exceeding (another one, for the case it was shifted from the lower to the upper end, but radius still exceeds)
      gridposition = gridstart[2] + GPZ * gridspacing;
      if ( atoms.pos(i)[2] + radii[i] > gridposition) {
	coordfail='Z';
	protruding = atoms.pos(i)[2] + radii[i] - gridposition + gridspacing/100;
	radii[i] = radii[i] - protruding;
	os << "# Needed to shrink the radius of atom " << i+1 << " on the " << coordfail << " axis since it exceeded the grid " << endl;
	os << "# New radius atom " << i+1  << " " << radii[i] << endl;
      };
      

      
      
  }// end of for loop
} // end of atomshift



void FDPoissonBoltzmann::increasegrid(ofstream &os){
  GPX += std::ceil(4/gridspacing);
  GPY += std::ceil(4/gridspacing);
  GPZ += std::ceil(4/gridspacing);
  os << "For NPBC, I increased the number of gridpoints to X = " << GPX << " Y = " << GPY << " Z = " << GPZ  << endl; 
}
  


void FDPoissonBoltzmann::gridcheck(ofstream &os){
  int size = atoms.size();
  double rad;

  for (int i=0; i < size; ++i) {
    rad = radii[i];
    
    os << "# GRIDCHECK for atom " << i+1 << " radius " << rad << " charge " << atoms.charge(i) <<
      " pos " << atoms.pos(i)[0] << " " << atoms.pos(i)[1] << " " << atoms.pos(i)[2] << endl;
    
    try{
      if (           ((atoms.pos(i))[0] + rad) > (gridstart[0] + GPX * gridspacing)
		     || ((atoms.pos(i))[0] - rad) < (gridstart[0])
		     || ((atoms.pos(i))[1]+ rad) > (gridstart[1] + GPY * gridspacing)
		     || ((atoms.pos(i))[1] - rad) < (gridstart[1])
		     || ((atoms.pos(i))[2]+ rad) > (gridstart[2] + GPZ * gridspacing)
		     || ((atoms.pos(i))[2]- rad) < (gridstart[2]) ) {
	
	
	
        
	os << "# ATOM EXTENDING GRID!" << endl;
	os << "# ATOM ID:     " << (i+1)  << endl;
	os << "# ATOM POS:    " << (atoms.pos(i))[0] << " " << (atoms.pos(i))[1] << " " << (atoms.pos(i))[2]<< endl;
	os << "# ATOM RADIUS: "  <<  rad << endl;
	os << "# GRID START   "  << gridstart[0] << " " << gridstart[1] << " " << gridstart[2] << endl;
	os << "# GRID END     "  <<  (gridstart[0] + GPX * gridspacing) << " " << (gridstart[1] + GPY * gridspacing) << " " << (gridstart[2] + GPZ * gridspacing) << endl;
        
        
	throw gromos::Exception("FDPoissonBolzmann","Atom extending grid. Exiting ...");
        
      } // if endif
    } // end of try
    
    catch (const gromos::Exception &e){
      cerr << e.what() << endl;
      exit(1);
    }
  } // loop over atoms 
} // end of gridcheck

int FDPoissonBoltzmann:: index(int x, int y, int z) {
  //convert 3D to 1D array index

  if (pbc) {
    return (x + (GPX) * (y + (GPY) * (z%GPZ)));
  } else{
    return (x + (GPX) * (y + (GPY) * (z)));
  }
}

        
// void  FDPoissonBoltzmann::setpermittivity(double solute, double solvent) {
//        epssolute  = solute;
//        epssolvent = solvent;
// }


void FDPoissonBoltzmann::getgridcenter(double& X, double& Y, double& Z) {
  X=gridcenter[0];
  Y=gridcenter[1];
  Z=gridcenter[2];
}
void FDPoissonBoltzmann::getgridstart(double& X, double& Y, double& Z) {
  X=gridstart[0];
  Y=gridstart[1];
  Z=gridstart[2];
}
//double *FDPoissonBoltzmann::getgridcenter() { //double& X, double& Y, double& Z) {
//  return gridcenter;
//}
//   X=gridcenter[0];
//    Y=gridcenter[1];
//    Z=gridcenter[2];
//  }
//double *FDPoissonBoltzmann::getgridstart() { //double& X, double& Y, double& Z) {
//  return gridstart;
//}
//   X=gridstart[0];
//   Y=gridstart[1];
//   Z=gridstart[2];
//  }


void FDPoissonBoltzmann::setboundarySolvent() {
  // SETS THE BOUNDARY TO SOLVENT.


  double epss4PI = epssolvent/(4 * gmath::physConst.get_pi());

  for (int j=0; j < GPY; ++j) {
    for (int i=0; i < GPX; ++i) {


      epsCgrid[(index(i,j,0))] += epss4PI;
      epsCgrid[(index(i,j,GPZ-1))] += epss4PI;
    }
  }

  for (int k=0; k < GPZ; ++k) {
    for (int i=0; i < GPX; ++i) {
      epsCgrid[(index(i,0,k))] += epss4PI;
      epsCgrid[(index(i,GPY-1,k))] += epss4PI;
    }
  }

  for (int k=0; k < GPZ; ++k) {
    for (int j=0; j < GPY; ++j) {
      epsCgrid[(index(0,j,k))] += epss4PI;
      epsCgrid[(index(GPX-1,j,k))] += epss4PI;
    }
  }

}


void FDPoissonBoltzmann::radiusboundaryEPSI( std::vector<double>    &  epsgrid) {

  int size=atoms.size();

  double epsi4PI = epssolute/(4 * ppp.getPI());
  double epss4PI = epssolvent/(4 * ppp.getPI());


  //set grid values to that of the solvent permittivity

  for (int i=0; i < GPXGPYGPZ; ++i) epsgrid[i] = epss4PI;


  //set extending points to 0 (NPBC) or epssolv (PBC)

  if (!pbc) {
    for (int k = 0; k < GPZ; ++k) {
      for (int j = 0; j < GPY; ++j) {

	epsgrid[index(GPX-1,j,k)] = 0.0;

      }
    }
  }
  else {
    for (int k = 0; k < GPZ; ++k) {
      for (int j = 0; j < GPY; ++j) {
	epsgrid[index(GPX-1,j,k)] = epss4PI;
      }
    }
  }



  //do the grid...

  //first increment by solvent radius

  for (int i=0; i < size; ++i) {

    double range1 = (radii[i])/gridspacing;

    double ijk[3];
    //double fraction[3];

    double fuzoff = 0.5;
    int I, J, K;
    //convert to grid units
    for (int j=0; j < 3; ++j) {
      ijk[j] = (  (atoms.pos(i))[j] - gridstart[j])/gridspacing;
      //get the fractional remainder
      //fraction[j] = ijk[j] - (int) ijk[j];
    }

    int low1 = (int) (ijk[2] - range1 + 1.0);
    int high1 = (int) (ijk[2] + range1);


    for (K = max(1,low1); K <= min(GPZ,high1); ++K) {
      double arg = range1 * range1 - pow((ijk[2] - K),2);
      double range2 = sqrt(arg);
      if (arg < 0) range2 = -range2;
      int low2 = (int) (ijk[1] - range2 + 1);
      int high2 = (int) (ijk[1] + range2);


      for (J = max(1,low2); J <= min(GPY,high2); ++J) {
	arg = range2 * range2 - pow((ijk[1] - J),2);
	double range3 = sqrt(arg);
	if (arg < 0) range2 = -range2;

	if (range3 > 0.0) {
	  double low3 = ijk[0] - range3 + fuzoff;
	  double high3 = ijk[0] + range3 - fuzoff;



	  for ( I = max(1, (int) low3); I <= min(GPX-1, (int) high3); ++I)  {

                         

	    epsgrid[(index(I-1,J-1,K-1))] = epsi4PI;
	  }

	}
      }
    }
  }


}


void FDPoissonBoltzmann::radiusboundaryEPSJ( std::vector<double> &epsgrid) {
  int size=atoms.size();




  double epsi4PI = epssolute/(4 * ppp.getPI());
  double epss4PI = epssolvent/(4 * ppp.getPI());





  //set grid values to that of the solvent permittivity

  for (int i=0; i < GPXGPYGPZ; ++i) epsgrid[i] = epss4PI;

  //set the extending points to 0.0
  if (!pbc) {
    for (int k = 0; k < GPZ; ++k) {
      for (int i = 0; i < GPX; ++i) {
	epsgrid[index(i,GPY-1,k)] = 0.0;
      }
    }
  }
  else {
    for (int k = 0; k < GPZ; ++k) {
      for (int i = 0; i < GPX; ++i) {
	epsgrid[index(i,GPY-1,k)] = epss4PI;
      }
    }
  }
  //do the grid...

  //first increment by solvent radius

  for (int i=0; i < size; ++i) {


    double range1 = (radii[i])/gridspacing;

    double ijk[3];
    //double fraction[3];

    double fuzoff = 0.5;

    int I, J, K;

    //convert to grid units
    for (int j=0; j < 3; ++j) {
      ijk[j] = ((atoms.pos(i))[j] - gridstart[j])/gridspacing;
      //get the fractional remainder
      // fraction[j] = ijk[j] - (int) ijk[j];
    }

    int low1 = (int) (ijk[2] - range1 + 1.0);
    int high1 = (int) (ijk[2] + range1);



    for (K = max(1,low1); K <= min(GPZ,high1); ++K) {

      double arg = range1 * range1 - pow((ijk[2] - K),2);
      double range2 = sqrt(arg);
      if (arg < 0) range2 = -range2;
      int low2 = (int) (ijk[0] - range2 + 1);
      int high2 = (int) (ijk[0] + range2);

      for (I = max(1,low2); I <= min(GPX,high2); ++I) {
	arg = range2 * range2 - pow((ijk[0] - I),2);
	double range3 = sqrt(arg);
	if (arg < 0) range2 = -range2;

	if (range3 > 0.0) {
	  double low3 = ijk[1] - range3 + fuzoff;
	  double high3 = ijk[1] + range3 - fuzoff;

	  for ( J = max(1, (int) low3); J <= min(GPY-1, (int) high3); ++J) epsgrid[(index(I-1,J-1,K-1))] = epsi4PI;

	}
      }
    }
  }


}


void  FDPoissonBoltzmann::radiusboundaryEPSK( std::vector<double> & epsgrid) {

  double epsi4PI = epssolute/(4 * ppp.getPI());
  double epss4PI = epssolvent/(4 * ppp.getPI());

  int size=atoms.size();

  //set grid values to that of the solvent permittivity

  for (int i=0; i <  GPXGPYGPZ; ++i) epsgrid[i] = epss4PI;

  if (!pbc) {
    for (int j = 0; j < GPY; ++j) {
      for (int i = 0; i < GPX; ++i) {
	epsgrid[index(i,j,GPZ-1)] = 0.0;
      }
    }
  }
  else {
    for (int j = 0; j < GPY; ++j) {
      for (int i = 0; i < GPX; ++i) {
	epsgrid[index(i,j,GPZ-1)] = epss4PI;
      }
    }

  }


  //do the grid...

  //first increment by solvent radius

  for (int i=0; i < size; ++i) {


    double range1 = radii[i] /gridspacing;

    double ijk[3];
    //double fraction[3];

    double fuzoff = 0.5;

    int I, J, K;
    //convert to grid units
    for (int j=0; j < 3; ++j) {
      ijk[j] = (   (atoms.pos(i))[j] - gridstart[j])/gridspacing;
      //get the fractional remainder
      // fraction[j] = ijk[j] - (int) ijk[j];
    }

    int low1 = (int) (ijk[1] - range1 + 1.0);
    int high1 = (int) (ijk[1] + range1);



    for (J = max(1,low1); J <= min(GPY,high1); ++J) {

      double arg = range1 * range1 - pow((ijk[1] - J),2);
      double range2 = sqrt(arg);
      if (arg < 0) range2 = -range2;
      int low2 = (int) (ijk[0] - range2 + 1);
      int high2 = (int) (ijk[0] + range2);

      for (I = max(1,low2); I <= min(GPX,high2); ++I) {
	arg = range2 * range2 - pow((ijk[0] - I),2);
	double range3 = sqrt(arg);
	if (arg < 0) range2 = -range2;

	if (range3 > 0.0) {
	  double low3 = ijk[2] - range3 + fuzoff;
	  double high3 = ijk[2] + range3 - fuzoff;




	  for (K = max(1, (int) low3); K <= min(GPZ-1, (int) high3); ++K) epsgrid[(index(I-1,J-1,K-1))] = epsi4PI;



	}
      }
    }
  }


}



void FDPoissonBoltzmann::chargeGridtrilinear(ofstream &os) {
  //spread charges via trilinear interpolation over 8 closest grid points, i.e.
  //put partial charges at points:
  //(i,j,k),
  //(i+1,j,k),
  //(i+1,j+1,k),
  //(i,j+1,k),
  //(i,j+1,k+1),
  //(i+1,j,k+1),
  //(i,j,k+1),
  //(i+1,j+1,k+1)

  //the charge at a gridpoint will then be calculated as:
  // q(gridpoint) += q0 * (1-a) * (1-b) * (1-c)
  // where a,b and c are the fractional distances along the gridaxes from the actual atomic charge



  int size = 0;

  size = atoms.size();


  double charge;


  int I, J, K;

  double ijk[3] ;
  double fraction[3];

  //loop over all atoms
  for (int i=0; i < size; ++i) {
    //convert to grid units
    for (int j=0; j < 3; ++j) {
      ijk[j] = (   (atoms.pos(i))[j] - gridstart[j])/gridspacing;
      //get the fractional remainder
      fraction[j] = ijk[j] - (int) ijk[j];
    }

    charge = atoms.charge(i);

    //substract 1, because we are counting from 0
    I = (int) ijk[0] - 1;
    J = (int) ijk[1] - 1;
    K = (int) ijk[2] - 1;


    //assign to rhogrid
    rhogrid[(index(I,J,K))]       += charge * (1-fraction[0]) * (1-fraction[1]) * (1-fraction[2])/gridspacing;
    rhogrid[(index(I+1,J,K))]     += charge * (fraction[0]) * (1-fraction[1]) * (1-fraction[2])/gridspacing;
    rhogrid[(index(I,J+1,K))]     += charge * (1-fraction[0]) * (fraction[1]) * (1-fraction[2])/gridspacing;
    rhogrid[(index(I+1,J+1,K))]   += charge * (fraction[0]) * (fraction[1]) * (1-fraction[2])/gridspacing;
    rhogrid[(index(I,J,K+1))]     += charge * (1-fraction[0]) * (1-fraction[1]) * (fraction[2])/gridspacing;
    rhogrid[(index(I+1,J,K+1))]   += charge * (fraction[0]) * (1-fraction[1]) * (fraction[2])/gridspacing;
    rhogrid[(index(I,J+1,K+1))]   += charge * (1-fraction[0]) * (fraction[1]) * (fraction[2])/gridspacing;
    rhogrid[(index(I+1,J+1,K+1))] += charge * (fraction[0]) * (fraction[1]) * (fraction[2])/gridspacing;




  } //end loop over all atoms

  // sum all charges;
  // this and the following code is used for adding
  // the homogenous neutralizing background charge density under PBC





  //subtract average charge if under PBC
  if ( pbc ) {
    os << "# @@@@@@@@@@@@@@@ pbc-> subtract aver " << endl;
    double	sumchrg = 0;
    // calculate total charge on grid
    for (int i = 0; i < GPXGPYGPZ; ++i) sumchrg += rhogrid[i] * gridspacing;
    sumchrg /= (GPXGPYGPZ*gridspacing);
    for (int i = 0; i < GPXGPYGPZ; ++i) rhogrid[i] -= sumchrg;//(gpxgpygpz*gridspacing);
  }




  // Maria check for debugging
  double znorm = 0;
  for (int i=0; i < GPXGPYGPZ; ++i) {
    znorm+=fabs(rhogrid[i]);
  }
  os << "# @chargegridtrilinear: znorm : " << znorm << endl;


} //end chargeGridtrilinear





/*double[]  FDPoissonBoltzmann::getpotential() {
  return phigrid;
  }
  double[]  FDPoissonBoltzmann::getgridstart() {
  return gridstart;
  }*/







void FDPoissonBoltzmann::DebyeHueckel(double gridspacing, double gridstart[3], 
				      double kappa,  std::vector<double> & rhogrid, ofstream &os) {
  //this routine adds a Debye-Hueckel potential to the boundary of the grid, i.e.
  //the boundary potential will be a sum of the contributions from the individual atom
  //charges as: PHI(M) = SUM_i frac{Q_I}{epsint * r_iM} * EXP{frac{-kappa}{r_iM}

      
  os << "# Entering FDPoissonBoltzmann::DebyeHueckel ..." << endl;
        


  int size = atoms.size();
        
  double PI4 = 4 * ppp.getPI();
  double gridspacing_sq = gridspacing * gridspacing;
  double TWOgridspacing = 2 * gridspacing;



  //loop over atoms
  for (int i=0; i < size; ++i) {

    double factor = atoms.charge(i)/( PI4* (1 + kappa * radii[i]));

    //do the the k=0 face


    int K=0, J, I;
    double dvec[3];
    double dvec_sq[3];
    for (int j=0; j <3; ++j) {
      dvec[j] = gridstart[j] -(atoms.pos(i))[j];
      dvec_sq[j] = dvec[j] * dvec[j];
    }


    for (J=0; J < GPY; ++J) {
      dvec_sq[1] += TWOgridspacing * dvec[1] + gridspacing_sq;
      dvec[1]    += gridspacing;
      //BEWARE: this is not a programming error! (Mika)
      dvec_sq[0] =(atoms.pos(i))[0] *(atoms.pos(i))[0];
      double dvecc = dvec[0];
      double dvecc_sq = dvecc * dvecc;
      double dvec_sqYplusZ = dvec_sq[1] + dvec_sq[2];
      for (I=0; I < GPX; ++I) {
	dvecc_sq += TWOgridspacing * dvecc + gridspacing_sq;
	dvecc    += gridspacing;
	double dist = sqrt(dvecc_sq + dvec_sqYplusZ);
	rhogrid[(index(I,J,K))] += factor * exp(kappa * (radii[i] - dist))/dist;
      }
    }

    //do the the k=gridpointsK-1 face

    K = GPZ-1;
    for (int j=0; j < 3; ++j) {
      dvec[j] = gridstart[j] -(atoms.pos(i))[j];
      dvec_sq[j] = dvec[j] * dvec[j];
    }
    dvec[2] = gridstart[2] -(atoms.pos(i))[2] + (GPZ+1) * gridspacing;
    dvec_sq[2] = dvec[2] * dvec[2];

    for (J=0; J < GPY; ++J) {
      dvec_sq[1] += TWOgridspacing * dvec[1] + gridspacing_sq;
      dvec[1]    += gridspacing;
      //BEWARE: this is not a programming error! (Mika)
      dvec_sq[0] = (atoms.pos(i))[0] *(atoms.pos(i))[0];
      double dvecc = dvec[0];
      double dvecc_sq = dvecc * dvecc;
      double dvec_sqYplusZ = dvec_sq[1] + dvec_sq[2];
      for (I=0; I < GPX; ++I) {
	dvecc_sq += TWOgridspacing * dvecc + gridspacing_sq;
	dvecc    += gridspacing;
	double dist = sqrt(dvecc_sq + dvec_sqYplusZ);

	rhogrid[(index(I,J,K))] += factor * exp(kappa * (radii[i] - dist))/dist;
      }
    }

        

    //do the I = 0 face
    I = 0;
    for (int j=0; j < 3 ;++j){
      dvec[j] = gridstart[j] -(atoms.pos(i))[j];
      dvec_sq[j] = dvec[j] * dvec[j];
    }



    for (K=0; K < GPZ; ++K) {
      dvec_sq[2] += TWOgridspacing * dvec[2] + gridspacing_sq;
      dvec[2]    += gridspacing;
      //BEWARE: this is not a programming error! (Mika)
      dvec_sq[1] =(atoms.pos(i))[1] *(atoms.pos(i))[1];
      double dvecc = dvec[1];
      double dvecc_sq = dvecc * dvecc;
      double dvec_sqZplusX = dvec_sq[2] + dvec_sq[0];
      for (J=0; J < GPY; ++J) {
	dvecc_sq += TWOgridspacing * dvecc + gridspacing_sq;
	dvecc    += gridspacing;
	double dist = sqrt(dvecc_sq + dvec_sqZplusX);
	rhogrid[(index(I,J,K))] += factor * exp(kappa * (radii[i] - dist))/dist;
      }
    }
    //do the I = gridpointsI-1 face
    I = GPX-1;
    for (int j=0; j < 3; ++j){
      dvec[j] = gridstart[j] -(atoms.pos(i))[j];
      dvec_sq[j] = dvec[j] * dvec[j];
    }
    dvec[0] = gridstart[0] -(atoms.pos(i))[0] + (GPX + 1) * gridspacing;
    dvec_sq[0] = dvec[0] * dvec[0];

    for (K=0; K < GPZ; ++K) {
      dvec_sq[2] += TWOgridspacing * dvec[2] + gridspacing_sq;
      dvec[2]    += gridspacing;
      //BEWARE: this is not a programming error! (Mika)
      dvec_sq[1] =(atoms.pos(i))[1] *(atoms.pos(i))[1];
      double dvecc = dvec[1];
      double dvecc_sq = dvecc * dvecc;
      double dvec_sqZplusX = dvec_sq[2] + dvec_sq[0];
      for (J=0; J < GPY; ++J) {
	dvecc_sq += TWOgridspacing * dvecc + gridspacing_sq;
	dvecc    += gridspacing;
	double dist = sqrt(dvecc_sq + dvec_sqZplusX);
	rhogrid[(index(I,J,K))] += factor * exp(kappa * (radii[i] - dist))/dist;
      }
    }


    //do the J = 0 face
    J = 0;
    for (int j=0; j <3;++j) {
      dvec[j] = gridstart[j] -(atoms.pos(i))[j];
      dvec_sq[j] = dvec[j] * dvec[j];
    }

    for (K=0; K < GPZ; ++K) {
      dvec_sq[2] += TWOgridspacing * dvec[2] + gridspacing_sq;
      dvec[2]    += gridspacing;
      //BEWARE: this is not a programming error!  (Mika)
      dvec_sq[0] =(atoms.pos(i))[0]*(atoms.pos(i))[0];
      double dvecc = dvec[0];
      double dvecc_sq = dvecc * dvecc;
      double dvec_sqZplusY = dvec_sq[2] + dvec_sq[1];
      for (I=0; I < GPX; ++I) {
	dvecc_sq += TWOgridspacing * dvecc + gridspacing_sq;
	dvecc    += gridspacing;
	double dist = sqrt(dvecc_sq + dvec_sqZplusY);
	rhogrid[(index(I,J,K))] += factor * exp(kappa * (radii[i] - dist))/dist;
      }
    }


    //do the J = gridpointsJ-1 face
    J = GPY - 1;
    for (int j=0; j < 3; ++j) {
      dvec[j] = gridstart[j] -(atoms.pos(i))[j];
      dvec_sq[j] = dvec[j] * dvec[j];
    }


    dvec[1] = gridstart[1] -(atoms.pos(i))[1] + (GPY + 1) * gridspacing;
    dvec_sq[1] = dvec[1] * dvec[1];

    for (K=0; K < GPZ; ++K) {
      dvec_sq[2] += TWOgridspacing * dvec[2] + gridspacing_sq;
      dvec[2]    += gridspacing;
      //BEWARE: this is not a programming error!  (Mika)
      dvec_sq[0] =(atoms.pos(i))[0]*(atoms.pos(i))[0];
      double dvecc = dvec[0];
      double dvecc_sq = dvecc * dvecc;
      double dvec_sqZplusY = dvec_sq[2] + dvec_sq[1];
      for (I=0; I < GPX; ++I) {
	dvecc_sq += TWOgridspacing * dvecc + gridspacing_sq;
	dvecc    += gridspacing;
	double dist = sqrt(dvecc_sq + dvec_sqZplusY);
	rhogrid[(index(I,J,K))] += factor * exp(kappa * (radii[i] - dist))/dist;
      }
    }
  } //end loop over atoms charges


} //end  DebyeHueckel




