#include <iostream>
#include <stdio.h>
#include <string>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <cassert>
#include "../args/Arguments.h"
#include "../args/BoundaryParser.h"
#include "../args/GatherParser.h"
#include "../bound/Boundary.h"
#include "../gio/Ginstream.h"
#include "../gio/InG96.h"
#include "../gcore/System.h"
#include "../gcore/Molecule.h"
#include "../gcore/MoleculeTopology.h"
#include "../gcore/Solvent.h"
#include "../gcore/SolventTopology.h"
#include "../gcore/AtomTopology.h"
#include "../gmath/Vec.h"

#include "AtomSpecifier.h"
#include "Hbondcalc.h"
#include "Hbond.h"
#include "Hbond3c.h"
#include "Neighbours.h"


using namespace args;
using namespace gio;
using namespace gcore;
using namespace bound;
using namespace std;

using gcore::System;
using args::Arguments;
using utils::Hbondcalc;
using utils::AtomSpecifier;


void Hbondcalc::readinmasses(std::string filename)
{
  
  //new read in stuff, using the Ginstream...
  Ginstream nf(filename);
  vector<string> buffer;
  nf.getblock(buffer);

  if(buffer[0]!="HYDROGENMASS")
    throw gromos::Exception("Hbondcalc","Mass file does not contain a HYDROGENMASS block!");
  
  istringstream is; 
  double mass = 0;
  
  //read in the hydrogen masses...
  for(unsigned int j=1; j< buffer.size()-1; j++){
    is.clear();
    is.str(buffer[j]);
    is >> mass;
    d_mass_hydrogens.push_back(mass);
  }

  //get the ACCEPTORMASS block
  nf.getblock(buffer);
  
  if(buffer[0]!="ACCEPTORMASS")
    throw gromos::Exception("Hbondcalc","Mass file does not contain a ACCEPTORMASS block!");

   //read in the acceptor masses...
  for(unsigned int j=1; j< buffer.size()-1; j++){
    is.clear();
    is.str(buffer[j]);
    is >> mass;
    d_mass_acceptors.push_back(mass);
  }
  /*
  for (unsigned int i=0; i < d_mass_hydrogens.size(); ++i) cout << "H: " << d_mass_hydrogens[i] << endl;
  for (unsigned int i=0; i < d_mass_acceptors.size(); ++i) cout << "A: " << d_mass_acceptors[i] << endl;
  */
}//end readinmasses

 


void Hbondcalc::determineAtoms()
{
  // we allow for the definition of four groups of atoms
  // 1. Donor atoms of group A
  // 2. Acceptor atoms of group A
  // 3. Donor atoms of group B
  // 4. Acceptor atoms of group B
  
  // everything can contain solvent. So to find out how many of those we have,
  // read a frame
  readframe();
  
  Arguments::const_iterator iter=d_args -> lower_bound("DonorAtomsA");
  Arguments::const_iterator to=d_args -> upper_bound("DonorAtomsA");
  
  for(;iter!=to;iter++){
    string spec=iter->second.c_str();
    d_donors.addSpecifier(spec);
  }
  
  iter=d_args -> lower_bound("AcceptorAtomsA");
  to=d_args -> upper_bound("AcceptorAtomsA");
  
  for(;iter!=to;iter++){
    string spec=iter->second.c_str();
    d_acceptors.addSpecifier(spec);
  }

  //sort them and find the atoms bound to the donor
  d_donors.sort();
  d_acceptors.sort();

  int m, a;
  for (int i=0; i < d_donors.size(); ++i) {
    m = d_donors.mol(i);
    a = d_donors.atom(i);
    if(m<0){
      int j=a % d_sys->sol(0).topology().numAtoms();
      Neighbours neigh(*d_sys,0,j,0);
      d_bound.addAtomStrict(-1,a-j+neigh[0]);
    }
    else{
      Neighbours neigh(*d_sys,m,a); 
      d_bound.addAtomStrict(m, neigh[0]);
    }
  }
  // store how many acceptors and donor we have in A
  d_num_A_donors=d_donors.size();
  d_num_A_acceptors=d_acceptors.size();
  
  // if there is no B specified, we take A
  if(d_args->count("DonorAtomsB")<=0 && d_args->count("AcceptorAtomsB")<=0){
    for(int i=0; i<d_num_A_donors; i++){
      d_donors.addAtomStrict(d_donors.mol(i), d_donors.atom(i));
      d_bound.addAtomStrict(d_bound.mol(i), d_bound.atom(i));
    }
    for(int i=0; i<d_num_A_acceptors; i++){
      d_acceptors.addAtomStrict(d_acceptors.mol(i), d_acceptors.atom(i));
    }
  }
  else{
    AtomSpecifier donor_B(*d_sys);
    AtomSpecifier bound_B(*d_sys);
    AtomSpecifier acceptor_B(*d_sys);
    
    iter=d_args -> lower_bound("DonorAtomsB");
    to=d_args -> upper_bound("DonorAtomsB");
    for(;iter!=to;iter++){
      string spec=iter->second.c_str();
      donor_B.addSpecifier(spec);
    }
    iter=d_args -> lower_bound("AcceptorAtomsB");
    to=d_args -> upper_bound("AcceptorAtomsB");
    for(;iter!=to; iter++){
      string spec=iter->second.c_str();
      acceptor_B.addSpecifier(spec);
    }


    
    // and sort these as well and find the hydrogens bound to donor_B
    donor_B.sort();
    acceptor_B.sort();
    for (int i=0; i < donor_B.size(); ++i) {
      m = donor_B.mol(i);
      a = donor_B.atom(i);
      if(m<0){
	int j=a % d_sys->sol(0).topology().numAtoms();
	Neighbours neigh(*d_sys,0,j,0);
	bound_B.addAtomStrict(-1,a-j+neigh[0]);
      }
      else{
	Neighbours neigh(*d_sys,m,a); 
	bound_B.addAtomStrict(m, neigh[0]);
      }
    }

    // copy them into the d_donors, d_bound, d_acceptors
    for(int i=0; i<donor_B.size(); i++){
      d_donors.addAtomStrict(donor_B.mol(i), donor_B.atom(i));
      d_bound.addAtomStrict(bound_B.mol(i), bound_B.atom(i));
    }
    for(int i=0; i<acceptor_B.size(); i++){
      d_acceptors.addAtomStrict(acceptor_B.mol(i), acceptor_B.atom(i));
    }

  }
} //end Hbondcalc::determineAtoms()


void Hbondcalc::determineAtomsbymass()
{
  bool keep=false;
  
  //donors
  for (int i=0; i < d_donors.size(); ++i) {
    keep=false;
    for (unsigned int j=0; j < d_mass_hydrogens.size(); ++j) {      
      if (d_donors.mass(i) == d_mass_hydrogens[j]) {
	keep=true;
	break;
      }
    }
    if(!keep){ 
      d_donors.removeAtom(i);
      d_bound.removeAtom(i);
      if(i<d_num_A_donors) d_num_A_donors--;
      i--;	
    }
  }
  
  // acceptors 
  for (int i=0; i < d_acceptors.size(); ++i) {
    
    keep=false;
    for (unsigned int j=0; j < d_mass_acceptors.size(); ++j) {
      if (d_acceptors.mass(i) == d_mass_acceptors[j]) {
	keep=true;
	break;
      }
    }
    if(!keep){
      d_acceptors.removeAtom(i);
      if(i<d_num_A_acceptors) d_num_A_acceptors--;
      i--;
    }
  }
} //end Hbondcalc::determineAtomsbymass()

void Hbondcalc::init()
{
  d_pbc = BoundaryParser::boundary(*d_sys, *d_args);  
  
  d_frames = 0, d_numHB = 0, d_numHB3c=0;
}

void Hbondcalc::calc()
{
  ++d_frames; d_numHB = 0;
  
  // loop over possible hydrogen bonds
  // first A -> B
  for(int i=0; i<d_num_A_donors; i++){
    for(int j=d_num_A_acceptors; j<d_acceptors.size(); j++){
      // check if j is bound to i
      if(d_bound.atom(i)!=d_acceptors.atom(j) ||
	 d_bound.mol(i) !=d_acceptors.mol(j)){
	calculate_single(i,j);
      }
    }
  }
  // and B->A
  for(int i=d_num_A_donors; i<d_donors.size(); i++){
    for(int j=0; j< d_num_A_acceptors; j++){
      // check if j is bound to i
      if(d_bound.atom(i)!=d_acceptors.atom(j) ||
	 d_bound.mol(i) !=d_acceptors.mol(j)){
	calculate_single(i,j);

      }
    }
  }
  
  d_time += d_dt;        

} //end Hbondcalc::calc()

void Hbondcalc::calc3c()
{
  // This does the normal calculation, including the search for 
  // 3 center hydrogen bonds.
  ++d_frames; d_numHB = 0, d_numHB3c = 0;
  double d2_1=0, d2_2=0;
  
  // loop over possible hydrogen bonds
  // first A -> B
  for(int i=0; i<d_num_A_donors; i++){
    for(int j=d_num_A_acceptors; j<d_acceptors.size(); j++){
      // check if j is bound to i
      if(d_bound.atom(i)!=d_acceptors.atom(j) ||
	 d_bound.mol(i) !=d_acceptors.mol(j)){

	// for the 3c HB we calculate the distance here
	*d_acceptors.coord(j) = 
	  d_pbc->nearestImage(d_donors.pos(i), d_acceptors.pos(j),
			      d_sys -> box());
	d2_1 = (d_donors.pos(i) - d_acceptors.pos(j)).abs2();
	if(d2_1 <= d_maxdist2){
	  // a bit sad that we will calculate the distance in there again
	  calculate_single(i,j);
	}

	if(d2_1 <= d_maxdist3c2){
	  // go into a second loop
	  for(int k=j+1; k<d_acceptors.size(); k++){
	    // check if k is bound to i
	    if(d_bound.atom(i)!=d_acceptors.atom(k) ||
	       d_bound.mol(i) !=d_acceptors.mol(k)){
	      *d_acceptors.coord(k) = 
		d_pbc->nearestImage(d_donors.pos(i), d_acceptors.pos(k),
				    d_sys->box());
	      d2_2 = (d_donors.pos(i) - d_acceptors.pos(k)).abs2();
	      
	      if(d2_2 <= d_maxdist3c2){
		calculate_single3c(i,j,k,d2_1,d2_2);
	      }
	    }
	  }
	}
      }
    }
  }
  // and B->A
  for(int i=d_num_A_donors; i<d_donors.size(); i++){
    for(int j=0; j< d_num_A_acceptors; j++){
      // check if j is bound to i
      if(d_bound.atom(i)!=d_acceptors.atom(j) ||
	 d_bound.mol(i) !=d_acceptors.mol(j)){
	
	// for the 3c HB we calculate the distance here
	*d_acceptors.coord(j) = 
	  d_pbc->nearestImage(d_donors.pos(i), d_acceptors.pos(j),
			      d_sys -> box());
	d2_1 = (d_donors.pos(i) - d_acceptors.pos(j)).abs2();
	if(d2_1 <= d_maxdist2){
	  // a bit sad that we will calculate the distance in there again
	  calculate_single(i,j);
	}
	if(d2_1 <= d_maxdist3c2){
	  // go into a second loop
	  for(int k=j+1; k<d_num_A_acceptors; k++){
	    // check if k is bound to i
	    if(d_bound.atom(i)!=d_acceptors.atom(k) ||
	       d_bound.mol(i) !=d_acceptors.mol(k)){
	      *d_acceptors.coord(k) = 
		d_pbc->nearestImage(d_donors.pos(i), d_acceptors.pos(k),
				    d_sys->box());
	      d2_2 = (d_donors.pos(i) - d_acceptors.pos(k)).abs2();

	      if(d2_2 <= d_maxdist3c2){
		calculate_single3c(i,j,k,d2_1,d2_2);
	      }
	    }
	  }
	}
      }
    }
  }
  
  d_time += d_dt;        

} //end Hbondcalc::calc3c()

void Hbondcalc::calc_native()
{
  ++d_frames; d_numHB=0;
  // loop over hydrogen bonds that we already have only
  std::map<int,Hbond>::const_iterator iter=d_hbonds.begin(), to=d_hbonds.end();
  for(; iter!=to;++iter){
    int i=iter->first/d_acceptors.size();
    int j=iter->first%d_acceptors.size();
    if(d_bound.atom(i)!=d_acceptors.atom(j) ||
       d_bound.mol(i) !=d_acceptors.mol(j)){
      calculate_single(i,j);
    }
  }
  d_time += d_dt;
}

void Hbondcalc::clear()
{
  d_frames--;
  tstime.resize(0);
  tsnum.resize(0);
  
   // loop over hydrogen bonds that we already have only
  std::map<int,Hbond>::iterator iter=d_hbonds.begin(), to=d_hbonds.end();
  for(; iter!=to;++iter){
    iter->second.clear();
  }
}


void Hbondcalc::calculate_single3c(int i, int j, int k, double d2_1, double d2_2){

  // the distances already match!!

  gmath::Vec tmpA, tmpB, tmpC, p1, p2, p3;
  double angle1=0, angle2=0, angle3=0;
  double dihedral=0;
  
  // calculate angle 1
  *d_acceptors.coord(j) = 
    d_pbc->nearestImage(d_donors.pos(i), d_acceptors.pos(j),
			d_sys->box());
  tmpA = d_acceptors.pos(j) - d_donors.pos(i);
  *d_bound.coord(i) = 
    d_pbc->nearestImage(d_donors.pos(i), d_bound.pos(i), 
			d_sys -> box());
  tmpB = d_bound.pos(i) - d_donors.pos(i);                    
  angle1 = acos((tmpA.dot(tmpB))/(tmpA.abs()*tmpB.abs()))*180/M_PI;

  // calculate angle 2  
  *d_acceptors.coord(k) = 
    d_pbc->nearestImage(d_donors.pos(i), d_acceptors.pos(k),
			d_sys->box());
  tmpA = d_acceptors.pos(k) - d_donors.pos(i);
  tmpB = d_bound.pos(i) - d_donors.pos(i);                    
  angle2 = acos((tmpA.dot(tmpB))/(tmpA.abs()*tmpB.abs()))*180/M_PI;
  
  if (angle1 >= d_minangle3c && angle2 >= d_minangle3c) { 
    // first test is passed, calculate angle3 as well
    // everything is already at the nearest image to the hydrogen
    tmpA = d_acceptors.pos(j) - d_donors.pos(i);
    tmpB = d_acceptors.pos(k) - d_donors.pos(i);
    angle3 = acos((tmpA.dot(tmpB))/(tmpA.abs()*tmpB.abs()))*180/M_PI;
    
    if((angle1 + angle2 + angle3) >= d_minanglesum3c){
      // second test is passed, calculate the dihedral
      tmpA = d_bound.pos(i) - d_acceptors.pos(j);
      tmpB = d_donors.pos(i) - d_acceptors.pos(k);
      tmpC = d_acceptors.pos(k) - d_acceptors.pos(j);

      p1 = tmpA.cross(tmpC);
      p2 = tmpB.cross(tmpC);

      dihedral = acos((p1.dot(p2))/(p1.abs()*p2.abs()))*180/M_PI;

      p3 = p1.cross(p2);
      if (p3.dot(tmpC)<0)
	dihedral = -dihedral;   
      
      if(dihedral <= d_maxdihedral3c ||
	 dihedral >= -d_maxdihedral3c){
	
	// we found a hydrogen bond!
	int asize=d_acceptors.size();
	int index=i*asize*asize + j*asize + k;
      
	d_hbonds3c[index].setIndices(i,i,j,k);
	d_hbonds3c[index].adddistances(sqrt(d2_1), sqrt(d2_2));
	d_hbonds3c[index].addangles(angle1, angle2);
	d_hbonds3c[index].addanglesum(angle1+angle2+angle3);
	d_hbonds3c[index].adddihedral(dihedral);
	d_hbonds3c[index].addnum();
	++d_numHB3c;
	tstime3c.push_back(d_time);
	tsnum3c.push_back(index);
      }
    }
  }
}

void Hbondcalc::calculate_single(int i, int j){

  gmath::Vec tmpA, tmpB;
  double angle=0, distance2=0;
  
  //check whether we meet the distance
  *d_acceptors.coord(j) = 
    d_pbc->nearestImage(*d_donors.coord(i), *d_acceptors.coord(j),
			d_sys -> box());
  distance2 = (*d_donors.coord(i) - *d_acceptors.coord(j)).abs2();

  if (distance2 <= d_maxdist2) {
    tmpA = *d_acceptors.coord(j) - *d_donors.coord(i);
    *d_bound.coord(i) = 
      d_pbc->nearestImage(*d_donors.coord(i), *d_bound.coord(i), 
			  d_sys -> box());
    tmpB = *d_bound.coord(i) - *d_donors.coord(i);                    
    angle = acos((tmpA.dot(tmpB))/(tmpA.abs()*tmpB.abs()))*180/M_PI;
    if (angle >= d_minangle) { 
      // we found a hydrogen bond!
      int index=i*d_acceptors.size()+j;
      
      d_hbonds[index].setIndices(i,i,j);
      d_hbonds[index].adddistance(sqrt(distance2));
      d_hbonds[index].addangle(angle);
      d_hbonds[index].addnum();
      ++d_numHB;
      tstime.push_back(d_time);
      tsnum.push_back(index);
    }
  }
}

void Hbondcalc::printstatistics()
{

   int count = 0;
   int i_d, i_a;
   
   vector<int> totnum, realnum;
   map<int, Hbond>::const_iterator it=d_hbonds.begin();
   map<int, Hbond>::const_iterator to=d_hbonds.end();
   
   for(int i;it!=to; ++it, ++i){

     d_hbond = it->second;
    
     if (d_hbond.num() > 0) {
       ++count;
       totnum.push_back(it->first); realnum.push_back(count);
       d_hbond.calcmean();
       i_d=d_hbond.don_ind();
       i_a=d_hbond.ac_ind();
       
       std::cout << setw(3) << count;
       if(d_donors.mol(i_d)<0) std::cout << setw(8) << " ";
       else std::cout << setw(8) << d_donors.mol(i_d)+1;
       std::cout << setw(4) << d_donors.resnum(i_d)+1
	 	 << setw(4) << d_donors.resname(i_d)
		 << setw(2) << "-";
       if(d_acceptors.mol(i_a)<0) std::cout << setw(4) << " ";
       else std::cout << setw(4) << d_acceptors.mol(i_a)+1;
       std::cout << setw(4) << d_acceptors.resnum(i_a)+1
	 	 << setw(4) << d_acceptors.resname(i_a)
		 << setw(6) << d_bound.atom(i_d)+1 
		 << setw(4) << d_bound.name(i_d) 
		 << setw(2) << "-"
		 << setw(6) << d_donors.atom(i_d)+1  
		 << setw(4) << d_donors.name(i_d) 
		 << setw(2) << "-"
		 << setw(6) << d_acceptors.atom(i_a)+1 
		 << setw(4) << d_acceptors.name(i_a);
       std::cout.precision(3); 
       std::cout << setw(8) << d_hbond.meandist();
       std::cout.precision(3);
       std::cout << setw(8) << d_hbond.meanangle();
       std::cout.precision(0);
       std::cout << setw(8) << d_hbond.num();
       std::cout.setf(ios::floatfield, ios_base::fixed);
       std::cout.precision(2);             
       std::cout << setw(8) << ((d_hbond.num()/ (double) d_frames)*100)
		 << endl;       
     } // if end
   }
   
   //sort the Hbts.out according to the output
   std::vector<int>::const_iterator iter;
   for (unsigned int i=0; i < tsnum.size(); ++i) {
     int find = tsnum[i];
     iter = std::find(totnum.begin(), totnum.end(), find);
     
     timeseriesHB.precision(10);
     timeseriesHB << setw(10) << tstime[i] << " ";
     timeseriesHB << setw(10) << realnum[iter - totnum.begin()] << endl;
   }
   

} //end printstatistics

void Hbondcalc::printstatistics3c()
{

   int count = 0;
   int i_d, i_a1, i_a2;
   
   vector<int> totnum, realnum;
   map<int, Hbond3c>::const_iterator it=d_hbonds3c.begin();
   map<int, Hbond3c>::const_iterator to=d_hbonds3c.end();
   
   for(int i;it!=to; ++it, ++i){

     utils::Hbond3c hb3c = it->second;
    
     if (hb3c.num() > 0) {
       ++count;
       totnum.push_back(it->first); realnum.push_back(count);
       hb3c.calcmean();
       i_d =hb3c.don_ind();
       i_a1=hb3c.ac1_ind();
       i_a2=hb3c.ac2_ind();
       
       std::cout << setw(3) << count;
       if(d_donors.mol(i_d)<0) std::cout << setw(8) << " ";
       else std::cout << setw(8) << d_donors.mol(i_d)+1;
       std::cout << setw(4) << d_donors.resnum(i_d)+1
	 	 << setw(4) << d_donors.resname(i_d)
		 << setw(2) << "-";
       if(d_acceptors.mol(i_a1)<0) std::cout << setw(4) << " ";
       else std::cout << setw(4) << d_acceptors.mol(i_a1)+1;
       std::cout << setw(4) << d_acceptors.resnum(i_a1)+1
	 	 << setw(4) << d_acceptors.resname(i_a1);
       std::cout << setw(6) << d_bound.atom(i_d)+1 
		 << setw(4) << d_bound.name(i_d) 
		 << setw(2) << "-"
		 << setw(6) << d_donors.atom(i_d)+1  
		 << setw(4) << d_donors.name(i_d) 
		 << setw(2) << "-"
		 << setw(6) << d_acceptors.atom(i_a1)+1 
		 << setw(4) << d_acceptors.name(i_a1);
       
       std::cout.precision(3); 
       std::cout << setw(8) << hb3c.meandist_don_a1();
       std::cout.precision(3);
       std::cout << setw(8) << hb3c.meanangle_b_don_a1();

       // get the occurrence of this HB as two centered HB
       int index=d_acceptors.size() * i_d + i_a1;
       int occur1=d_hbonds[index].num();
       
       std::cout.precision(0);
       std::cout << setw(8) << occur1;
       std::cout.setf(ios::floatfield, ios_base::fixed);
       std::cout.precision(2);             
       std::cout << setw(8) << ((occur1/ (double) d_frames)*100);
       std::cout.precision(3);
       std::cout << setw(10) << hb3c.meanangle_sum();
       std::cout.precision(3);
       std::cout << setw(8) << hb3c.meandihedral();
       std::cout.precision(0);
       std::cout << setw(8) << hb3c.num();
       std::cout.setf(ios::floatfield, ios_base::fixed);
       std::cout.precision(2);             
       std::cout << setw(8) << ((hb3c.num()/ (double) d_frames)*100)
		 << endl;       
       // and the second line
       std::cout << setw(19) << " "
		 << setw(2) <<  "\\";
       if(d_acceptors.mol(i_a2)<0) std::cout << setw(4) << " ";
       else std::cout << setw(4) << d_acceptors.mol(i_a2)+1;
       std::cout << setw(4) << d_acceptors.resnum(i_a2)+1
	 	 << setw(4) << d_acceptors.resname(i_a2)
		 << setw(22) << " "
		 << setw(2) << "\\"
		 << setw(6) << d_acceptors.atom(i_a1)+1 
		 << setw(4) << d_acceptors.name(i_a1);

       std::cout.precision(3);
       std::cout << setw(8) << hb3c.meandist_don_a2();
       std::cout.precision(3);
       std::cout << setw(8) << hb3c.meanangle_b_don_a2();
       // get the occurrence of this HB as two centered HB
       index=d_acceptors.size() * i_d + i_a2;
       int occur2=d_hbonds[index].num();
       
       std::cout.precision(0);
       std::cout << setw(8) << occur2;
       std::cout.setf(ios::floatfield, ios_base::fixed);
       std::cout.precision(2);             
       std::cout << setw(8) << ((occur2/ (double) d_frames)*100);
     } // if end
   }
   
   //sort the Hbts.out according to the output
   std::vector<int>::const_iterator iter;
   for (unsigned int i=0; i < tsnum3c.size(); ++i) {
     int find = tsnum3c[i];
     iter = std::find(totnum.begin(), totnum.end(), find);
     
     timeseriesHB3c.precision(10);
     timeseriesHB3c << setw(10) << tstime3c[i] << " ";
     timeseriesHB3c << setw(10) << realnum[iter - totnum.begin()] << endl;
   }
   

} //end printstatistics3c
    
void Hbondcalc::setmaxdist(double i){
  d_maxdist2 = i*i;
}

void Hbondcalc::setminangle(double i){
  d_minangle = i;
}

void Hbondcalc::setmaxdist3c(double i){
  d_maxdist3c2 = i*i;
}
void Hbondcalc::setminangle3c(double i){
  d_minangle3c = i;
}
void Hbondcalc::setminanglesum3c(double i){
  d_minanglesum3c = i;
}
void Hbondcalc::setmaxdihedral3c(double i)
{
  d_maxdihedral3c = i;
}

void Hbondcalc::settime(double i, double j)
{
  d_time = i;
  d_dt = j;
}

void Hbondcalc::opents(string fi1, string fi2)
{
  timeseriesHB.open(fi1.c_str());
  timeseriesHBtot.open(fi2.c_str());
}
void Hbondcalc::opents3c(string fi1, string fi2)
{
  timeseriesHB3c.open(fi1.c_str());
  timeseriesHB3ctot.open(fi2.c_str());
}

   
void Hbondcalc::readframe()
{
  
  InG96 icc;
  
  try{
    d_args -> check("ref",1);
    Arguments::const_iterator iterr=d_args -> lower_bound("ref");
    icc.open((iterr->second).c_str());
  }
  catch(const Arguments::Exception &){
    d_args -> check("traj",1);
    Arguments::const_iterator iterr=d_args -> lower_bound("traj");
    icc.open((iterr->second).c_str());
  }
  icc.select("ALL");
  icc >> *d_sys;
  icc.close();
  
}

void Hbondcalc::writets()
{
  
  timeseriesHBtot.precision(6);
  timeseriesHBtot << setw(10) << d_time-d_dt;
  timeseriesHBtot.precision(5);
  timeseriesHBtot << setw(10) << d_numHB << endl;
  if(d_do3c){
    timeseriesHB3ctot.precision(6);
    timeseriesHB3ctot << setw(10) << d_time-d_dt;
    timeseriesHB3ctot.precision(5);
    timeseriesHB3ctot << setw(10) << d_numHB3c << endl;
  }
}


Hbondcalc::Hbondcalc(gcore::System &sys, args::Arguments &args)
{
  d_sys=&sys;
  d_args=&args;
  d_donors = AtomSpecifier(sys);
  d_bound  = AtomSpecifier(sys);
  d_acceptors = AtomSpecifier(sys);

  if(args.count("threecenter")>=0) d_do3c=true;
  else d_do3c=false;
  
  //open timeseries file
  opents("Hbts.out", "Hbnumts.out");
  if(d_do3c){
    opents3c("Hb3cts.out", "Hb3cnumts.out");
  }
}
