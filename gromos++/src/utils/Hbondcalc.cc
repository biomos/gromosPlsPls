#include <iostream>
#include <stdio.h>
#include <string>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include "AtomSpecifier.h"
#include "Hbondcalc.h"
#include "Hbond.h"
#include "Neighbours.h"
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


using namespace args;
using namespace gio;
using namespace gcore;
using namespace bound;
using namespace std;

using gcore::System;
using args::Arguments;
using utils::Hbondcalc;
using utils::AtomSpecifier;


void Hbondcalc::readinmasses(std::string fi)
{
  char buffer[245];
  ifstream massfile;
  massfile.open(fi.c_str());

  string dum;

  while (!massfile.eof()) {
    
     massfile.getline(buffer,100);
     dum = string(buffer);
     string hdum = dum;
     string adum = dum;

     if (dum == "HYDROGENMASS") {
         massfile.getline(buffer,100);
         hdum = string(buffer);
       while (hdum != "END") { 	 
         d_mass_hydrogens.push_back(atof(hdum.c_str()));
         massfile.getline(buffer,100);
         hdum = string(buffer);
       }
     }

     if (dum == "ACCEPTORMASS") {
       massfile.getline(buffer,100);
       adum = string(buffer);
       while (adum != "END") { 	 
         d_mass_acceptors.push_back(atof(adum.c_str()));
         massfile.getline(buffer,100);
         adum = string(buffer);
       }
     }
  
  } //end while massfile


} //end readinmasses


void Hbondcalc::determineAtoms()
{
 Arguments::const_iterator iter=d_args -> lower_bound("SoluteDonorAtoms");
 Arguments::const_iterator to=d_args -> upper_bound("SoluteDonorAtoms");

 for(;iter!=to;iter++){
   string spec=iter->second.c_str();
   d_donors.addSpecifier(spec);
 }

 iter=d_args -> lower_bound("SoluteAcceptorAtoms");
 to=d_args -> upper_bound("SoluteAcceptorAtoms");

 for(;iter!=to;iter++){
   string spec=iter->second.c_str();
   d_acceptors.addSpecifier(spec);
 }

  
 //sort them, create the list of atoms that the hydrogens are bound to
 d_donors.sort();
 d_acceptors.sort();

 int m, a;
 for (int i=0; i < (d_donors).size(); ++i) {
  m = (d_donors).mol(i);
  a = (d_donors).atom(i);
  Neighbours neigh(*d_sys,m,a); 
  d_donors_bound_to.addAtomStrict(m, neigh[0]);
 }
 
    
 try{
     d_args -> check("SolventAcceptorAtoms");
     { 

      //read initial frame
      readframe();
     
      iter=d_args -> lower_bound("SolventAcceptorAtoms");
      to=d_args -> upper_bound("SolventAcceptorAtoms");

      for(;iter!=to;iter++){
       string spec=iter->second.c_str();
       d_acceptors_solv.addSpecifier(spec);
      }
   
     }
   
   }
   catch(Arguments::Exception e){
   }


   try{
   d_args -> check("SolventDonorAtoms");
   { 

     //read initial frame
     readframe();

    iter=d_args -> lower_bound("SolventDonorAtoms");
    to=d_args -> upper_bound("SolventDonorAtoms");

    for(;iter!=to;iter++){
     string spec=iter->second.c_str();
     d_donors_solv.addSpecifier(spec);
    }
         
    for (int i=0; i < (d_donors_solv).size(); ++i) {
     m = (d_donors_solv).mol(i);
     a = (d_donors_solv).atom(i);
     string name = (d_donors_solv).name(a);
     for (int j=0;j< d_sys -> sol(0).topology().numAtoms();++j) {
       if (name ==  d_sys -> sol(0).topology().atom(j).name()) {
        Neighbours neigh(*d_sys,0,j,0);
        d_donors_bound_to_solv.addAtomStrict(-1, a + (neigh[0] - j));
       }
     }      
    }
   
   }

   }
   catch(Arguments::Exception e){
   }
 


} //end Hbondcalc::determineAtoms()


void Hbondcalc::determineAtomsbymass()
{

  int m, a;
  //check whether specified atoms are indeed hydrogens
  
  AtomSpecifier donor_tmp(*d_sys);
  AtomSpecifier donor_bound_to_tmp(*d_sys);
  AtomSpecifier acceptor_tmp(*d_sys);
  AtomSpecifier donor_solv_tmp(*d_sys);
  AtomSpecifier donor_bound_to_solv_tmp(*d_sys);
  AtomSpecifier acceptor_solv_tmp(*d_sys);


  //solute donors
  for (int i=0; i < (d_donors).size(); ++i) {
      m = (d_donors).mol(i);
      a = (d_donors).atom(i);
    for (int j=0; j < (int) (d_mass_hydrogens).size(); ++j) {      
      if (d_sys -> mol(m).topology().atom(a).mass() == d_mass_hydrogens[j]) { 
        donor_tmp.addAtomStrict(m,a);
        int tmpM = d_donors_bound_to.mol(i);
        int tmpA = d_donors_bound_to.atom(i);
        donor_bound_to_tmp.addAtomStrict(tmpM,tmpA);
      }
    }
  }
  
  //check whether specified atoms are indeed acceptor atoms
  for (int i=0; i < (d_acceptors).size(); ++i) {
      m = (d_acceptors).mol(i);
      a = (d_acceptors).atom(i);
    for (int j=0; j < (int) (d_mass_acceptors).size(); ++j) {
      if (d_sys -> mol(m).topology().atom(a).mass() == d_mass_acceptors[j]) {
        acceptor_tmp.addAtomStrict(m,a);
      }
    }
  }


  try{
   d_args -> check("SolventDonorAtoms");
   { 
     
    //solvent donors
    for (int i=0; i < (d_donors_solv).size(); ++i) {
      m = (d_donors_solv).mol(i);
      a = (d_donors_solv).atom(i);
      
      //map to solvent topology
      string name = (d_donors_solv).name(a);
      for (int j=0;j< d_sys -> sol(0).topology().numAtoms();++j) {
       if (name ==  d_sys -> sol(0).topology().atom(j).name()) { 
	 
        for (int k=0; k < (int) (d_mass_hydrogens).size(); ++k) {      
         if (d_sys -> sol(0).topology().atom(j).mass() == d_mass_hydrogens[k]) {
          donor_solv_tmp.addAtomStrict(-1,a);
          int tmpA = d_donors_bound_to_solv.atom(i);
          donor_bound_to_solv_tmp.addAtomStrict(-1,tmpA);
	 }
	}
       }
      }
    }
  

   }

   }
   catch(Arguments::Exception e){
   }

  try{
   d_args -> check("SolventAcceptorAtoms");
   { 
    //check whether specified atoms are indeed acceptor atoms; this is for solvent
    for (int i=0; i < (d_acceptors_solv).size(); ++i) {
      m = (d_acceptors_solv).mol(i);
      a = (d_acceptors_solv).atom(i);
      //map to solvent topology
      string name = (d_acceptors_solv).name(a);
      for (int j=0;j< d_sys -> sol(0).topology().numAtoms();++j) {
       if (name ==  d_sys -> sol(0).topology().atom(j).name()) { 
        for (int k=0; k < (int) (d_mass_acceptors).size(); ++k) {      
         if (d_sys -> sol(0).topology().atom(j).mass() == d_mass_acceptors[k]) {
          acceptor_solv_tmp.addAtomStrict(-1,a);
	 }
	}
       }
      }
    }
  
   }
   }
   catch(Arguments::Exception e){
   } 

  
  d_donors = donor_tmp;
  d_donors_bound_to = donor_bound_to_tmp;
  d_acceptors = acceptor_tmp;
  d_donors_solv = donor_solv_tmp;
  d_donors_bound_to_solv = donor_bound_to_solv_tmp;
  d_acceptors_solv = acceptor_solv_tmp;

     
} //end Hbondcalc::determineAtomsbymass()

void Hbondcalc::calcHintra_native_init()
{
 calcHintra_init();
 readframe();
 calcHintra();

 vector<Hbond> HBtmp;
 for (int i=0; i < (int) d_hbonds.size(); ++i) {
   d_hbond = d_hbonds[i];
   if (d_hbond.num() > 0) {
    d_hbond.clear();
    HBtmp.push_back(d_hbond);
   }
 }

 d_hbonds = HBtmp;


} //end Hbondcalc::calcHintra_native_init()

void Hbondcalc::calcHinter_native_init()
{
 calcHinter_init();
 readframe();
 calcHinter();

 vector<Hbond> HBtmp;
 for (int i=0; i < (int) d_hbonds.size(); ++i) {
   d_hbond = d_hbonds[i];
   if (d_hbond.num() > 0) {
    d_hbond.clear();
    HBtmp.push_back(d_hbond);
   }
 }

 d_hbonds = HBtmp;
 

} // Hbondcalc::calcHinter_native_init()

void Hbondcalc::calcHintra_init()
{

  d_pbc = BoundaryParser::boundary(*d_sys, *d_args);  

  d_frames = 0, d_numHB = 0;

  for (int i=0; i < (d_donors).size(); ++i) {
   for (int j=0; j < (d_acceptors).size(); ++j) {
     //check whether we are in the same molecule      
     if ((d_donors).mol(i) == (d_acceptors).mol(j)) {
       //remove atoms that are bound to the respective hydrogen
       if ((d_donors_bound_to).atom(i) != (d_acceptors).atom(j)) {
         d_hbond.setHBond(d_donors.mol(i),d_donors.atom(i), d_donors_bound_to.mol(i), d_donors_bound_to.atom(i),
                     d_acceptors.mol(j), d_acceptors.atom(j));
         d_hbond.setindices(i, i, j);
	 d_hbonds.push_back(d_hbond);
       }
     }
   }
  }

} //end Hbondcalc::calcHintra_init()


void Hbondcalc::calcHinter_init()
{

  d_omit_self_species = false;
  vector<string> arg;
  vector<int> start, stop;
  try{
   d_args -> check("molrange");
   {
   Arguments::const_iterator iter=d_args -> lower_bound("molrange");
   Arguments::const_iterator to=d_args -> upper_bound("molrange");

   for(;iter!=to;iter++){
    string spec=iter->second.c_str();
    arg.push_back(spec);
   }

   d_nummol = atoi((arg[0]).c_str());
   
   for (int i=1; i < (int) arg.size(); ++i) {
     string tmp = arg[i];
     int mid = tmp.find("-");
     string s1 = tmp.substr(0,mid);
     string s2 = tmp.substr(mid+1,tmp.size());
     start.push_back(atoi(s1.c_str())-1);
     stop.push_back(atoi(s2.c_str())-1);
   }
  
   d_omit_self_species = true;

   if ((int) start.size() != d_nummol) throw Hbondcalc::Exception(" Inconsistent input for @molrange! Number of molecules and range mismatch!.\n");
   if ((int) stop.size() != d_nummol) throw Hbondcalc::Exception(" Inconsistent input for @molrange! Number of molecules and range mismatch!.\n");
  
   }
   
  }
   catch(Arguments::Exception e){
   }
 
 
  d_pbc = BoundaryParser::boundary(*d_sys, *d_args);  

  d_frames = 0, d_numHB = 0;  

  for (int i=0; i < (d_donors).size(); ++i) {
   for (int j=0; j < (d_acceptors).size(); ++j) {
     //check whether we are NOT in the same molecule      
     if ((d_donors).mol(i) != (d_acceptors).mol(j)) {
       //check for identical residue names
       if(d_omit_self_species) {
	 int m1 = (d_donors).mol(i);
         int m2 = (d_acceptors).mol(j);
         int molnumdon=0;
         int molnumac=0;
         for (int k=0; k < (int) start.size(); ++k) {
	   if ((m1 >= start[k]) && (m1 <= stop[k])) molnumdon = k;
           if ((m2 >= start[k]) && (m2 <= stop[k])) molnumac = k;
	 }
         if (molnumdon != molnumac) {
          if ((d_donors_bound_to).atom(i) != (d_acceptors).atom(j)) {
           d_hbond.setHBond(d_donors.mol(i),d_donors.atom(i), d_donors_bound_to.mol(i), d_donors_bound_to.atom(i),
                     d_acceptors.mol(j), d_acceptors.atom(j));
           d_hbond.setindices(i, i, j);
	   d_hbonds.push_back(d_hbond);
	  }
	 }
       }
       else {
       //remove atoms that are bound to the respective hydrogen
        if ((d_donors_bound_to).atom(i) != (d_acceptors).atom(j)) {
         d_hbond.setHBond(d_donors.mol(i),d_donors.atom(i), d_donors_bound_to.mol(i), d_donors_bound_to.atom(i),
                     d_acceptors.mol(j), d_acceptors.atom(j));
         d_hbond.setindices(i, i, j);
	 d_hbonds.push_back(d_hbond);
	}
       }
     }
   }
  }

} //end Hbondcalc::calcHinter_init()


void Hbondcalc::calcHsolusolv_init()
{
 d_pbc = BoundaryParser::boundary(*d_sys, *d_args);  

  d_frames = 0, d_numHB = 0;

  for (int i=0; i < (d_donors).size(); ++i) {
   for (int j=0; j < (d_acceptors_solv).size(); ++j) {
     dist.push_back(0.0);
     ang.push_back(0.0);
     num.push_back(0);
   }
  }
   
  for (int i=0; i < (d_donors_solv).size(); ++i) {
   for (int j=0; j < (d_acceptors).size(); ++j) {
     dist.push_back(0.0);
     ang.push_back(0.0);
     num.push_back(0);
   }
  } 
  

} //end Hbondcalc::calcHsolusolv_init()

void Hbondcalc::calcHsolvsolv_init()
{
 d_pbc = BoundaryParser::boundary(*d_sys, *d_args);  

  d_frames = 0, d_numHB = 0;

  for (int i=0; i < (d_donors_solv).size(); ++i) {
   for (int j=0; j < (d_acceptors_solv).size(); ++j) {
    if ((d_donors_bound_to_solv).atom(i) != (d_acceptors_solv).atom(j)) {
     dist.push_back(0.0);
     ang.push_back(0.0);
     num.push_back(0);
    }
   }
  }

} //end Hbondcalc::calcHsolvsolv_init()
 
void Hbondcalc::calcHintra()
{

  ++d_frames; d_numHB = 0;

  gmath::Vec tmpA, tmpB;
  double angle=0, distance=0;
  //Hbond d_hbond;
  for (int i=0; i < (int) d_hbonds.size(); ++i) {
    d_hbond = d_hbonds[i];
    //check whether we meet the distance
    *d_acceptors.coord(d_hbond.ac_ind()) = d_pbc->nearestImage(*d_donors.coord(d_hbond.don_ind()),*d_acceptors.coord(d_hbond.ac_ind()),d_sys -> box());
    distance = (*d_donors.coord(d_hbond.don_ind()) - *d_acceptors.coord(d_hbond.ac_ind())).abs();
    if (distance <= d_maxdist) {
     tmpA = (*d_acceptors.coord(d_hbond.ac_ind()) - *d_donors.coord(d_hbond.don_ind()));
     *d_donors_bound_to.coord(d_hbond.b_ind()) = d_pbc->nearestImage(*d_donors.coord(d_hbond.don_ind()),*d_donors_bound_to.coord(d_hbond.b_ind()), d_sys -> box());
     tmpB = (*d_donors_bound_to.coord(d_hbond.b_ind()) - *d_donors.coord(d_hbond.don_ind()));                    
     angle = acos((tmpA.dot(tmpB))/(sqrt(tmpA.dot(tmpA))*(sqrt(tmpB.dot(tmpB)))))*180/3.1416;
      if (angle >= d_minangle) { 
       d_hbond.adddistance(distance);
       d_hbond.addangle(angle);
       d_hbond.addnum();
       d_hbonds[i] = d_hbond;
       ++d_numHB;
       tstime.push_back(d_time);
       tsnum.push_back(i);

      }
    }
  }
  
  //write time series
  writets();

  d_time += d_dt;        

} //end Hbondcalc::calcintra()

void Hbondcalc::calcHinter()
{
  //loop is identical, so just call Hbondcalc::calcHintra()
  calcHintra();
}

void Hbondcalc::calcHsolusolv()
{
  //this extra loop is need because of potential memory problems when pushing back
  //Hbond objects into a std::vector
  //in the future this will be replaced by a call to calcHintra()...
    ++d_frames; d_numHB = 0; int numm = -1;


  gmath::Vec tmpA, tmpB;
  double angle=0, distance=0;
  //Hbond d_hbond;
  for (int i=0; i < (d_donors).size(); ++i) {
   for (int j=0; j < (d_acceptors_solv).size(); ++j) {
     ++numm;
      //check whether we meet the distance
    *d_acceptors_solv.coord(j) = d_pbc->nearestImage(*d_donors.coord(i),*d_acceptors_solv.coord(j),d_sys -> box());
    distance = (*d_donors.coord(i) - *d_acceptors_solv.coord(j)).abs();
    if (distance <= d_maxdist) {
     tmpA = (*d_acceptors_solv.coord(j) - *d_donors.coord(i));
     *d_donors_bound_to.coord(i) = d_pbc->nearestImage(*d_donors.coord(i),*d_donors_bound_to.coord(i), d_sys -> box());
     tmpB = (*d_donors_bound_to.coord(i) - *d_donors.coord(i));                    
     angle = acos((tmpA.dot(tmpB))/(sqrt(tmpA.dot(tmpA))*(sqrt(tmpB.dot(tmpB)))))*180/3.1416;
      if (angle >= d_minangle) {
	dist[numm] = distance;
        ang[numm] = angle;
        num[numm] +=1;
       ++d_numHB;
       
      }
    }
   }
  }
  

  for (int i=0; i < (d_donors_solv).size(); ++i) {
   for (int j=0; j < (d_acceptors).size(); ++j) {
     ++numm;
      //check whether we meet the distance
    *d_acceptors.coord(j) = d_pbc->nearestImage(*d_donors_solv.coord(i),*d_acceptors.coord(j),d_sys -> box());
    distance = (*d_donors_solv.coord(i) - *d_acceptors.coord(j)).abs();
    if (distance <= d_maxdist) {
     tmpA = (*d_acceptors.coord(j) - *d_donors_solv.coord(i));
     *d_donors_bound_to_solv.coord(i) = d_pbc->nearestImage(*d_donors_solv.coord(i),*d_donors_bound_to_solv.coord(i), d_sys -> box());
     tmpB = (*d_donors_bound_to_solv.coord(i) - *d_donors_solv.coord(i));                    
     angle = acos((tmpA.dot(tmpB))/(sqrt(tmpA.dot(tmpA))*(sqrt(tmpB.dot(tmpB)))))*180/3.1416;
      if (angle >= d_minangle) { 
	dist[numm] = distance;
        ang[numm] = angle;
        num[numm] +=1;
       ++d_numHB;
       
      }
    }
   }
  }

  //write time series
  writets();
        

  d_time += d_dt;     

} //end Hbondcalc::calcHsolusolv()

void Hbondcalc::calcHsolvsolv()
{
  //this extra loop is need because of potential memory problems when pushing back
  //Hbond objects into a std::vector
  //in the future this will be replaced by a call to calcHintra()...
    ++d_frames; d_numHB = 0; int numm = -1;


  gmath::Vec tmpA, tmpB;
  double angle=0, distance=0;
  //Hbond d_hbond;
  for (int i=0; i < (d_donors_solv).size(); ++i) {
   for (int j=0; j < (d_acceptors_solv).size(); ++j) {
     //check that we are not in the same molecule!
    if ((d_donors_bound_to_solv).atom(i) != (d_acceptors_solv).atom(j)) {
     ++numm;
      //check whether we meet the distance
     *d_acceptors_solv.coord(j) = d_pbc->nearestImage(*d_donors_solv.coord(i),*d_acceptors_solv.coord(j),d_sys -> box());
     distance = (*d_donors_solv.coord(i) - *d_acceptors_solv.coord(j)).abs();
     if (distance <= d_maxdist) {
      tmpA = (*d_acceptors_solv.coord(j) - *d_donors_solv.coord(i));
      *d_donors_bound_to_solv.coord(i) = d_pbc->nearestImage(*d_donors_solv.coord(i),*d_donors_bound_to_solv.coord(i), d_sys -> box());
      tmpB = (*d_donors_bound_to_solv.coord(i) - *d_donors_solv.coord(i));                    
      angle = acos((tmpA.dot(tmpB))/(sqrt(tmpA.dot(tmpA))*(sqrt(tmpB.dot(tmpB)))))*180/3.1416;
       if (angle >= d_minangle) {
	dist[numm] = distance;
        ang[numm] = angle;
        num[numm] +=1;
        ++d_numHB;       
       }
     }
    }
   }
  }
  //write time series
  writets();

  d_time += d_dt;     

} //end Hbondcalc::calcHsolvsolv()
 

 void Hbondcalc::printstatistics()
 {

   int count = 0;
   vector<int> totnum, realnum;
   for (int i=0; i < (int) d_hbonds.size(); ++i) {

    d_hbond = d_hbonds[i];
    
    if (d_hbond.num() > 0) {
      ++count;
     totnum.push_back(i); realnum.push_back(count);
     d_hbond.calcmean();

     std::cout << count
	       << setw(8) << d_hbond.Hmol()+1 
               << setw(4) << d_sys -> mol(d_hbond.Hmol()).topology().resNum(d_hbond.Hatom())+1
 	                  << d_sys -> mol(d_hbond.Hmol()).topology().resName(d_sys -> mol(d_hbond.Hmol()).topology().resNum(d_hbond.Hatom()))
	       << setw(2) << "-" 
               << setw(4) << d_hbond.Acmol()+1 
               << setw(4) << d_sys -> mol(d_hbond.Acmol()).topology().resNum(d_hbond.Acatom())+1
	                  << d_sys -> mol(d_hbond.Acmol()).topology().resName(d_sys -> mol(d_hbond.Acmol()).topology().resNum(d_hbond.Acatom()))
               << setw(6) << d_hbond.Htoatom()+1 
               << setw(4) << d_sys -> mol(d_hbond.Htomol()).topology().atom(d_hbond.Htoatom()).name() 
               << setw(2) << "-"
               << setw(6) << d_hbond.Hatom()+1  
               << setw(4) << d_sys -> mol(d_hbond.Hmol()).topology().atom(d_hbond.Hatom()).name() 
               << setw(2) << "-"
               << setw(6) << d_hbond.Acatom()+1 
               << setw(4) << d_sys -> mol(d_hbond.Acmol()).topology().atom(d_hbond.Acatom()).name();
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
   for (int i=0; i < (int) tsnum.size(); ++i) {
     int find = tsnum[i];
     iter = std::find(totnum.begin(), totnum.end(), find);

     timeseriesHB.precision(6);
     timeseriesHB << setw(10) << tstime[i];
     timeseriesHB.precision(5);
     timeseriesHB << setw(10) << realnum[iter - totnum.begin()] << endl;
   }
   

 } //end printstatistics
    
 void Hbondcalc::printstatistics_solusolv() 
{
  int n=0;
  int c=0;
    for (int i=0; i < (d_donors).size(); ++i) {
     for (int j=0, m=0; j < (d_acceptors_solv).size(); ++j, ++m) {
       ++c;
 
      if (num[c] > 0) {
       ++n;

       int solvatom = 0;
       string name = (d_acceptors_solv).name(j);
       for (int k=0;k< d_sys -> sol(0).topology().numAtoms();++k) {
        if (name ==  d_sys -> sol(0).topology().atom(k).name()) {
	  solvatom = k;
	}
       }

       std::cout << n  
	       << setw(8) << d_donors.mol(i)  
               << setw(4) << d_sys -> mol(d_donors.mol(i)).topology().resNum(d_donors.atom(i))+1
		          << d_sys -> mol(d_donors.mol(i)).topology().resName(d_sys -> mol(d_donors.mol(i)).topology().resNum(d_donors.atom(i)))
	       << setw(2) << "-" 
               << setw(4) << c
		 << setw(4) //<< d_sys -> sol(0).topology().resNum(solvatom)+1
	 //<< d_sys -> sol(0).topology().resName(d_sys -> sol(0).topology().resNum(solvatom))
               << setw(6) << d_donors_bound_to.atom(i)+1 
               << setw(4) << d_sys -> mol(d_donors_bound_to.mol(i)).topology().atom(d_donors_bound_to.atom(i)).name() 
               << setw(2) << "-"
               << setw(6) << d_donors.atom(i)+1  
               << setw(4) << d_sys -> mol(d_donors.mol(i)).topology().atom(d_donors.atom(i)).name() 
               << setw(2) << "-"
               << setw(6) << d_acceptors_solv.atom(j)+1 
               << setw(4) << d_sys -> sol(0).topology().atom(solvatom).name();
     std::cout.precision(3); 
     std::cout << setw(8) << dist[c]/(double) num[c];///d_hbond.meandist();
     std::cout.precision(3);
     std::cout << setw(8) << ang[c]/(double) num[c];//d_hbond.meanangle();
     std::cout.precision(0);
     std::cout << setw(8) << num[c];
     std::cout.setf(ios::floatfield, ios_base::fixed);
     std::cout.precision(2);             
     std::cout << setw(8) << ((num[c]/ (double) d_frames)*100)
                          << endl;
     
      }
     }
    }


       int cc = 0;
    for (int i=0, m=0; i < (d_donors_solv).size(); ++i, ++m) {
     for (int j=0; j < (d_acceptors).size(); ++j) {
       ++c; ++cc;
 
      if (num[c] > 0) {
       ++n;

       int solvatom = 0;
       string name = (d_donors_solv).name(i);
       for (int k=0;k< d_sys -> sol(0).topology().numAtoms();++k) {
        if (name ==  d_sys -> sol(0).topology().atom(k).name()) {
	  solvatom = k;
	}
       }

       int solvatomb = 0;
       name = (d_donors_bound_to_solv).name(i);
       for (int k=0;k< d_sys -> sol(0).topology().numAtoms();++k) {
        if (name ==  d_sys -> sol(0).topology().atom(k).name()) {
	  solvatomb = k;
	}
       }

       std::cout << n  
	       << setw(8) << cc
		 << setw(4) //<< d_sys -> sol(0).topology().resNum(solvatom)+1
	 //<< d_sys -> sol(0).topology().resName(d_sys -> sol(0).topology().resNum(solvatom))
	       << setw(2) << "-" 
               << setw(4) << d_acceptors.mol(j)
	       << setw(4) << d_sys -> mol(d_acceptors.mol(j)).topology().resNum(d_acceptors.atom(j))+1
	                  << d_sys -> mol(d_acceptors.mol(j)).topology().resName(d_sys -> mol(d_acceptors.mol(j)).topology().resNum(d_acceptors.atom(j)))
               << setw(6) << d_donors_bound_to_solv.atom(i)+1 
               << setw(4) << d_sys -> sol(0).topology().atom(solvatomb).name() 
               << setw(2) << "-"
               << setw(6) << d_donors_solv.atom(i)+1  
               << setw(4) << d_sys -> sol(0).topology().atom(solvatom).name() 
               << setw(2) << "-"
               << setw(6) << d_acceptors.atom(j)+1 
               << setw(4) << d_sys -> mol(d_acceptors.mol(j)).topology().atom(d_acceptors.atom(j)).name();
     std::cout.precision(3); 
     std::cout << setw(8) << dist[c]/(double) num[c];///d_hbond.meandist();
     std::cout.precision(3);
     std::cout << setw(8) << ang[c]/(double) num[c];//d_hbond.meanangle();
     std::cout.precision(0);
     std::cout << setw(8) << num[c];
     std::cout.setf(ios::floatfield, ios_base::fixed);
     std::cout.precision(2);             
     std::cout << setw(8) << ((num[c]/ (double) d_frames)*100)
                          << endl;
     
      }
     }
    }


} //end print

 void Hbondcalc::setmaxdist(double i)
    {
      d_maxdist = i;
    }

 void Hbondcalc::setminangle(double i)
 {
  d_minangle = i;
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
  timeseriesHBtot << setw(10) << d_time;
  timeseriesHBtot.precision(5);
  timeseriesHBtot << setw(10) << d_numHB << endl;

}


Hbondcalc::Hbondcalc(gcore::System &sys, args::Arguments &args)
{
  d_sys=&sys;
  d_args=&args;
  d_donors = AtomSpecifier(sys);
  d_donors_bound_to = AtomSpecifier(sys);
  d_acceptors = AtomSpecifier(sys);
  d_donors_solv = AtomSpecifier(sys);
  d_donors_bound_to_solv = AtomSpecifier(sys); 
  d_acceptors_solv = AtomSpecifier(sys);

  //open timeseries file
  opents("Hbts.out", "Hbnumts.out");

}
