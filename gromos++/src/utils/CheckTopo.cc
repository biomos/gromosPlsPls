#include <iostream>
#include <string>
#include <sstream>
#include <set>
#include "CheckTopo.h"
#include "../gcore/System.h"
#include "../gcore/Molecule.h"
#include "../gcore/MoleculeTopology.h"
#include "../gcore/AtomTopology.h"
#include "../gcore/Bond.h"
#include "../gcore/Angle.h"
#include "../gcore/Improper.h"
#include "../gcore/Dihedral.h"
#include "../gcore/Exclusion.h"
#include "../utils/Neighbours.h"
#include "../gromos/Exception.h"

using namespace gcore;
using namespace std;
using namespace utils;

using utils::CheckTopo;
CheckTopo::CheckTopo(const gcore::System &sys, int m)
{
  if(m<sys.numMolecules()){
    d_mt = &sys.mol(m).topology();
    d_chargePrecision = 5;
  }
  else
    throw gromos::Exception("CheckTopo", "Too few molecules in topology");
}

int CheckTopo::checkBonds()
{
  int num_errors_before=d_error.size();

  // loop over the atoms
  for(int a=0; a< d_mt->numAtoms(); a++){
    Neighbours neigh(*d_mt, a);
      
    //check for multiple bonds make sets for the exclusions
    for(unsigned int i=0; i< neigh.size(); i++){
      //check for multiple bonds
      if(neigh[i]>a){
	for(unsigned int j=i+1; j< neigh.size(); j++){
	  if(neigh[i]==neigh[j]){
	    ostringstream os;
	    os << "More than one bond connecting atoms " << a+1 << " and " 
	       << neigh[i]+1;
	    d_error.push_back(os.str());
	    
	  }
	}
      }
    }
  }
  return d_error.size() - num_errors_before;
}

int CheckTopo::checkAngles()
{
  int num_errors_before=d_error.size();

  // loop over the atoms
  for(int a=0; a< d_mt->numAtoms(); a++){
    Neighbours neigh(*d_mt, a);
    
    if(neigh.size()>1){
      for(unsigned int i=0; i< neigh.size(); i++){
	for(unsigned int j=i+1; j< neigh.size(); j++){
	  int b,c;
	  if (neigh[i]<neigh[j]){b=neigh[i]; c=neigh[j];}
	  else                  {b=neigh[j]; c=neigh[i];}
	  AngleIterator ai(*d_mt);
	  int count=0;
	  for(;ai;++ai)
	    if(ai()[0]==b&&ai()[1]==a&&ai()[2]==c) count++;
	  if(count==0){
	    ostringstream os;
	    os << "No angle in topology for atoms " 
	       << b+1 << "-" << a+1 << "-" << c+1;
	    d_error.push_back(os.str());
	  }
	  if(count>1){
	    ostringstream os;
	    os << "More than one angle in topology "
	       << "for atoms " << b+1 << "-" << a+1 << "-" << c+1;
	    d_error.push_back(os.str());
	  }
	}
      }
    }
  }
  return d_error.size() - num_errors_before;
}

int CheckTopo::checkImpropers()
{
  int num_errors_before=d_error.size();

  // loop over the atoms
  for(int a=0; a< d_mt->numAtoms(); a++){
    Neighbours neigh(*d_mt, a);
    if(neigh.size()==3){
      int at[4]={a,neigh[0], neigh[1], neigh[2]};
      ImproperIterator ii(*d_mt);
      int count=0;
      for(;ii;++ii){
	int found[4]={0,0,0,0};
	for(int i=0;i<4;i++)
	  for(int j=0; j<4; j++)
	    if(ii()[i]==at[j]) found[j]=1;
	if(found[0]&&found[1]&&found[2]&&found[3]) count++;
      }
      if(count==0){
	ostringstream os;
	os << "No improper dihedral in topology " 
	   << "for atoms " << at[0]+1 << "-" << at[1]+1 << "-" << at[2]+1 
	   << "-" << at[3]+1;
	d_error.push_back(os.str());
	}
      if(count>1){
	ostringstream os;
	os << "More than one dihedral in topology "
	   << "for atoms " << at[0]+1 << "-" << at[1]+1 << "-" << at[2]+1
	   << "-" << at[3]+1;
	d_error.push_back(os.str());
      }
    }
  }
  return d_error.size() - num_errors_before;
}

int CheckTopo::checkChargeGroups()
{
  int num_errors_before=d_error.size();
  double chrg_precision=1.0;
  for(int i=0; i<d_chargePrecision; i++) chrg_precision/=10.0;
  
  double charge=0.0;
  int chargerest, chargegroup=0;
  
  // loop over the atoms
  for(int a=0; a< d_mt->numAtoms(); a++){ 
    charge+=d_mt->atom(a).charge();
    if(d_mt->atom(a).chargeGroup()){
      chargegroup++;
      chargerest=int(charge/chrg_precision)%int(1.0/chrg_precision);
      if(chargerest){
	ostringstream os;
	os << "Non-integer valued charge in charge "
	   << "group " << chargegroup << ".\n"
	   << "Ends at atom " << a+1 << " : " 
	   << chrg_precision*chargerest
	   << endl;
	d_error.push_back(os.str());
      }
      charge=0.0;
    }
  }
  return d_error.size() - num_errors_before;
}

void CheckTopo::setChargePrecision(int i)
{
  d_chargePrecision = i;
}

int CheckTopo::checkExclusions()
{
  int num_errors_before=d_error.size();

  // a set of residue numbers that should be aromatic
  set<int> aromatic;
  // loop over the atoms
  for(int a=0; a< d_mt->numAtoms(); a++){
    Neighbours n12(*d_mt, a);

    // ex13 will contain all atoms that are 1,2 or 1,3 neighbours. Atoms
    // that are in ex14 BUT NOT IN ex13 are 1,4 neighbours.
    set<int> ex13;
    set<int> ex14;
    for(unsigned int i=0; i< n12.size(); i++){
      ex13.insert(n12[i]);
      Neighbours n13(*d_mt, n12[i]);
      for(unsigned int j=0; j<n13.size(); j++){
	ex13.insert(n13[j]);
	//and their neighbours
	Neighbours n14(*d_mt, n13[j]);
	for(unsigned int k=0; k< n14.size(); k++){
	  ex14.insert(n14[k]);	
	}
      }
    }
    // check for the standard exclusion rules
    // ignore14 will contain all 1,4 neighbours that are in the exlcusion list
    // we will not complain about those missing in the 1,4 exclusions.
    set<int> ignore14;
    for(int i=0; i< d_mt->atom(a).exclusion().size(); i++){
      int j=d_mt->atom(a).exclusion().atom(i);
      if(j>a && ex13.count(j)==0){
	if(ex14.count(j)){
	  ignore14.insert(j);
	  aromatic.insert(d_mt->resNum(a));
	}
	else{
	  ostringstream os;
	  os << "Atom " << j+1 << " is in the exclusion list of atom " 
	     << a+1 << " but is not a 1,2 or 1,3 neighbour";
	  d_error.push_back(os.str());
	}
      }
    }
    for(set<int>::const_iterator iter=ex13.begin(), to=ex13.end(); 
	iter!=to; ++iter){
      if(*iter>a){
	int found =0;
	for(int i=0; i<d_mt->atom(a).exclusion().size();i++)
	  if(*iter==d_mt->atom(a).exclusion().atom(i)) 
	    found=1;
	if(!found){
	  ostringstream os;
	  os << "Atom " << *iter + 1 << " is a 1,2 or 1,3 neighbour of atom "
	     << a+1 << " but is not in its exclusion list" ;
	  d_error.push_back(os.str());
	}
      }
    }
    // check for 1,4 exclusions
    for(int i=0; i< d_mt->atom(a).exclusion14().size();i++){
      int j=d_mt->atom(a).exclusion14().atom(i);
      if(j>a && ex14.count(j)==0){
	ostringstream os;
	os << "Atom " << j+1 << " is in the 1,4 exclusion list of atom " 
	   << a+1 << " but is not a 1,4 neighbour";
	d_error.push_back(os.str());
      }
    }
    for(set<int>::const_iterator iter=ex14.begin(), to=ex14.end(); 
	iter!=to; ++iter){
      if(*iter>a && ex13.count(*iter)==0){
	int found =0;
	for(int i=0; i<d_mt->atom(a).exclusion14().size(); i++)
	  if(*iter==d_mt->atom(a).exclusion14().atom(i)) 
	    found=1;
	if(!found && ignore14.count(*iter)==0){
	  ostringstream os;
	  os << "Atom " << *iter + 1 << " is a 1,4 neighbour of atom "
	     << a+1 << " but is not in its 1,4 exclusion list";
	  d_error.push_back(os.str());
	}
      }
    } 
  }
  
  if(aromatic.size()){
    ostringstream os;
    os << "Residue";
    if(aromatic.size()>1) os << "s";
    for(set<int>::const_iterator iter=aromatic.begin(),
	  to=aromatic.end(); iter!=to; ++iter){
      os << " " << (*iter) + 1 << d_mt->resName(*iter);
    }
    if(aromatic.size()>1) os << " have ";
    else os << " has ";
    os << "1,4 neighbours in the exclusion lists.\n";
    os << "If no aromatic systems are involved"
       << " there might be a problem.";
    d_error.push_back(os.str());
  }
  return d_error.size() - num_errors_before;
}

int CheckTopo::checkAll()
{
  int num_errors_before=d_error.size();
  checkBonds();
  checkAngles();
  checkImpropers();
  checkChargeGroups();
  checkExclusions();
  return d_error.size() - num_errors_before;
}


void CheckTopo::clearErrors()
{
  d_error.resize(0);
}
    
int CheckTopo::chargePrecision()
{
  return d_chargePrecision;
}

int CheckTopo::numErrors()
{
  return d_error.size();
}

std::string CheckTopo::error(int i)
{
  if(i < int(d_error.size()))
    return d_error[i];
  else
    throw gromos::Exception("CheckTopo", 
			    "Trying to access an error that does noet exist");
  
}
