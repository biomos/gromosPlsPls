#include "Energy.h"
#include "../gcore/System.h"
#include "../gcore/Molecule.h"
#include "../gcore/MoleculeTopology.h"
#include "../gcore/Solvent.h"
#include "../gcore/SolventTopology.h"
#include "../gcore/AtomTopology.h"
#include "../gcore/Exclusion.h"
#include "../gcore/AtomPair.h"
#include "../gcore/LJType.h"
#include "../gcore/Bond.h"
#include "../gcore/BondType.h"
#include "../gcore/Angle.h"
#include "../gcore/AngleType.h"
#include "../gcore/Improper.h"
#include "../gcore/ImproperType.h"
#include "../gcore/Dihedral.h"
#include "../gcore/DihedralType.h"
#include "../gcore/GromosForceField.h"
#include "../bound/Boundary.h"
#include "AtomSpecifier.h"
#include "SimplePairlist.h"
#include "PropertyContainer.h"
#include "Property.h"
#include "../gmath/Vec.h"

using namespace gcore;
using namespace std;
using namespace utils;
//using utils::Energy;
namespace utils
{
Energy::Energy(gcore::System &sys, gcore::GromosForceField &gff,
               bound::Boundary &pbc)
{
  d_sys=&sys;
  d_gff=&gff;
  d_pbc=&pbc;
  d_as = new utils::AtomSpecifier(sys);
  d_pc = new utils::PropertyContainer(sys);
  d_soft= new utils::AtomSpecifier(sys);
  d_lam=0.0;
  d_alj=0.0;
  d_ac=0.0;
  d_eps=1.0;
  d_kap=0.0;
  d_cut=1.4;
}
void Energy::calc()
{
  calcNb();
  calcCov();
}
void Energy::calcNb()
{
  // set some properties to zero
  d_p_vdw = 0.0; d_p_el = 0.0;

  // define some variables that we will need
  int sf=0;
  gmath::Vec d, dd;
  double qq, cuts, d1, d2, d6, drf, c6=0, c12=0, vdw, el;
  double cut3= d_cut*d_cut*d_cut;
  double l2alj = d_lam*d_lam*d_alj;
  double l2ac = d_lam*d_lam*d_ac;
  double crf = ((2-2*d_eps)*(1+d_kap*d_cut)-d_eps*(d_kap*d_kap*d_cut*d_cut)) /
  	       ((1+2*d_eps)*(1+d_kap*d_cut)+d_eps*(d_kap*d_kap*d_cut*d_cut));
  double dirf = (1 - 0.5 * crf) / d_cut;

  // create a pairlist
  SimplePairlist pl(*d_sys, *d_pbc, d_cut);
  
  // loop over the atoms
  for(int i=0;  i<d_as->size(); i++){
    int mi=d_as->mol(i);
    int ai=d_as->atom(i);
    int iaci=d_as->iac(i);
    double qi=d_as->charge(i);
    int sft=0;
    gmath::Vec vi=d_sys->mol(mi).pos(ai);
    // calculate the coordinates of the charge-group we are in
    gmath::Vec chgrp1=calcChgrp(i);
    // check if this atom is soft
    if(d_soft->findAtom(mi, ai)!=-1) sft=1;
    // set the arrays for this atom to zero
    d_vdw_m[i]=0.0; d_el_m[i]=0.0;
    d_vdw_s[i]=0.0; d_el_s[i]=0.0;

    // set and calculate the pairlist
    pl.clear();
    pl.setAtom(mi,ai);
    pl.calcCgb();
    pl.removeExclusions();
    //cout << pl.size() << " elements in the pairlist" << endl;
    
    
    // now, loop over the pairlist
    for(int j=0; j<pl.size(); j++){
      int mj=pl.mol(j);
      int aj=pl.atom(j);

      // determine parameters
      gcore::LJType lj(d_gff->ljType(AtomPair(iaci, pl.iac(j))));
      if(d_third[i].count(aj)&&mj==mi){
	c6=lj.cs6(); c12=lj.cs12();
      }
      else {
	c6=lj.c6(); c12=lj.c12();
      }
      qq = qi * pl.charge(j);
      // now, we calculate the distance between atoms
      dd=d_pbc->nearestImage(vi, *pl.coord(j), d_sys->box());

      d1=(vi-dd).abs();
      d2=d1*d1;
      d6=d2*d2*d2;
 
      // check if we have a soft atom
      sf=0;
      if(sft || d_soft->findAtom(mj,aj)!=-1) sf=1;
      if(sf){
	if(c6!=0&&c12!=0) d6+=l2alj*c12/c6;
	cuts=l2ac+d_cut*d_cut;
	cuts=cuts*sqrt(cuts);
	drf = 1/sqrt(l2ac+d2) - 0.5*crf*d2/cuts - dirf;
      }
      else 
	drf=1/d1 - 0.5*crf*d2/cut3 - dirf;
 
      vdw=(c12/d6-c6)/d6;
      el=qq*drf*d_gff->fpepsi();

      // finally, check if atom a was also in d_as
      if(d_as->findAtom(mj,aj)!=-1){
	d_p_vdw += 0.5 * vdw;
	d_p_el  += 0.5 * el;
      }
      // and store the energy in the correct array
      if(mj<0){
	d_vdw_s[i] += vdw; d_el_s[i] += el;
      }
      else{
	d_vdw_m[i] += vdw; d_el_m[i] += el;
      }
    }
  }
}
void Energy::calcCov()
{
  // first calculate the values for all properties
  d_pc->calc();

  // loop over properties
  for(unsigned int i=0;i<d_pc->size();i++){
    switch((*d_pc)[i]->atoms().size()){
      case 2:
	d_cov[i]=calcBond((*d_pc)[i]->getValue(), d_covpar[i]);
        break;
      case 3:
	d_cov[i]=calcAngle((*d_pc)[i]->getValue(), d_covpar[i]);
	break;
    case 4:
      	if(d_covpar[i][0] == -1000.0)
	  d_cov[i]=calcImproper((*d_pc)[i]->getValue(), d_covpar[i]);
	else
	  d_cov[i]=calcDihedral((*d_pc)[i]->getValue(), d_covpar[i]);
	break;
    }
  }
}
void Energy::calcPair(int i, int j, double &vdw, double &el)
{
  double qq, cuts, d, d1, d2, d6, drf, c6=0, c12=0;
  double cut3=d_cut*d_cut*d_cut;
  double l2alj=d_lam*d_lam*d_alj;
  double l2ac=d_lam*d_lam*d_ac;
  double crf = ((2-2*d_eps)*(1+d_kap*d_cut)-d_eps*(d_kap*d_kap*d_cut*d_cut)) /
               ((1+2*d_eps)*(1+d_kap*d_cut)+d_eps*(d_kap*d_kap*d_cut*d_cut));
  double dirf = (1.0-0.5*crf)/d_cut;
  
  int ai=d_as->atom(i);
  int aj=d_as->atom(j);
  int mi=d_as->mol(i);
  int mj=d_as->mol(j);
  int soft=0;
  gmath::Vec dd;
  gmath::Vec chgrp1=calcChgrp(i);
  gmath::Vec chgrp2=calcChgrp(j);
  // check if one of the atoms is soft
  if(d_soft->findAtom(mi,ai)!=-1 ||
     d_soft->findAtom(mj,aj)!=-1) soft=1;

  // calculate the distances between the chargegroups
  chgrp2=d_pbc->nearestImage(chgrp1, chgrp2, d_sys->box());
  d=(chgrp2-chgrp1).abs2();
  if(d<=d_cut*d_cut){
    if(mi!=mj||!d_ex[i].count(aj)){
      //determine parameters
      gcore::LJType lj(d_gff->ljType(AtomPair(d_as->iac(i), d_as->iac(j))));
      qq = d_as->charge(i) * d_as->charge(j);

      // check third neighbour
      if(d_third[i].count(aj)&&mj==mi){
        c6=lj.cs6(); c12=lj.cs12();
      }
      else {
        c6=lj.c6(); c12=lj.c12();
      }
      // now, we calculate the distance between atoms
      dd=d_pbc->nearestImage(*(d_as->coord(i)),
			     *(d_as->coord(j)),
			     d_sys->box());
      d1=(*d_as->coord(i)-dd).abs();
      d2=d1*d1;
      d6=d2*d2*d2;
      if(soft){
        if(c6!=0&&c12!=0) d6+=l2alj*c12/c6;
        cuts=l2ac+d_cut*d_cut;
	cuts=cuts*sqrt(cuts);
	drf = 1/sqrt(l2ac+d2) - 0.5*crf*d2/cuts - dirf;
      }
      else 
	drf=1/d1 - 0.5*crf*d2/cut3 - dirf;

      vdw=(c12/d6-c6)/d6;
      el=qq*drf*d_gff->fpepsi();
    }
  }  
}

int Energy::setAtoms(utils::AtomSpecifier &as)
{    
  d_ex.resize(0);
  d_third.resize(0);
  d_vdw_m.resize(0);
  d_vdw_s.resize(0);
  d_el_m.resize(0);
  d_el_s.resize(0);
  
  d_as=&as;
  // for all specified atoms, determine all excluded atoms and all third 
  // neighbours

  for(int i=0; i<d_as->size();i++){
    std::set<int> ex, third;
    int m=d_as->mol(i);
    int a=d_as->atom(i);
    if(m>=0){
      
      // first find atoms (<a) from which a is excluded
      for(int ai=0; ai<a; ai++)
	for(int e=0; e<d_sys->mol(m).topology().atom(ai).exclusion().size(); e++)
	  if(a==d_sys->mol(m).topology().atom(ai).exclusion().atom(e))
	    ex.insert(ai);
      // now, we add the exclusions of a itself
      for(int e=0; e<d_sys->mol(m).topology().atom(a).exclusion().size(); e++)
	ex.insert(d_sys->mol(m).topology().atom(a).exclusion().atom(e));
      // and a is excluded of itself
      ex.insert(a);
      
      // first find atoms (<a) which have a as third neighbour
      for(int ai=0; ai<a; ai++)
	for(int e=0; e<d_sys->mol(m).topology().atom(ai).exclusion14().size();
	    e++)
	  if(a==d_sys->mol(m).topology().atom(ai).exclusion14().atom(e))
	    third.insert(ai);
      // now, we add the third neighbours of a itself
      for(int e=0; e<d_sys->mol(m).topology().atom(a).exclusion14().size(); e++)
	third.insert(d_sys->mol(m).topology().atom(a).exclusion14().atom(e));
    }
    
    // add things to the necessary vectors
    d_ex.push_back(ex);
    d_third.push_back(third);
	       
    d_vdw_m.push_back(0.0);
    d_vdw_s.push_back(0.0);
    d_el_m.push_back(0.0);
    d_el_s.push_back(0.0);
  }
  return d_as->size();
}

int Energy::setProperties(utils::PropertyContainer &pc)
{
  d_pc=&pc;
  for(unsigned int i=0;i<d_pc->size();i++){
    gmath::Vec temp;
    switch(pc[i]->atoms().size()){
      case 2:{
	int t=findBond(*pc[i]);
        temp[0]=d_gff->bondType(t).b0();
        temp[1]=d_gff->bondType(t).fc();
      }
      break;
      case 3:{
	int t=findAngle(*pc[i]);
	temp[0]=d_gff->angleType(t).t0();
        temp[1]=d_gff->angleType(t).fc();
      }
      break;
      case 4:{
	int t=findDihedral(*pc[i]);
	// Dirty fix to deal with both impropers and dihedrals:
        // if t < 0 this means it is an improper, with types counting
	// -1, -2, ...
	// in the function calc, the improper is recognised by
	// the first parameter-value of -1000 (useless in case of torsion)
	if(t>=0){
	  temp[0]=d_gff->dihedralType(t).pd();
	  temp[1]=d_gff->dihedralType(t).np();
          temp[2]=d_gff->dihedralType(t).fc();
	}
	else {
	  t = -1*(t+1);
	  temp[0]=-1000.0;
	  temp[1]=d_gff->improperType(t).q0();
	  temp[2]=d_gff->improperType(t).fc();
	}
      }
      break; 
      default:
        throw Energy::Exception(
	   " number of atoms in this property is unknown: "+ pc[i]->toTitle());    
    }
    d_covpar.push_back(temp);
    d_cov.push_back(0.0);

  }
  return d_pc->size();
}

gmath::Vec Energy::calcChgrp(int i){
  gmath::Vec chgrp(0.0, 0.0, 0.0);
  int mi=d_as->mol(i), ai=d_as->atom(i);
  if(mi<0){
    int nsa=d_sys->sol(0).topology().numAtoms();
    int solv=ai/nsa;
    solv*=nsa;
    return d_sys->sol(0).pos(solv);
  }
  
  int begin=ai-1, end=ai; 
  if(ai>0)
    for(begin=ai-1; 
	begin>=0&&d_sys->mol(mi).topology().atom(begin).chargeGroup()!=1; 
	begin--);
  for(end=ai;
      d_sys->mol(mi).topology().atom(end).chargeGroup()!=1;
      end++);

  // charge group goes from begin+1 to end
  for(int k=begin+1; k<=end; k++)
    chgrp+=d_sys->mol(mi).pos(k);
  return chgrp/(end-begin);
}

int Energy::findBond(utils::Property &pp){
  int m, a, b,f=0;
  if(pp.mols()[0]==pp.mols()[1])
      m=pp.mols()[0];
  else
    throw Energy::Exception(
       " Covalent interactions are always within one molecule: "+pp.toTitle());
  if(pp.atoms()[0]<pp.atoms()[1]){
    a=pp.atoms()[0];
    b=pp.atoms()[1];
  }
  else{
    a=pp.atoms()[1];
    b=pp.atoms()[0];
  }
  BondIterator bi(d_sys->mol(m).topology());
  while(bi&&f==0)
      if(bi()[0]==a && bi()[1]==b) f=1;
    else ++bi;
  if(bi) return bi().type();
  else
    throw Energy::Exception(
	" Bond not found in topology: "+pp.toTitle());
}
int Energy::findAngle(utils::Property &pp){
  int m, a, b,c,f=0;
  if(pp.mols()[0]==pp.mols()[1]&&pp.mols()[0]==pp.mols()[2])
      m=pp.mols()[0];
  else
    throw Energy::Exception(
       " Covalent interactions are always within one molecule: "+pp.toTitle());
  if(pp.atoms()[0]<pp.atoms()[2]){
    a=pp.atoms()[0];
    b=pp.atoms()[1];
    c=pp.atoms()[2];
  }
  else{
    a=pp.atoms()[2];
    b=pp.atoms()[1];
    c=pp.atoms()[0];
  }
  AngleIterator ai(d_sys->mol(m).topology());
  while(ai&&f==0)
      if(ai()[0]==a && ai()[1]==b && ai()[2]==c) f=1;
    else ++ai;
  if(ai) return ai().type();
  else
    throw Energy::Exception(
	" Angle not found in topology: "+pp.toTitle());
}
int Energy::findDihedral(utils::Property &pp){
  int m, a, b, c, d, f=0;
  if(pp.mols()[0]==pp.mols()[1] &&
     pp.mols()[0]==pp.mols()[2] &&
     pp.mols()[0]==pp.mols()[3])
      m=pp.mols()[0];
  else
    throw Energy::Exception(
       " Covalent interactions are always within one molecule: "+pp.toTitle());
  if(pp.atoms()[1]<pp.atoms()[2]){
    a=pp.atoms()[0];
    b=pp.atoms()[1];
    c=pp.atoms()[2];
    d=pp.atoms()[3];
  }
  else{
    a=pp.atoms()[3];
    b=pp.atoms()[2];
    c=pp.atoms()[1];
    d=pp.atoms()[0];
  }
  DihedralIterator di(d_sys->mol(m).topology());
  while(di&&f==0)
    if(di()[0]==a && di()[1]==b && di()[2]==c && di()[3]==d) f=1;
    else ++di;
  if(di) return di().type();
  else {
    //Maybe we have an improper
    ImproperIterator ii(d_sys->mol(m).topology());
    while(ii&&f==0) 
      if(ii()[0]==a && ii()[1]==b && ii()[2]==c && ii()[3]==d) f=1;
      else ++ii;
    if(ii) return -1*(ii().type()+1);
    else
      throw Energy::Exception(
	" (improper) Dihedral not found in topology: "+pp.toTitle());
  }
}
double Energy::calcBond(double val, gmath::Vec par)
{
  double diff=val*val-par[0]*par[0];
  return 0.25*par[1]*diff*diff;
}
double Energy::calcAngle(double val, gmath::Vec par)
{
  val=val*pi/180.0;
  double t0=par[0]*pi/180.0;
  double diff=cos(val) - cos(t0);
  return 0.5*par[1]*diff*diff;
}
double Energy::calcDihedral(double val, gmath::Vec par)
{
  val=val*pi/180.0;
  return par[2]*(1+par[0]*cos(par[1]*val));
}
double Energy::calcImproper(double val, gmath::Vec par)
{
  if(val>180.0) val=val-360;
  double diff=val-par[1];
  return 0.5*par[2]*diff*diff;
}
double Energy::cov()
{
  double e=0.0;
  for(unsigned int i=0;i< d_pc->size();i++)
      e+=this->cov(i);
  return e;
}
double Energy::vdw()
{
  double e = 0.0;
  for(int i=0; i< d_as->size();i++)
      e+=this->vdw(i);
  return e - d_p_vdw;
}
double Energy::el()
{
  double e=0.0;
  for(int i=0; i< d_as->size();i++)
      e+=this->el(i);
  return e - d_p_el;
}
double Energy::vdw(int i)
{
  assert(i < d_as->size());
  return d_vdw_m[i]+d_vdw_s[i];
}
double Energy::el(int i)
{
  assert(i < d_as->size());
  return d_el_m[i] + d_el_s[i];
}
double Energy::cov(int i)
{
  assert(i < int(d_pc->size()));
  return d_cov[i];
}
double Energy::vdw_m(int i)
{
  assert(i<d_as->size());
  return d_vdw_m[i];
}
double Energy::vdw_s(int i)
{
  assert(i<d_as->size());
  return d_vdw_s[i];
}
double Energy::el_m(int i)
{
  assert(i<d_as->size());
  return d_el_m[i];
}
double Energy::el_s(int i)
{
  assert(i<d_as->size());
  return d_el_s[i];
}
}

