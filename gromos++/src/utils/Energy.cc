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
  d_nkt=0.0;
  d_eps=1.0;
  d_kap=0.0;
  d_cut=1.4;
}
void Energy::calc()
{
  // first, get the non-bonded interactions
  
  // set some properties to zero
  d_p_vdw = 0.0;
  d_p_el = 0.0;

  // define some variables that we will need
  int ch_b=0, ch_e=0, sf=0, pair=0, s=0;
  gmath::Vec d, dd;
  double qq, alc=0, cuts, d1, d2, d6, drf, c6=0, c12=0, vdw, el;
  double cut3= d_cut*d_cut*d_cut;
  double al2 = d_lam*d_lam*d_alj;
  double crf = ((2-2*d_eps)*(1+d_kap*d_cut)-d_eps*(d_kap*d_kap*d_cut*d_cut)) /
  	       ((1+2*d_eps)*(1+d_kap*d_cut)+d_eps*(d_kap*d_kap*d_cut*d_cut));
  int na=d_sys->sol(0).topology().numAtoms();
  int tna=d_sys->sol(0).numCoords();
  int at=0;
  
  // loop over the atoms
  for(int i=0;  i<d_as->size(); i++){
    int mi=d_as->mol(i);
    int ai=d_as->atom(i);
    int sft=0;
    gmath::Vec vi=d_sys->mol(mi).pos(ai);
    // calculate the coordinates of the charge-group we are in
    gmath::Vec chgrp1=calcChgrp(i);
    // check if this atom is soft
    for(s=0; s< d_soft->size();s++)
      if(mi==d_soft->mol(s)&&ai==d_soft->atom(s)) sft=1;
    // set the arrays for this atom to zero
    d_vdw_m[i]=0.0; d_el_m[i]=0.0;
    d_vdw_s[i]=0.0; d_el_s[i]=0.0;

    // first, we do the interactions with solute
    // loop over the molecules
    for(int m=0; m<d_sys->numMolecules(); m++){
      // loop over the charge groups  
      ch_b=0; ch_e=0;
      while(ch_b<d_sys->mol(m).numAtoms()){

	gmath::Vec chgrp2(0.0,0.0,0.0);
        for(ch_e=ch_b;
	    d_sys->mol(m).topology().atom(ch_e).chargeGroup()!=1;
            ch_e++)
	  chgrp2+=d_sys->mol(m).pos(ch_e);
        chgrp2+=d_sys->mol(m).pos(ch_e);
	chgrp2/=(ch_e-ch_b+1);
	// calculate the distance chgrp1 - chgrp2
        chgrp2=d_pbc->nearestImage(chgrp1, chgrp2, d_sys->box());
        d=chgrp2-chgrp1;
 
        if(d.abs2()<=d_cut*d_cut){

	  // charge group is within cut-off, loop over atoms
	  for(int a=ch_b; a<=ch_e; a++){
	    // check if excluded
	    if(m!=mi||!d_ex[i].count(a)){
	      // determine parameters
              gcore::LJType lj(d_gff->ljType(AtomPair(
		  d_sys->mol(m).topology().atom(a).iac(),
		  d_sys->mol(mi).topology().atom(ai).iac())));
              qq=d_sys->mol(m).topology().atom(a).charge() *
		  d_sys->mol(mi).topology().atom(ai).charge();
	      // check third neighbour
              if(d_third[i].count(a)&&m==mi){
		c6=lj.cs6(); c12=lj.cs12();
	      }
	      else {
		c6=lj.c6(); c12=lj.c12();
	      }
	      // now, we calculate the distance between atoms
              dd=d_pbc->nearestImage(vi, d_sys->mol(m).pos(a),
				     d_sys->box());
              d1=(vi-dd).abs();
	      d2=d1*d1;
	      d6=d2*d2*d2;
 
	      // check if we have a soft atom
	      sf=0;
              if(!sft)
	        for(s=0; s<d_soft->size()&&!sf;s++)
		  if(a==d_soft->atom(s)&&m==d_soft->mol(s)) sf=1;
              if(sft||sf){
		if(c6!=0&&c12!=0) d6+=al2*c12/c6;
		
                alc=d_lam*qq*d_gff->fpepsi()/d_nkt;
		alc=alc*alc;
                cuts=alc+d_cut*d_cut;
		cuts=cuts*sqrt(cuts);
		drf = 1/sqrt(alc+d2) - 0.5*crf*d2/cuts - (1-0.5*crf)/d_cut;
	      }
	      else 
		drf=1/d1 - 0.5*crf*d2/cut3 - (1-0.5*crf)/d_cut; 
              vdw=(c12/d6-c6)/d6;
              el=qq*drf*d_gff->fpepsi();
              // finally, check if atom a was also in d_as
              pair=0;
              for(s=0;s<d_as->size()&&!pair;s++)
		if(a==d_as->atom(s)&&m==d_as->mol(s)){
                  pair=1;
                  d_p_vdw+=0.5*vdw;
                  d_p_el+=0.5*el;
		}
	      // and store the energy in the correct array
              d_vdw_m[i]+=vdw;
              d_el_m[i]+=el;
	    }
	  }
        }
	ch_b=ch_e+1;
      }
    }
    // Now, we have to do the same for the solvent. In GROMOS the solvent
    // charge group is centered on the first atom.
    for(int a=0; a<tna; a+=na){
      d=d_pbc->nearestImage(chgrp1, d_sys->sol(0).pos(a), d_sys->box())
            - chgrp1;

      //check distance
      if(d.abs2()<=d_cut*d_cut){

	for(at=0; at<na;at++){
	  gcore::LJType lj(d_gff->ljType(gcore::AtomPair(
	      d_sys->mol(mi).topology().atom(ai).iac(),
              d_sys->sol(0).topology().atom(at).iac())));
	  qq=d_sys->mol(mi).topology().atom(ai).charge() *
	      d_sys->sol(0).topology().atom(at).charge();
          dd=d_pbc->nearestImage(vi, d_sys->sol(0).pos(a+at),
				 d_sys->box());
	  d1=(vi-dd).abs();
          d2=d1*d1;
	  d6=d2*d2*d2;
          c6=lj.c6();
          c12=lj.c12();
          if(sft){
	    if(c6!=0&&c12!=0) d6+=al2*c12/c6;
            alc=d_lam*qq*d_gff->fpepsi()/d_nkt;
	    alc=alc*alc;
            cuts=alc+d_cut*d_cut;
	    cuts=cuts*sqrt(cuts);
	    drf = 1/sqrt(alc+d2) - 0.5*crf*d2/cuts - (1-0.5*crf)/d_cut;
	  }
	  else 
	    drf=1/d1 - 0.5*crf*d2/cut3 - (1-0.5*crf)/d_cut; 

	  d_vdw_s[i]+=(c12/d6 - c6)/d6;
          d_el_s[i]+=qq*drf*d_gff->fpepsi();
        }
      }
    }
  }  

  // next, do the covalent interactions
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
  double qq, alc=0, cuts, d, d1, d2, d6, drf, c6=0, c12=0;
  double cut3=d_cut*d_cut*d_cut;
  double al2=d_lam*d_lam*d_alj;
  double crf = ((2-2*d_eps)*(1+d_kap*d_cut)-d_eps*(d_kap*d_kap*d_cut*d_cut)) /
               ((1+2*d_eps)*(1+d_kap*d_cut)+d_eps*(d_kap*d_kap*d_cut*d_cut));
  int ai=d_as->atom(i);
  int aj=d_as->atom(j);
  int mi=d_as->mol(i);
  int mj=d_as->mol(j);
  int soft=0;
  gmath::Vec dd;
  gmath::Vec chgrp1=calcChgrp(i);
  gmath::Vec chgrp2=calcChgrp(j);
  // check if one of the atoms is soft
  for(int s=0; s< d_soft->size(); s++)
    if((mi==d_soft->mol(s) && ai==d_soft->atom(s)) ||
       (mj==d_soft->mol(s) && aj==d_soft->atom(s)))
      soft=1;
  // calculate the distances between the chargegroups
  chgrp2=d_pbc->nearestImage(chgrp1, chgrp2, d_sys->box());
  d=(chgrp2-chgrp1).abs2();
  if(d<=d_cut*d_cut){
    if(mi!=mj||!d_ex[i].count(aj)){
      //determine parameters
      gcore::LJType lj(d_gff->ljType(AtomPair(
            d_sys->mol(mj).topology().atom(aj).iac(),
            d_sys->mol(mi).topology().atom(ai).iac())));
      qq=d_sys->mol(mj).topology().atom(aj).charge() *
         d_sys->mol(mi).topology().atom(ai).charge();
      // check third neighbour
      if(d_third[i].count(aj)&&mj==mi){
        c6=lj.cs6(); c12=lj.cs12();
      }
      else {
        c6=lj.c6(); c12=lj.c12();
      }
      // now, we calculate the distance between atoms
      dd=d_pbc->nearestImage(d_sys->mol(mi).pos(ai),
			     d_sys->mol(mj).pos(aj),
			     d_sys->box());
      d1=(d_sys->mol(mi).pos(ai)-dd).abs();
      d2=d1*d1;
      d6=d2*d2*d2;
      if(soft){
        if(c6!=0&&c12!=0) d6+=al2*c12/c6;
        alc=d_lam*qq*d_gff->fpepsi()/d_nkt;
	alc=alc*alc;
        cuts=alc+d_cut*d_cut;
	cuts=cuts*sqrt(cuts);
	drf = 1/sqrt(alc+d2) - 0.5*crf*d2/cuts - (1-0.5*crf)/d_cut;
      }
      else 
	drf=1/d1 - 0.5*crf*d2/cut3 - (1-0.5*crf)/d_cut; 

      vdw=(c12/d6-c6)/d6;
      el=qq*drf*d_gff->fpepsi();
    }
  }  
}

int Energy::setAtoms(utils::AtomSpecifier &as)
{
  d_as=&as;
  // for all specified atoms, determine all excluded atoms and all third 
  // neighbours

  for(int i=0; i<d_as->size();i++){
    std::set<int> ex, third;
    int m=d_as->mol(i);
    int a=d_as->atom(i);
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

