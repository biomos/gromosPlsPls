#include <iostream>
#include <stdio.h>
#include <string>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include "AtomSpecifier.h"
#include "Dssp.h"
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
using utils::Dssp;
using utils::AtomSpecifier;



void Dssp::determineAtoms()
{  
  for(int m=0; m< d_sys->numMolecules(); m++){
    for(int a=1; a< d_sys->mol(m).numAtoms(); a++){
      if(d_sys->mol(m).topology().atom(a-1).name()=="N" &&
	 d_sys->mol(m).topology().atom(a  ).name()=="H"){
	d_H.addAtom(m,a); d_N.addAtom(m,a-1);
      }
      if(d_sys->mol(m).topology().atom(a-1).name()=="C" &&
	 d_sys->mol(m).topology().atom(a  ).name()=="O"){
	d_O.addAtom(m,a); d_C.addAtom(m,a-1);
      }
    }
  } 
} //end Dssp::determineAtoms()


void Dssp::calcHintra_init()
{
  d_pbc = BoundaryParser::boundary(*d_sys, *d_args);  
  d_pbc -> gather();
  d_CA.addSpecifier("a:CA");
} //end Dssp::calcHintra_init()

void Dssp::calcHb_Kabsch_Sander()
{
  acc_res.clear();
  don_res.clear();

  double rON=0, rCH=0, rOH=0, rCN=0;
  double q1=0.42, q2=0.20, f=33.2, E=0, cutoff=-0.5;

  for (int i=0; i < (int) d_O.size(); ++i) {
    for (int j=0; j < (int) d_H.size(); ++j) {
      rON = (*d_N.coord(j) - *d_O.coord(i)).abs();
      rCH = (*d_H.coord(j) - *d_C.coord(i)).abs();
      rOH = (*d_H.coord(j) - *d_O.coord(i)).abs();
      rCN = (*d_N.coord(j) - *d_C.coord(i)).abs();
      E = q1*q2*(1/rON + 1/rCH - 1/rOH - 1/rCN)*f;
      if ((E < cutoff) && (abs((d_sys->mol(0).topology().resNum(d_O.atom(i)) - 
      d_sys->mol(0).topology().resNum(d_H.atom(j)))) > 1)){
	acc_res.push_back(d_sys->mol(0).topology().resNum(d_O.atom(i)));
	don_res.push_back(d_sys->mol(0).topology().resNum(d_H.atom(j)));
      }
    }
  }
} //end Dssp::calcHb_Kabsch_Sander()

void Dssp::calc_Helices()
{  
  vector<int> turn3, turn4, turn5;
  vector<int> helix3_tmp, helix4_tmp, helix5_tmp;

  helix3.clear();
  helix4.clear();
  helix5.clear();
  Helix.clear();
  Turn.clear();

  // loop over existing hydrogen bonds to identify elementary H-bond patterns
  // begin with three types of H-bonded turns
  // fill up Turn - may contain duplicates if one residue H-bonds to more than one 
  // other residue

  for (int i=0; i < (int) acc_res.size(); ++i) {     
     if (don_res[i] == acc_res[i] + 3) {
       turn3.push_back(acc_res[i]);
       Turn.push_back(acc_res[i]);
     }
     if (don_res[i] == acc_res[i] + 4) {
       turn4.push_back(acc_res[i]);
       Turn.push_back(acc_res[i]);
     }
     if (don_res[i] == acc_res[i] + 5) {
       turn5.push_back(acc_res[i]);
       Turn.push_back(acc_res[i]);
     }
  }

  // see if turns form helices
  for (int i=0; i < (int) turn3.size(); ++i) {
    if (turn3[i] == (turn3[i+1] - 1)) {
      helix3_tmp.push_back(turn3[i]+1);
      helix3_tmp.push_back(turn3[i]+2);
      helix3_tmp.push_back(turn3[i]+3);
    }
  }
  for (int i=0; i < (int) turn4.size(); ++i) {
    if (turn4[i] == (turn4[i+1] - 1)) {
      helix4_tmp.push_back(turn4[i]+1);
      helix4_tmp.push_back(turn4[i]+2);
      helix4_tmp.push_back(turn4[i]+3);
      helix4_tmp.push_back(turn4[i]+4);
    }
  }
  for (int i=0; i < (int) turn5.size(); ++i) {
    if (turn5[i] == (turn5[i+1] - 1)) {
      helix5_tmp.push_back(turn5[i]+1);
      helix5_tmp.push_back(turn5[i]+2);
      helix5_tmp.push_back(turn5[i]+3);
      helix5_tmp.push_back(turn5[i]+4);
      helix5_tmp.push_back(turn5[i]+5);
    }
  }
  // remove "duplicates", this will also sort them
  // also fill up Helix vector to help filter results
  for (int i = 0; i < d_sys->mol(0).topology().numRes(); ++i) {
    if (helix3_tmp.size() > 0) {
      for (int j=0; j < (int) helix3_tmp.size(); ++j) {
	if (helix3_tmp[j] == i ) {
	  helix3.push_back(helix3_tmp[j]);
	  Helix.push_back(helix3_tmp[j]);
	  break;
	}
      }
    }
    if (helix4_tmp.size() > 0) {
      for (int j=0; j < (int) helix4_tmp.size(); ++j) {
	if (helix4_tmp[j] == i ) {
	  helix4.push_back(helix4_tmp[j]);
	  Helix.push_back(helix4_tmp[j]);
	  break;
	}
      }
    }
    if (helix5_tmp.size() > 0) {
      for (int j=0; j < (int) helix5_tmp.size(); ++j) {
	if (helix5_tmp[j] == i ) {
	  helix5.push_back(helix5_tmp[j]);
	  Helix.push_back(helix5_tmp[j]);
	  break;
	}
      }
    }
  }
} // end Dssp::calc_Helices()

void Dssp::calc_Betas()
{
  vector<int> p_bridge_tmp, ap_bridge_tmp, p_bridge_tmp2, ap_bridge_tmp2;
  vector<int> bridge_tmp, extended_tmp;

  bridge.clear();
  extended.clear();
  Beta.clear();

  // identify single beta bridges (parallel or antiparallel)
  // loop over Hbonds i and j (they are _not_ residue numbers)
  for (int i=0; i < (int) acc_res.size(); ++i) {
    for (int j=0; j < (int) acc_res.size(); ++j) {
      if ((don_res[i] == acc_res[j]) && (acc_res[i] == (don_res[j]-2))) {
	if (abs(acc_res[i]+1 - don_res[i]) > 2) {
	  p_bridge_tmp.push_back(acc_res[i]+1);
	  p_bridge_tmp.push_back(don_res[i]);
	}
      }
      if ((acc_res[i] == (don_res[j]-2)) && (don_res[i] == (acc_res[j]+2))) {
	if (abs(acc_res[i]+1 - don_res[i]-1) > 2) {
	  ap_bridge_tmp.push_back(acc_res[i]+1);
	  ap_bridge_tmp.push_back(don_res[i]-1);
	}
      }
      if ((don_res[i] == acc_res[j]) && (acc_res[i] == don_res[j])) {
	if (abs(don_res[i] - acc_res[i]) > 2) {
	  ap_bridge_tmp.push_back(don_res[i]);
	  ap_bridge_tmp.push_back(acc_res[i]);
	}
      }
    }
  }
  // remove "duplicates" also for the beta-bridges, this will also sort them
  for (int i=0; i < d_sys->mol(0).topology().numRes(); ++i) {
    if (p_bridge_tmp.size() > 0 ) {  
      for (int j=0; j < (int) p_bridge_tmp.size(); ++j) {
	if (p_bridge_tmp[j] == i ) {
	  p_bridge_tmp2.push_back(p_bridge_tmp[j]);
	  break;
	}
      }
    }
    if (ap_bridge_tmp.size() > 0 ) {
      for (int j=0; j < (int) ap_bridge_tmp.size(); ++j) {
	if (ap_bridge_tmp[j] == i ) {
	  ap_bridge_tmp2.push_back(ap_bridge_tmp[j]);
	  break;
	}
      }
    }
  }
  // isolated bridge or extended strand?
  for (int i=0; i < (int) p_bridge_tmp2.size(); ++i) {    
    if (p_bridge_tmp2[i+1] == (p_bridge_tmp2[i] + 1)) {
      extended_tmp.push_back(p_bridge_tmp2[i]);
      extended_tmp.push_back(p_bridge_tmp2[i+1]);
    }
    else if ((p_bridge_tmp2[i] != (p_bridge_tmp2[i-1] +1)) || (i == 0)){
      bridge_tmp.push_back(p_bridge_tmp2[i]);
    }
  }
  for (int i=0; i < (int) ap_bridge_tmp2.size(); ++i) {
    if (ap_bridge_tmp2[i+1] == (ap_bridge_tmp2[i] + 1)) {
      extended_tmp.push_back(ap_bridge_tmp2[i]);
      extended_tmp.push_back(ap_bridge_tmp2[i+1]);
    }
    else if ((ap_bridge_tmp2[i] != (ap_bridge_tmp2[i-1] +1)) || (i == 0)) {
      bridge_tmp.push_back(ap_bridge_tmp2[i]);
    }
  }
  // remove duplicates, fill up Beta vector
  for (int i=0; i < d_sys->mol(0).topology().numRes(); ++i) {
    if (extended_tmp.size() > 0 ) {
      for (int j=0; j < (int) extended_tmp.size(); ++j) {
	if (extended_tmp[j] == i ) {
	  extended.push_back(extended_tmp[j]);
	  Beta.push_back(extended_tmp[j]);
	  break;
	}
      }
    }
    if (extended_tmp.size() > 0 ) {
      for (int j=0; j < (int) bridge_tmp.size(); ++j) {
	if (bridge_tmp[j] == i ) {
	  bridge.push_back(bridge_tmp[j]);
	  Beta.push_back(bridge_tmp[j]);
	  break;
	}
      }
    }
  }
} //end Dssp::calc_Betas()

void Dssp::calc_Bends()
{
  gmath::Vec tmpA, tmpB; 
  double angle = 0.0;
  Bend.clear();
  d_pbc -> gather();

  for (int i = 2; i < (int) d_CA.size()-2; ++i) {
    tmpA = (*d_CA.coord(i-2) - *d_CA.coord(i));
    tmpB = (*d_CA.coord(i+2) - *d_CA.coord(i));
    angle = acos((tmpA.dot(tmpB))/(sqrt(tmpA.dot(tmpA))*
				   (sqrt(tmpB.dot(tmpB)))))*180/3.1416;
    if (angle < 110) {
      Bend.push_back(i);
    }
  }
} //end Dssp::calc_Bends()

void Dssp::filter_SecStruct()
{
  vector<int>::iterator iter;

  helix.clear();
  turn.clear();

  // ad hoc priority rules:
  // Beta > Helix > turn > bend

  // remove duplicates in Turn and Helix ==> turn, helix
  for (int i=0; i < d_sys->mol(0).topology().numRes(); ++i) {
    for (int j=0; j < (int) Turn.size(); ++j) {
      if (Turn[j] == i ) {
	turn.push_back(i);
	break;
      }
    }
    for (int j=0; j < (int) Helix.size(); ++j) {
     if (Helix[j] == i ) {
	helix.push_back(i);
	break;
      }
    } 
  }
  // first remove (many) turn-residues that form a Helix from turn
  // also remove (many) Bend-residues that are part of a Helix from Bend
  for (int i = 0; i < (int) helix.size(); ++i) {
    for (iter=turn.begin(); iter!=turn.end(); ++iter){
      if (*iter == helix[i]) {
	turn.erase(iter);
	--iter;
      }
    }
    for (iter=Bend.begin(); iter!=Bend.end(); ++iter){
      if (*iter == helix[i]) { 
	Bend.erase(iter);
	--iter;
      }
    }
  }
  // there should be no Bend-residues in turn-residues
  for (int i = 0; i < (int) turn.size(); ++i) {
    if (Bend.size() > 0) {
      for (iter=Bend.begin(); iter!=Bend.end(); ++iter){
	if (*iter == turn[i]) { 
	  Bend.erase(iter);
	  --iter;
	}
      }
    }
  }
  // if helix-residue is part of Beta, remove helix-residue 
  for (int i = 0; i < (int) Beta.size(); ++i) {
    if (helix.size() > 0) {
      if (helix3.size() > 0) {
	for (iter=helix3.begin(); iter!=helix3.end(); ++iter){
	  if (*iter == Beta[i]) {
	    helix3.erase(iter);
	    --iter;
	  }
	}
      }
      if (helix4.size() > 0) {
	for (iter=helix4.begin(); iter!=helix4.end(); ++iter){
	  if (*iter == Beta[i]) {
	    helix4.erase(iter);
	    --iter;
	  }
	}
      }
      if (helix5.size() > 0) {
	for (iter=helix5.begin(); iter!=helix5.end(); ++iter){
	  if (*iter == Beta[i]) {
	    helix5.erase(iter);
	    --iter;
	  }
	}
      }
    }
    // if turn-residue is part of Beta, remove Turn-residue 
    if (turn.size() > 0) {
      for (iter=turn.begin(); iter!=turn.end(); ++iter){
	if (*iter == Beta[i]) {
	  turn.erase(iter);
	  --iter;
	}
      }
    }
    // if Bend-residue is part of Beta, remove Bend-residue
    if (Bend.size() > 0) {
      for (iter=Bend.begin(); iter!=Bend.end(); ++iter){
	if (*iter == Beta[i]) {
	  Bend.erase(iter);
	  --iter;
	}
      }
    }
  }
  // helices may not overlap - priority rules:
  // 5-helix > 4-helix > 3-helix
  for (int i = 0; i < (int) helix5.size(); ++i) {
    for (iter=helix4.begin(); iter!=helix4.end(); ++iter){
      if (*iter == helix5[i]) {
	helix4.erase(iter);
	--iter;
      }
    }
    for (iter=helix3.begin(); iter!=helix3.end(); ++iter){
      if (*iter == helix5[i]) {
	helix3.erase(iter);
	--iter;
      }
    }
  }
  for (int i = 0; i < (int) helix4.size(); ++i) {
    for (iter=helix3.begin(); iter!=helix3.end(); ++iter){
      if (*iter == helix4[i]) {
	helix3.erase(iter);
	--iter;
      }
    }
  }
} //end Dssp::filter_SecStruct()

void Dssp::writeToFiles(int n)
{
  for (int i=0; i < (int) helix3.size(); ++i) {
    timeseries3Helix << setw(10) << d_time << setw(10) << helix3[i]+1<< endl;
  }
  for (int i=0; i < (int) helix4.size(); ++i) {
    timeseries4Helix << setw(10) << d_time << setw(10) << helix4[i]+1<< endl;
  }
  for (int i=0; i < (int) helix5.size(); ++i) {
    timeseries5Helix << setw(10) << d_time << setw(10) << helix5[i]+1<< endl;
  }
  for (int i=0; i < (int) turn.size(); ++i) {
    timeseriesTurn << setw(10) << d_time << setw(10) << turn[i]+1<< endl;
  }
  for (int i=0; i < (int) extended.size(); ++i) {
    timeseriesBStrand << setw(10) << d_time << setw(10) << extended[i]+1 << endl;
  }
  for (int i=0; i < (int) bridge.size(); ++i) {
    timeseriesBBridge << setw(10) << d_time << setw(10) << bridge[i]+1 << endl;
  }
  for (int i=0; i < (int) Bend.size(); ++i) {
    timeseriesBend << setw(10) << d_time << setw(10) << Bend[i]+1 << endl;
  }

  d_time += d_dt*n;        

} //end Dssp::writeToFiles()


 void Dssp::settime(double i, double j)
 {
  d_time = i;
  d_dt = j;
 }

 void Dssp::opents(string fi1, string fi2, string fi3, string fi4, string fi5, string fi6, string fi7)
{
  timeseriesTurn.open(fi1.c_str());
  timeseries3Helix.open(fi2.c_str());
  timeseries4Helix.open(fi3.c_str());
  timeseries5Helix.open(fi4.c_str());
  timeseriesBBridge.open(fi5.c_str());
  timeseriesBStrand.open(fi6.c_str());
  timeseriesBend.open(fi7.c_str());
}

 void Dssp::readframe()
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

Dssp::Dssp(gcore::System &sys, args::Arguments &args)
{
  d_sys=&sys;
  d_args=&args;
  d_O = AtomSpecifier(sys);
  d_H = AtomSpecifier(sys);
  d_C = AtomSpecifier(sys);
  d_N = AtomSpecifier(sys);  
  d_CA = AtomSpecifier(sys); 
  opents("Turn.out", "3-Helix.out", "4-Helix.out", 
	 "5-Helix.out", "Beta-Bridge.out", "Beta-Strand.out", 
	 "Bend.out");
}