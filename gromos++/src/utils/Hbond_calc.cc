#include <cassert>
#include <algorithm>
#include <vector>
#include <iterator>
#include <string>

#include "../args/Arguments.h"
#include "../bound/Boundary.h"
#include "../args/BoundaryParser.h"
#include "../gio/InG96.h"
#include "../gio/Ginstream.h"
#include "Neighbours.h"
#include "Hbond_calc.h"
#include "Hbond.h"
#include "CubeSystem.hcc"

using namespace bound;
using namespace gio;
using namespace args;

using utils::HB_calc;
using args::Arguments;

void HB_calc::setval(gcore::System& _sys, args::Arguments& _args) {
  sys = &_sys;
  args = &_args;

  pbc = BoundaryParser::boundary(_sys, _args);

  donors = AtomSpecifier(_sys);
  bound = AtomSpecifier(_sys);
  acceptors = AtomSpecifier(_sys);

  determineAtoms();
  //check for massfile
  string mfile;
  if (_args.count("massfile") > 0) {
    Arguments::const_iterator iter = _args.lower_bound("massfile");
    if (iter != _args.upper_bound("massfile")) {
      mfile = (iter->second.c_str());
    }
    readinmasses(mfile);
    determineAtomsbymass();
  }
  //all atoms that are not donors or acceptors have been removed, now populate donAB,...accB vectors with atomspecifier numbers:
  set_subatomspecs();

/*  cout << "DONORS " << donors.size() << ", " << num_A_donors << endl;
  for(int i=0; i< donors.size();++i)
    cout << i <<  " " << donors.mol(i) << ":" << donors.resnum(i) << " " << donors.resname(i) << ", " << donors.name(i)<< " " << donors.atom(i) << endl;
cout << "ACCEPTORS " << acceptors.size() << ", " << num_A_acceptors<< endl;
  for(int i=0; i< acceptors.size();++i)
    cout << i <<  " " << acceptors.mol(i) << ":" << acceptors.resnum(i) << " " << acceptors.resname(i) << ", " << acceptors.name(i) <<" " << acceptors.atom(i) << endl;
    */
}//end HB_calc::setval()

void HB_calc::readinmasses(std::string filename) {

  //new read in stuff, using the Ginstream...
  Ginstream nf(filename);
  vector<string> buffer;

  istringstream is;
  double mass = 0;

  while(!nf.getblock(buffer).eof()){
          //read in the hydrogen masses...
        if(buffer[0] == "HYDROGENMASS"){
            for (unsigned int j = 1; j < buffer.size() - 1; j++) {
                is.clear();
                is.str(buffer[j]);
                is >> mass;
                mass_hydrogens.push_back(mass);
            }

        }
        else if(buffer[0] == "ACCEPTORMASS"){
            for (unsigned int j = 1; j < buffer.size() - 1; j++) {
                is.clear();
                is.str(buffer[j]);
                is >> mass;
                mass_acceptors.push_back(mass);
            }

        }
        else if(buffer[0] == "DONORMASS"){ //optional. otherwise all atoms attached to an H will be used
            for (unsigned int j = 1; j < buffer.size() - 1; j++) {
                is.clear();
                is.str(buffer[j]);
                is >> mass;
                mass_donors.push_back(mass);
            }
        }
        else
            throw gromos::Exception("hbond", "Block " + buffer[0] + " not known!");
  }

  if (mass_hydrogens.empty())
    throw gromos::Exception("hbond", "HYDROGENMASS block in mass file is empty!");

  if (mass_acceptors.empty())
    throw gromos::Exception("hbond", "ACCEPTORMASS block in mass file is empty!");

}//end HB_calc::readinmasses()

void HB_calc::determineAtoms() {
  // we allow for the definition of four groups of atoms
  // 1. Donor atoms of group A
  // 2. Acceptor atoms of group A
  // 3. Donor atoms of group B
  // 4. Acceptor atoms of group B

  // everything can contain solvent. So to find out how many of those we have,
  // read a frame
  readframe();

  Arguments::const_iterator iter = args -> lower_bound("DonorAtomsA");
  Arguments::const_iterator to = args -> upper_bound("DonorAtomsA");

  for (; iter != to; iter++) {
    string spec = iter->second.c_str();
    donors.addSpecifier(spec);
  }

  iter = args -> lower_bound("AcceptorAtomsA");
  to = args -> upper_bound("AcceptorAtomsA");

  for (; iter != to; iter++) {
    string spec = iter->second.c_str();
    acceptors.addSpecifier(spec);
  }

  //sort them and find the atoms bound to the donor
  donors.sort();
  acceptors.sort();

  int m, a;
  for (int i = 0; i < donors.size(); ++i) {
    m = donors.mol(i);
    a = donors.atom(i);

    if (m < 0) {
      int j = a % sys->sol(0).topology().numAtoms();
      Neighbours neigh(*sys, 0, j, 0);
      //if we cannot find a neighbour: quit.
      if(neigh.empty())
            throw gromos::Exception("hbond","DonorAtomsA: " + donors.resname(i) + " has no bonding partner!");
      bound.addAtomStrict(-1, a - j + neigh[0]);
    } else {
      Neighbours neigh(*sys, m, a);
      if(neigh.empty()) //if we dont have a neighbour: quit
            throw gromos::Exception("hbond","DonorAtomsA: " + donors.resname(i) + " has no bonding partner!");
      bound.addAtomStrict(m, neigh[0]);
      //bound.addSolventType();
    }
  }

  // store how many acceptors and donor we have in A
  num_A_donors = donors.size();
  num_A_acceptors = acceptors.size();

  // if there is no B specified, we take A
  if (args->count("DonorAtomsB") <= 0 && args->count("AcceptorAtomsB") <= 0) {
    for (int i = 0; i < num_A_donors; i++) {
      donors.addAtomStrict(donors.mol(i), donors.atom(i));
      bound.addAtomStrict(bound.mol(i), bound.atom(i));
      if(donors.mol(i) < 0){
        donors.addSolventType();
      }
    }
    for (int i = 0; i < num_A_acceptors; i++) {
      acceptors.addAtomStrict(acceptors.mol(i), acceptors.atom(i));
      if(acceptors.mol(i) < 0)
        acceptors.addSolventType();
    }

  } else {
    AtomSpecifier donor_B(*sys);
    AtomSpecifier bound_B(*sys);
    AtomSpecifier acceptor_B(*sys);

    iter = args -> lower_bound("DonorAtomsB");
    to = args -> upper_bound("DonorAtomsB");

    for (; iter != to; iter++) {
      string spec = iter->second.c_str();
      donor_B.addSpecifier(spec);
    }
    iter = args -> lower_bound("AcceptorAtomsB");
    to = args -> upper_bound("AcceptorAtomsB");

    for (; iter != to; iter++) {
      string spec = iter->second.c_str();
      acceptor_B.addSpecifier(spec);
    }
    // and sort these as well and find the hydrogens bound to donor_B
    donor_B.sort();
    acceptor_B.sort();

    for (int i = 0; i < donor_B.size(); ++i) {
      m = donor_B.mol(i);
      a = donor_B.atom(i);
      if (m < 0) {
        int j = a % sys->sol(0).topology().numAtoms();
        Neighbours neigh(*sys, 0, j, 0);
        if(neigh.empty()) //if we dont have a neighbour
            throw gromos::Exception("hbond","DonorAtomsB: " + donor_B.resname(i) + " has no bonding partner!");
        bound_B.addAtomStrict(-1, a - j + neigh[0]);
      } else {
        Neighbours neigh(*sys, m, a);
        if(neigh.empty()) //if we dont have a neighbour
            throw gromos::Exception("hbond","DonorAtomsB: " + donor_B.resname(i) + " has no bonding partner!");
        bound_B.addAtomStrict(m, neigh[0]);
      }
    }

    // copy them into the d_donors, d_bound, d_acceptors
    for (int i = 0; i < donor_B.size(); ++i) {
      donors.addAtomStrict(donor_B.mol(i), donor_B.atom(i));
      if(donors.mol(i) < 0)
        donors.addSolventType();
    }
    for (int i =0 ;i < bound_B.size() ;++i ){
      bound.addAtomStrict(bound_B.mol(i), bound_B.atom(i));
      //if(bound.mol(i) < 0)
        //bound.addSolventType();
    }

    for (int i = 0; i < acceptor_B.size(); ++i) {
      acceptors.addAtomStrict(acceptor_B.mol(i), acceptor_B.atom(i));
      if(acceptors.mol(i) < 0)
        acceptors.addSolventType();
    }
  }
} //end HB_calc::determineAtoms()

void HB_calc::set_subatomspecs(){
    accAB.clear(); //holds all atoms that occur in AcceptorsA and B
    accA.clear(); //holds atoms that only occur in AcceptorsA
    accB.clear(); //holds atoms that only occur in AcceptorsB
    donAB.clear(); //holds all atoms that occur in DonorsA and B
    donA.clear(); //holds atoms that only occur in DonorsA
    donB.clear(); //holds atoms that only occur in DonorsB

    /*the following donor-acceptor combinations MUST BE USED in hbond calculation to ensure that all donor-acceptor pairs are included:
    (1) donAB & (accA + accB + accAB)
    (2) donA & (accB + accAB)
    (3) donB & (accA + accAB)
    */
    int i=0;
    for(int j=num_A_donors; j<donors.size(); ++j){ //go through all atoms in donorsB. search for donAB or donB atoms
        for(; i<num_A_donors; ++i){ //and search all atoms in donA
            if(donors.atom(i) == donors.atom(j) && donors.mol(i) == donors.mol(j))
                break;
        }
        if(i == num_A_donors){ //donB was not found in donA
            donB.push_back(j); //keep donorB atom and store in donB
            i=0; //reset donorA atom to 0
        }
        else{
            donAB.push_back(i); //donorB was found in donorA: use donorA atom and store in donAB, which holds atoms that are present in donorA and donorB
            ++i; //increase donorA atom by one, since the atomspec is ordererd, the next atom in donorB could be this atom in donorA
        }
        //increase j by 1 and start again
    }
    for(int j=0; j<num_A_donors; ++j){ //search donA for donA atoms
    	for(i=0; i<donAB.size(); ++i) //all atoms that are present in donorsA and donorsB are already in donAB, so every donorA atom that is not in donAB is stored on donA
    		if(donors.atom(donAB[i]) == donors.atom(j) && donors.mol(donAB[i]) == donors.mol(j))
    		     break;
    	if(i == donAB.size()){//donA atom was not found in donAB
    		donA.push_back(j);
  		}
    }
    i=0;
    //do the same for acceptors
    for(int j=num_A_acceptors; j<acceptors.size(); ++j){
        for(; i<num_A_acceptors; ++i){ //an search all atoms in accA
            if(acceptors.atom(i) == acceptors.atom(j) && acceptors.mol(i) == acceptors.mol(j))
                break;
        }
        if(i == num_A_acceptors){ //accB was not found in accA
            accB.push_back(j); //keep j
            i=0; //reset to 0
        }
        else{
            accAB.push_back(i); //otherwise use i
            ++i; //increase by one, since the atomspec is ordererd, the next atom in accB could be this in accA
        }
    }
    for(int j=0; j<num_A_acceptors; ++j){ //search accA for accA atoms
        for(i=0; i<accAB.size(); ++i)
            if(acceptors.atom(accAB[i]) == acceptors.atom(j) && acceptors.mol(accAB[i]) == acceptors.mol(j))
                 break;
        if(i == accAB.size())//accA atom was not found in accAB
            accA.push_back(j);
   }
/*
    cout << "=====\ndonAB:" << endl;
    for(int k=0; k<donAB.size();++k)
        cout << donors.resnum(donAB[k]) << endl;
	cout << "=====\ndonA:" << endl;
    for(int k=0; k<donA.size();++k)
        cout << donors.resnum(donA[k]) << endl;
	cout << "=====\ndonB:" << endl;
    for(int k=0; k<donB.size();++k)
        cout << donors.resnum(donB[k]) << endl;
    cout << "=====\naccAB:" << endl;
    for(int k=0; k<accAB.size();++k)
        cout << acceptors.resnum(accAB[k]) << endl;
	cout << "=====\naccA:" << endl;
    for(int k=0; k<accA.size();++k)
            cout << acceptors.resnum(accA[k]) << endl;
    cout << "=====\naccB:" << endl;
    for(int k=0; k<accB.size();++k)
            cout << acceptors.resnum(accB[k]) << endl;

*/
    if(donors.size() != 2*donAB.size() + donA.size() + donB.size())
        throw gromos::Exception("hbond","a problem with donors occured");
    if(acceptors.size() != 2*accAB.size() + accA.size() + accB.size())
        throw gromos::Exception("hbond","a problem with acceptors occured");
}

void HB_calc::determineAtomsbymass() {
  bool keep = false;
  //donors
  for (int i = 0; i < donors.size(); ++i) {
    keep = false;
    for (unsigned int j = 0; j < mass_hydrogens.size(); ++j) {
      if (donors.mass(i) == mass_hydrogens[j]) {
        keep = true;
        break;
      }
    }
    if(keep && !mass_donors.empty()){ //if we keep the H, we need to check if the donor was in the DONORMASS block (if something was specified in this block, otherwise keep by default)
        keep = false;
        for (unsigned int j = 0; j < mass_donors.size(); ++j) {
            if (bound.mass(i) == mass_donors[j]) {
                keep = true;
                break;
            }
        }
    }
    if (!keep) {
      donors.removeAtom(i);
      bound.removeAtom(i);
      if (i < num_A_donors) --num_A_donors;
      --i;
    }
  }
  //check if donorsA and B is empty:
  if(donors.empty())
    throw gromos::Exception("hbond","DonorAtomsA and DonorAtomsB do not contain any atoms. Please check your mass file.");

  // acceptors
  for (int i = 0; i < acceptors.size(); ++i) {
    keep = false;
    for (unsigned int j = 0; j < mass_acceptors.size(); ++j) {
      if (acceptors.mass(i) == mass_acceptors[j]) {
        keep = true;
        break;
      }
    }
    if (!keep) {
      acceptors.removeAtom(i);
      if (i < num_A_acceptors) --num_A_acceptors;
      --i;
    }
  }
  //check if acceptorA and B is empty:
  if(acceptors.empty())
    throw gromos::Exception("hbond","AcceptorAtomsA and AcceptorAtomsB do not contain any atoms. Please check your mass file.");
} //end HB_calc::determineAtomsbymass()

void HB_calc::readframe() {

  InG96 icc;

  try {
    args -> check("ref", 1);
    Arguments::const_iterator iterr = args -> lower_bound("ref");
    icc.open((iterr->second).c_str());
  } catch (const Arguments::Exception &) {
    args -> check("traj", 1);
    Arguments::const_iterator iterr = args -> lower_bound("traj");
    icc.open((iterr->second).c_str());
  }
  icc.select("ALL");
  icc >> *sys;
  icc.close();
}//end HB_calc::readframe()

void HB_calc::set_reduce(){
  if(reduce){
      int i, size, start;
      if(donors.numSolventAtoms()){ //if there was a solvent specified
          for(i=0; i < donors.size(); ++i) //find the first atom in the atomspecifier that denotes a solvent
            if(donors.mol(i) < 0)
                break;

          size = i; //store the first solvent atom of the atomspec
          start = donors.atom(i) % donors.numSolventAtoms(); //save the first atom number

          do{
              solv_donor.push_back(donors.atom(i) % donors.numSolventAtoms());  //store the solvent atom number (0,1,2,3..)
          }while( ++i < donors.size() && start != donors.atom(i) % donors.numSolventAtoms() ); //do this until the next solvent molecule starts

          solv_donor.push_back(size); //as last entry also store where the first solvent atom has been found.
      }

      if(acceptors.numSolventAtoms()){
          for(i=0; i < acceptors.size(); ++i)
            if(acceptors.mol(i) < 0)
                break;

          size = i; //store the first solvent atom of the atomspec
          start = acceptors.atom(i) % acceptors.numSolventAtoms();

          do{
              solv_acc.push_back(acceptors.atom(i) % acceptors.numSolventAtoms());  //store the solvent atom number (0,1,2,3..)
          }while( ++i < acceptors.size() && start != acceptors.atom(i) % acceptors.numSolventAtoms() );

          solv_acc.push_back(size);
      }
  }
}
