// args_ReferenceParser.cc

#include <cassert>
#include "ReferenceParser.h"
#include "Arguments.h"
#include "../fit/Reference.h"
#include "../gcore/System.h"
#include "../gcore/Molecule.h"

using gcore::System;
using fit::Reference;
using args::Arguments;

namespace args{

  // constructor
  ReferenceParser::ReferenceParser( 
				   gcore::System &sys,
				   const args::Arguments &args,
				   fit::Reference &ref)
    : 
    mySys(sys),
    myArgs(args),
    myRef(ref),
    myAdded(false)
  {
    std::vector<int> myMols;
  }

  void ReferenceParser::add_ref(){
  
    getMolecules();
    classes();
    atoms();

    // did we add anything at all?
    if(!myAdded){
      std::string errmsg = "Either \"class\" or \"atom\" "; 
      errmsg += "must be non-empty.\n";
      throw Arguments::Exception(errmsg);
    }
  }
  
  void ReferenceParser::getMolecules(){
    
    Arguments::const_iterator iter;
 
    // reset/initialize
    myMols.clear();

    // if no molecules supplied, add all
    if(myArgs.lower_bound("mol") == myArgs.upper_bound("mol"))
      for(int i = 0; i < mySys.numMolecules(); ++i)
        myMols.push_back(i);

    else{
      int thisMolNum = 0;

      // go through the args
      for(iter = myArgs.lower_bound("mol"); 
        iter != myArgs.upper_bound("mol");
        ++iter){

        thisMolNum = atoi(iter->second.c_str());

        // check
        if(thisMolNum > mySys.numMolecules()){
          std::string errmsg = "Supplied molecule index "; 
          errmsg += "is larger than the ";
          errmsg += "number of molecules in the system.";
          throw Arguments::Exception(errmsg);
        }

        // good
        myMols.push_back(thisMolNum - 1);
      }
    }
  }

  void ReferenceParser::classes(){

    Arguments::const_iterator iter;
    std::vector<int>::const_iterator mol;

    for(iter = myArgs.lower_bound("class"); 
      iter != myArgs.upper_bound("class");
      ++iter){
        for(mol = myMols.begin(); mol != myMols.end(); ++mol) 
          myRef.addClass(*mol,iter->second);
        myAdded = true;
    }
  }

  void ReferenceParser::atoms(){

    Arguments::const_iterator iter;
    int atomNum = 0;
    int molNum = 0;

    for(iter = myArgs.lower_bound("atoms"); 
      iter != myArgs.upper_bound("atoms");
      ++iter){

      atomNum = atoi(iter->second.c_str())-1;
      molNum = 0;
      while(atomNum >= mySys.mol(molNum).numAtoms()){
        atomNum -= mySys.mol(molNum).numAtoms();
        ++molNum;

        // check
        if(molNum == mySys.numMolecules()){
          std::string errmsg = "Supplied atom index is larger than the ";
          errmsg += "number of atoms in the system.";
          throw Arguments::Exception(errmsg);
        }
      }
      // good
      myRef.addAtom(molNum, atomNum);
      myAdded = true;
    }
  }
}
