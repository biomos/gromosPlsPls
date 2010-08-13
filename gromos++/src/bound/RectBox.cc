// bound_RectBox.cc

#include <cmath>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <set>
#include "../gmath/Vec.h"
#include "../gcore/System.h"
#include "../gcore/Molecule.h"
#include "../gcore/LJException.h"
#include "../gcore/MoleculeTopology.h"
#include "../gcore/Bond.h"
#include "../gcore/Solvent.h"
#include "../gcore/SolventTopology.h"
#include "../gcore/Box.h"
#include "../fit/PositionUtils.h"
#include "RectBox.h"

using namespace std;
using bound::RectBox;
using gmath::Vec;
using namespace gcore;
using namespace utils;

static Vec nim(const Vec &r1,const  Vec &r2, const Box &box){

  Vec diff=r2-r1;
  Vec a;
  const double kabs = box.K().abs();
  a[0] = diff[0] - kabs * rint(diff[0]/kabs);
  const double labs = box.L().abs();
  a[1] = diff[1] - labs * rint(diff[1]/labs);
  const double mabs = box.M().abs();
  a[2] = diff[2] - mabs * rint(diff[2]/mabs);

  Vec rec = r1 + a;
  
  return rec;
  
}


Vec RectBox::nearestImage(const Vec &v1, const Vec &v2, const Box &box)const{
  return nim(v1, v2, box);
}

void RectBox::nogather(){

}

// gather based on a general list
void RectBox::gatherlist(){
    //std::cout << "# gather with an atom list " << std::endl;
    if (!sys().hasBox) throw gromos::Exception("Gather problem",
                              "System does not contain Box block! Abort!");

  if (sys().box().K().abs() == 0 || sys().box().L().abs() == 0 || sys().box().M().abs() == 0)
    throw gromos::Exception("Gather problem",
			    "Box block contains element(s) of value 0.0! Abort!");

  // gather the first molecule
    Molecule &mol=sys().mol(0);
    mol.pos(0)=nim(reference(0),mol.pos(0),sys().box());
    for(int j=1;j<mol.numPos();++j)
      mol.pos(j)=nim(mol.pos(j-1),mol.pos(j),sys().box());
    
    // gather the rest molecules
    // check whether the molecule should be gathered according to an atom list
    
    for(int i=1; i<sys().numMolecules();++i){

        Molecule &mol=sys().mol(i);
        int m=sys().primlist[i][0];
        int n=sys().primlist[i][1];
        int o=sys().primlist[i][2];

        mol.pos(m)=nim(sys().mol(n).pos(o),mol.pos(m),sys().box());

        if(m==0){
            for(int j=1;j<mol.numPos();++j)
                mol.pos(j)=nim(mol.pos(j-1),mol.pos(j),sys().box());
        }else{
            for(int j=m-1;j>=0;--j){
                mol.pos(j)=nim(mol.pos(j+1),mol.pos(j),sys().box());}
            for(int j=m+1;j<mol.numPos();++j){
            mol.pos(j)=nim(mol.pos(j-1),mol.pos(j),sys().box());}
        }
    }

    // now calculate cog
    Vec cog(0.0, 0.0, 0.0);
    int natom=0;
    for(int i=1; i<sys().numMolecules();++i){
        Molecule &mol=sys().mol(i);
        if(mol.numAtoms()>8)
            for(int a=0; a<mol.numAtoms(); ++a){
                cog += mol.pos(a);
                natom+=1;
            }
    }
    cog /= double(natom);


    for(int i=1; i<sys().numMolecules();++i){
        Molecule &mol=sys().mol(i);
        if(mol.numAtoms()<=8){
            mol.pos(0)=nim(cog,mol.pos(0),sys().box());
            for(int j=1;j<mol.numPos();++j)
                mol.pos(j)=nim(mol.pos(j-1),mol.pos(j),sys().box());
        }
    }

  // do the solvent
  Solvent &sol=sys().sol(0);
  for(int i=0;i<sol.numPos();i+= sol.topology().numAtoms()){
    //sol.pos(i)=nim(reference(0),sol.pos(i),sys().box());
    sol.pos(i)=nim(cog,sol.pos(i),sys().box());
    for (int j=i+1;j < (i+sol.topology().numAtoms());++j){
      sol.pos(j)=nim(sol.pos(j-1),sol.pos(j),sys().box());
    }
  }
};

// gather in term of time
void RectBox::gathertime(){
    //std::cout << "# gather with an atom list " << std::endl;
    if (!sys().hasBox) throw gromos::Exception("Gather problem",
                              "System does not contain Box block! Abort!");

    if (sys().box().K().abs() == 0 || sys().box().L().abs() == 0 || sys().box().M().abs() == 0)
        throw gromos::Exception("Gather problem",
			    "Box block contains element(s) of value 0.0! Abort!");

    if(refSys().sol(0).numPos()!=sys().sol(0).numPos())
        throw gromos::Exception("Gather problem",
			    "Number of solvent atoms are not equal in reference and the current system! Abort!");

    //std::cout << "# now gather the system with respect to the previous frame " << endl;
    //std::cout << "# number of mol sys       " << sys().numMolecules() << endl;
    //std::cout << "# number of mol newrefSys " << refSys().numMolecules() << endl;
    for(int i=0;i<sys().numMolecules();++i){
        Molecule &mol=sys().mol(i);
        Molecule &refmol=refSys().mol(i);
        mol.pos(0)=nim(refmol.pos(0),mol.pos(0),sys().box());
        refmol.pos(0)=mol.pos(0);
        for(int j=1;j<mol.numPos();++j){
            mol.pos(j)=nim(mol.pos(j-1),mol.pos(j),sys().box());
            refmol.pos(j)=mol.pos(j);
        }
    }

    // correct for ions
    Vec cog(0.,0.,0.);
    int count=0;
    for(int i=0;i<sys().numMolecules();++i){
        Molecule &mol=sys().mol(i);
        if(mol.numPos()>8)
            for(int j=0;j<mol.numPos();++j){
                cog+=mol.pos(j);
                count+=1;
            }
    }
    cog/=double(count);

    for(int i=0;i<sys().numMolecules();++i){
        Molecule &mol=sys().mol(i);
        if(mol.numPos()<=8){
            Molecule &refmol=refSys().mol(i);
            mol.pos(0)=nim(cog,mol.pos(0),sys().box());
            refmol.pos(0)=mol.pos(0);
            for(int j=1;j<mol.numPos();++j){
                mol.pos(j)=nim(mol.pos(j-1),mol.pos(j),sys().box());
                refmol.pos(j)=mol.pos(j);
            }
        }
    }

    // do the solvent
    Solvent &sol=sys().sol(0);
    Solvent &refsol=refSys().sol(0);

    //std::cout << "# working on solv molecule " << sol.numPos() << endl;
    //std::cout << "# working on refsolv molecule " << refsol.numPos() << endl;
    if(refsol.numPos()==sol.numPos())
        for(int i=0;i<sol.numPos();i+= sol.topology().numAtoms()){
            //sol.pos(i)=nim(reference(0),sol.pos(i),sys().box());
            //std::cout << "# working on solv molecule " << i << endl;
            sol.pos(i)=nim(refsol.pos(i),sol.pos(i),sys().box());
            refsol.pos(i)=sol.pos(i);
            for (int j=i+1;j < (i+sol.topology().numAtoms());++j){
                //sol.pos(j)=nim(sol.pos(j-1),sol.pos(j),sys().box());
                sol.pos(j)=nim(sol.pos(j-1),sol.pos(j),sys().box());
                //sol.pos(j)=nim(refsol.pos(j),sol.pos(j),sys().box());
                //sol.pos(j)=nim(newreference(ref),sol.pos(j),sys().box());
                refsol.pos(j)=sol.pos(j);
            }
        }
        else{
            std::cout << "# solv num " << sol.numPos()
                    << " and refsolv num " << refsol.numPos()
                    << " are not equal. solv gathering based on cog : " << std::endl;

            for(int i=0;i<sol.numPos();i+= sol.topology().numAtoms()){
                sol.pos(i)=nim(cog,sol.pos(i),sys().box());
                for (int j=i+1;j < (i+sol.topology().numAtoms());++j){
                    sol.pos(j)=nim(sol.pos(j-1),sol.pos(j),sys().box());
//                    refsol.pos(j)=sol.pos(j);
                }
            }
        }

};

// gather the 1st frame based on an atom list, then the rest in term of time
// everytime we update the reference system, thus to avoid changing the code for an old state
void RectBox::gatherltime(){
    //std::cout << "# gather with an atom list " << std::endl;
    if (!sys().hasBox) throw gromos::Exception("Gather problem",
                              "System does not contain Box block! Abort!");

    if (sys().box().K().abs() == 0 || sys().box().L().abs() == 0 || sys().box().M().abs() == 0)
        throw gromos::Exception("Gather problem",
			    "Box block contains element(s) of value 0.0! Abort!");

    if(sys().primlist[0][0]==31415926){
        //std::cout << "# now gather the system with respect to the previous frame " << endl;
        //std::cout << "# number of mol sys       " << sys().numMolecules() << endl;
        //std::cout << "# number of mol newrefSys " << refSys().numMolecules() << endl;
        for(int i=0;i<sys().numMolecules();++i){
            Molecule &mol=sys().mol(i);
            Molecule &refmol=refSys().mol(i);
            mol.pos(0)=nim(refmol.pos(0),mol.pos(0),sys().box());
            refmol.pos(0)=mol.pos(0);
            for(int j=1;j<mol.numPos();++j){
                mol.pos(j)=nim(mol.pos(j-1),mol.pos(j),sys().box());
                refmol.pos(j)=mol.pos(j);
            }
        }

        // correct for ions
        Vec cog(0.,0.,0.);
        int count=0;
        for(int i=0;i<sys().numMolecules();++i){
            Molecule &mol=sys().mol(i);
            if(mol.numPos()>8)
                for(int j=0;j<mol.numPos();++j){
                    cog+=mol.pos(j);
                    count+=1;
                }
        }
        cog/=double(count);

        for(int i=0;i<sys().numMolecules();++i){
            Molecule &mol=sys().mol(i);
            if(mol.numPos()<=8){
                Molecule &refmol=refSys().mol(i);
                mol.pos(0)=nim(cog,mol.pos(0),sys().box());
                refmol.pos(0)=mol.pos(0);
                for(int j=1;j<mol.numPos();++j){
                    mol.pos(j)=nim(mol.pos(j-1),mol.pos(j),sys().box());
                    refmol.pos(j)=mol.pos(j);
                }
            }
        }

        // do the solvent
        Solvent &sol=sys().sol(0);
        Solvent &refsol=refSys().sol(0);

        if(refsol.numPos()==sol.numPos())
         for(int i=0;i<sol.numPos();i+= sol.topology().numAtoms()){
            //sol.pos(i)=nim(reference(0),sol.pos(i),sys().box());
            //std::cout << "# working on solv molecule " << i << endl;
            sol.pos(i)=nim(refsol.pos(i),sol.pos(i),sys().box());
            refsol.pos(i)=sol.pos(i);
            for (int j=i+1;j < (i+sol.topology().numAtoms());++j){
                //sol.pos(j)=nim(sol.pos(j-1),sol.pos(j),sys().box());
                sol.pos(j)=nim(sol.pos(j-1),sol.pos(j),sys().box());
                //sol.pos(j)=nim(refsol.pos(j),sol.pos(j),sys().box());
                //sol.pos(j)=nim(newreference(ref),sol.pos(j),sys().box());
                refsol.pos(j)=sol.pos(j);
            }
        }
        else{
            std::cout << "# solv num " << sol.numPos()
                    << " and refsolv num " << refsol.numPos()
                    << " are not equal. solv gathering based on cog : " << endl;

            for(int i=0;i<sol.numPos();i+= sol.topology().numAtoms()){
                sol.pos(i)=nim(cog,sol.pos(i),sys().box());
                for (int j=i+1;j < (i+sol.topology().numAtoms());++j){
                    sol.pos(j)=nim(sol.pos(j-1),sol.pos(j),sys().box());
//                    refsol.pos(j)=sol.pos(j);
                }
            }
        }
    }
    else {
        //std::cout << "# now gather the system based on the atom list " << endl;
        //RectBox::gatherlist();

        // gather the first molecule
        Molecule &mol=sys().mol(0);
        Molecule &refmol=refSys().mol(0);
        mol.pos(0)=nim(reference(0),mol.pos(0),sys().box());
        refmol.pos(0)=mol.pos(0);
        for(int j=1;j<mol.numPos();++j){
            mol.pos(j)=nim(mol.pos(j-1),mol.pos(j),sys().box());
            refmol.pos(j)=mol.pos(j);
        }

        // gather the rest molecules
        // check whether the molecule should be gathered according to an atom list
        for(int i=1; i<sys().numMolecules();++i){
            Molecule &mol=sys().mol(i);
            Molecule &refmol=refSys().mol(i);
            int m=sys().primlist[i][0];
            int n=sys().primlist[i][1];
            int o=sys().primlist[i][2];
            //if(n!=i-1 || o!=0)

            //    cout << "# mol " << i << " atom " << m << " ref mol " << n << " atom " << o << endl;
            mol.pos(m)=nim(sys().mol(n).pos(o),mol.pos(m),sys().box());
            if(m==0){
                for(int j=1;j<mol.numPos();++j){
                    mol.pos(j)=nim(mol.pos(j-1),mol.pos(j),sys().box());
                    refmol.pos(j)=mol.pos(j);
                }
            }else{
                for(int j=m-1;j>=0;--j){
                    mol.pos(j)=nim(mol.pos(j+1),mol.pos(j),sys().box());
                    refmol.pos(j)=mol.pos(j);
                }
                for(int j=m+1;j<mol.numPos();++j){
                    mol.pos(j)=nim(mol.pos(j-1),mol.pos(j),sys().box());
                    refmol.pos(j)=mol.pos(j);
                }
            }
        }

        // now calculate cog
        Vec cog(0.0, 0.0, 0.0);
        int natom=0;
        for(int i=1; i<sys().numMolecules();++i){
            Molecule &mol=sys().mol(i);
            if(mol.numAtoms()>8)
                for(int a=0; a<mol.numAtoms(); ++a){
                    cog += mol.pos(a);
                    natom+=1;
                }
        }

        cog /= double(natom);
        for(int i=1; i<sys().numMolecules();++i){
            Molecule &mol=sys().mol(i);
            //Molecule &refmol=refSys().mol(i);
            if(mol.numAtoms()<=8){
                mol.pos(0)=nim(cog,mol.pos(0),sys().box());
                //refmol.pos(0)=mol.pos(0);
                for(int j=1;j<mol.numPos();++j){
                    mol.pos(j)=nim(mol.pos(j-1),mol.pos(j),sys().box());
                    //refmol.pos(j)=mol.pos(j);
                }
            }
        }

        // do the solvent
        Solvent &sol=sys().sol(0);
        Solvent &refsol=refSys().sol(0);
        for(int i=0;i<sol.numPos();i+= sol.topology().numAtoms()){
            sol.pos(i)=nim(cog,sol.pos(i),sys().box());
            for (int j=i+1;j < (i+sol.topology().numAtoms());++j){
                sol.pos(j)=nim(sol.pos(j-1),sol.pos(j),sys().box());
                refsol.pos(i)=sol.pos(i);
            }
        }
    }
};

// gather based on a reference structure
void RectBox::gatherref(){
    if (!sys().hasBox) throw gromos::Exception("Gather problem",
                              "System does not contain Box block! Abort!");

    if (sys().box().K().abs() == 0 || sys().box().L().abs() == 0 || sys().box().M().abs() == 0)
        throw gromos::Exception("Gather problem",
			    "Box block contains element(s) of value 0.0! Abort!");

    if (sys().numMolecules() != refSys().numMolecules())
        throw gromos::Exception("Gather problem",
            "Number of SOLUTE  molecules in reference and frame are not the same.");
    if (sys().sol(0).numPos() != refSys().sol(0).numPos())
        throw gromos::Exception("Gather problem",
            "Number of SOLVENT molecules in reference and frame are not the same.");


  for (int i = 0; i < sys().numMolecules(); ++i) {
    Molecule &mol = sys().mol(i);
    Molecule &refMol = refSys().mol(i);
    for (int j = 0; j < mol.numPos(); ++j) {
      // gather to the reference
      mol.pos(j) = nim(refMol.pos(j), mol.pos(j), sys().box());
    }
  }
  // do the solvent
  Solvent &sol = sys().sol(0);
  Solvent &refSol = refSys().sol(0);

  Vec solcog(0.,0.,0.);
  for (int i=0; i<refSol.numPos(); i+=sol.topology().numAtoms()) {
    //sol.pos(i) = nim(refSol.pos(i), sol.pos(i), sys().box());
    for(int j=i;j< (i + sol.topology().numAtoms());++j){
        solcog+=refSol.pos(j);
    }
  }
  solcog /= refSol.numPos();

  for (int i=0; i<sol.numPos(); i+=sol.topology().numAtoms()) {
    //sol.pos(i) = nim(refSol.pos(i), sol.pos(i), sys().box());
    sol.pos(i) = nim(solcog, sol.pos(i), sys().box());
    for(int j=i+1;j< (i + sol.topology().numAtoms());++j){
        sol.pos(j) = nim(sol.pos(j-1),sol.pos(j),sys().box());
    }
  }
};

void RectBox::gatherrtime(){
  if (!sys().hasBox) throw gromos::Exception("Gather problem",
          "System does not contain Box block! Abort!");

  const gcore::Box & box = sys().box();

  if (sys().box().K().abs() == 0 || sys().box().L().abs() == 0 || sys().box().M().abs() == 0)
    throw gromos::Exception("Gather problem",
          "Box block contains element(s) of value 0.0! Abort!");


  if (sys().numMolecules() != refSys().numMolecules())
    throw gromos::Exception("Gather problem",
            "Number of molecules in reference and frame are not the same.");
  for (int i = 0; i < sys().numMolecules(); ++i) {
    Molecule &mol = sys().mol(i);
    Molecule &refMol = refSys().mol(i);
    for (int j=0; j<mol.numPos(); ++j) {
      // gather to the reference
      mol.pos(j) = nim(refMol.pos(j), mol.pos(j), box);
      // set the current frame as the reference for the next frame
      refMol.pos(j) = mol.pos(j);
    }
  }
  // do the solvent
  Solvent &sol = sys().sol(0);
  Solvent &refSol = refSys().sol(0);
  //std::cout << "# working on solv molecule " << sol.numPos() << endl;
  //std::cout << "# working on refsolv molecule " << refSol.numPos() << endl;
  if (sol.numPos() != refSol.numPos()) {
    throw gromos::Exception("Gather problem", "Number of solvent positions in "
            "reference and frame are not the same.");
  }
  for (int i=0; i<sol.numPos(); i+=sol.topology().numAtoms()) {
    sol.pos(i) = nim(refSol.pos(i), sol.pos(i), sys().box());
    refSol.pos(i) = sol.pos(i);
    //sol.pos(i) = nim(solcog, sol.pos(i), sys().box());
    for(int j=i+1;j< (i + sol.topology().numAtoms());++j){
        sol.pos(j) = nim(sol.pos(j-1),sol.pos(j),sys().box());
        refSol.pos(j) = sol.pos(j);
    }
  }
}

// gather based on bond connection
void RectBox::gatherbond(){
  if (!sys().hasBox) throw gromos::Exception("Gather problem",
                              "System does not contain Box block! Abort!");

  if (sys().box().K().abs() == 0 || sys().box().L().abs() == 0 || sys().box().M().abs() == 0)
    throw gromos::Exception("Gather problem",
			    "Box block contains element(s) of value 0.0! Abort!");

  for(int i=0; i<sys().numMolecules();++i){
    Molecule &mol=sys().mol(i);
    mol.pos(0)=nim(reference(i),mol.pos(0),sys().box());
    for(int j=1;j<mol.numPos();++j){

      //find a previous atom to which we are bound
      BondIterator bi(mol.topology());
      int k=0;
      for(;bi;++bi)
	if(bi()[1]==j) { k = bi()[0]; break; }
      mol.pos(j)=nim(mol.pos(k),mol.pos(j),sys().box());
    }
  }

  // do the solvent
  Solvent &sol=sys().sol(0);
  for(int i=0;i<sol.numPos();i+= sol.topology().numAtoms()){
    sol.pos(i)=nim(reference(0),sol.pos(i),sys().box());
    for (int j=i+1;j < (i+sol.topology().numAtoms());++j){
      sol.pos(j)=nim(sol.pos(j-1),sol.pos(j),sys().box());
    }
  }
};


void RectBox::gather(){

  if (!sys().hasBox) throw gromos::Exception("Gather problem",  
                              "System does not contain Box block! Abort!");

  if (sys().box().K().abs() == 0 || sys().box().L().abs() == 0 || sys().box().M().abs() == 0)
    throw gromos::Exception("Gather problem",  
			    "Box block contains element(s) of value 0.0! Abort!");  

  for(int i=0; i<sys().numMolecules();++i){
    Molecule &mol=sys().mol(i);
    mol.pos(0)=nim(reference(0),mol.pos(0),sys().box());
    for(int j=1;j<mol.numPos();++j)
      mol.pos(j)=nim(mol.pos(j-1),mol.pos(j),sys().box());
  }
  // do the solvent 
  Solvent &sol=sys().sol(0);
  for(int i=0;i<sol.numPos();i+= sol.topology().numAtoms()){
    sol.pos(i)=nim(reference(0),sol.pos(i),sys().box());  
    for (int j=i+1;j < (i+sol.topology().numAtoms());++j){
      sol.pos(j)=nim(sol.pos(j-1),sol.pos(j),sys().box());
    }
  }
}

void RectBox::gathergr(){

  if (!sys().hasBox) 
    throw gromos::Exception("Gather problem",  
			    "System does not contain Box block! Abort!");

   if (sys().box().K().abs() == 0 || sys().box().L().abs() == 0 || sys().box().M().abs() == 0)
     throw gromos::Exception("Gather problem",  
			     "Box block contains element(s) of value 0.0! Abort!");
   
   for(int i=0; i<sys().numMolecules();++i){
     Molecule &mol=sys().mol(i);
     mol.pos(0)=nim(reference(i),mol.pos(0),sys().box());
     for(int j=1;j<mol.numAtoms();++j)
       mol.pos(j)=nim(mol.pos(j-1),mol.pos(j),sys().box());
   }
}

void RectBox::gathermgr(){

  if (!sys().hasBox) 
    throw gromos::Exception("Gather problem",  
			    "System does not contain Box block! Abort!");

   if (sys().box().K().abs() == 0 || sys().box().L().abs() == 0 || sys().box().M().abs() == 0)
     throw gromos::Exception("Gather problem",  
			     "Box block contains element(s) of value 0.0! Abort!");
   
   const Vec centre(0.5 * sys().box().K().abs(), 0.5 * sys().box().L().abs(), 0.5 * sys().box().M().abs());

   for(int i=0; i<sys().numMolecules();++i){
     Molecule &mol=sys().mol(i);
     mol.pos(0)=nim(reference(i),mol.pos(0),sys().box());
     for(int j=1;j<mol.numAtoms();++j)
       mol.pos(j)=nim(mol.pos(j-1),mol.pos(j),sys().box());

     // now calculate cog
     Vec cog(0.0, 0.0, 0.0);
     for(int a=0; a<mol.numAtoms(); ++a)
       cog += mol.pos(a);
     cog /= mol.numAtoms();
     Vec cog_box = nim(centre, cog, sys().box());
     Vec trans = cog_box - cog;
     fit::PositionUtils::translate(mol, trans);
   }
}


void RectBox::coggather(){

    //std::cout << "# gather with cog " << std::endl;
  if (!sys().hasBox)
    throw gromos::Exception("Gather problem",  
			    "System does not contain Box block! Abort!");
  
  if (sys().box().K().abs() == 0 || sys().box().L().abs() == 0 || sys().box().M().abs() == 0)
    throw gromos::Exception("Gather problem",  
			    "Box block contains element(s) of value 0.0! Abort!");
  
  Molecule &mol=sys().mol(0);
  Solvent &sol=sys().sol(0);
  
  Vec ref(0.0,0.0,0.0);
  Vec cog;
  int atoms=0;
  
  // do mol(0) with respect to ref (0,0,0)
  mol.pos(0)=nim(ref,mol.pos(0),sys().box());
  for(int j=1;j<mol.numAtoms();++j){
    mol.pos(j)=nim(mol.pos(j-1),mol.pos(j),sys().box());}
  
  // calculate COG of mol(0)
  for (int i=0;i < mol.numAtoms(); i++) {
    cog = cog + mol.pos(i);
    ++atoms;
  }
  cog = (1.0/double(atoms))*cog;
  
  // do the rest of the molecules
  for(int i=1;i<sys().numMolecules();++i){
    Molecule &mol=sys().mol(i);
    mol.pos(0)=nim(cog,mol.pos(0),sys().box());      
    for(int j=1;j<mol.numPos();++j){
      mol.pos(j)=nim(mol.pos(j-1),mol.pos(j),sys().box());
    }
  }
  
  // do the solvent 
  for(int i=0;i<sol.numPos();i+= sol.topology().numAtoms()){
    sol.pos(i)=nim(cog,sol.pos(i),sys().box());  
    for (int j=i+1;j < (i+sol.topology().numAtoms());++j){
      sol.pos(j)=nim(sol.pos(j-1),sol.pos(j),sys().box());
    }
  }
  
} 

void RectBox::crsgather(){

  if (!sys().hasBox)
    throw gromos::Exception("Gather problem",
               "System does not contain Box block! Abort!");

  if (sys().box().K().abs() == 0 || sys().box().L().abs() == 0 || sys().box().M().abs() == 0)
    throw gromos::Exception("Gather problem",
                "Box block contains element(s) of value 0.0! Abort!");

  // Reconstruct the connectivity of the submolecules
  for(int i=0; i<sys().numMolecules();++i){
    Molecule &mol=sys().mol(i);
    mol.pos(0)=nim(reference(i),mol.pos(0),sys().box());
    for(int j=1;j<mol.numAtoms();++j)
      mol.pos(j)=nim(mol.pos(j-1),mol.pos(j),sys().box());
  }

  // Determine the positions of the centres of geometry of the gathered molecules 
  // and store them in vcog
  std::vector<Vec> vcog;
  for(int i=0; i<sys().numMolecules();++i){
    Vec cog(0.0,0.0,0.0);
    int numat=0;
    for (int j=0; j<sys().mol(i).numAtoms(); ++j) {
      cog = cog + sys().mol(i).pos(j);
      ++numat;
    }
    cog = (1.0/double(numat))*cog;
    vcog.push_back(cog);
  }

  // Gather nearest image of cog of molecule 1 w.r.t. origin
  // vcog[0]=nim((0.0,0.0,0.0),vcog[0],sys().box());

  // Now gather cog's w.r.t. cog of previous molecule
  // ocog: buffer to store the overall cog of already gathered cog's
  Vec ocog=vcog[0];
  for(int i=0; i<sys().numMolecules()-1;++i){
    Vec nimj=nim(vcog[i],vcog[i+1],sys().box());
    vcog[i+1]=nimj;
    ocog+=vcog[i+1];
  }

  // Now gather the atoms of the solute molecules with
  // the newly determined cog's of the respective molecule
  for(int i=0;i<sys().numMolecules();++i){
    Molecule &mol=sys().mol(i);
    for(int j=0;j<mol.numAtoms();++j){
      mol.pos(j)=nim(vcog[i],mol.pos(j),sys().box());
    }
  }

  // Gather the solvent molecules with ocog as a reference
  ocog=ocog/double(sys().numMolecules());
  Solvent &sol=sys().sol(0);
  for(int i=0;i<sol.numPos();i+=sol.topology().numAtoms()){
    sol.pos(i)=nim(ocog,sol.pos(i),sys().box());
    for(int j=i+1;j<(i+sol.topology().numAtoms());++j){
      sol.pos(j)=nim(sol.pos(j-1),sol.pos(j),sys().box());
    }
  }
}

void RectBox::seqgather(){

  if (!sys().hasBox)
    throw gromos::Exception("Gather problem",
               "System does not contain Box block! Abort!");

  if (sys().box().K().abs() == 0 || sys().box().L().abs() == 0 || sys().box().M().abs() == 0)
    throw gromos::Exception("Gather problem",
                "Box block contains element(s) of value 0.0! Abort!");

  // Reconstruct the connectivity of the submolecules
  for(int i=0; i<sys().numMolecules();++i){
    Molecule &mol=sys().mol(i);
    mol.pos(0)=nim(reference(i),mol.pos(0),sys().box());
    for(int j=1;j<mol.numAtoms();++j)
      mol.pos(j)=nim(mol.pos(j-1),mol.pos(j),sys().box());
  }

  // Determine the positions of the centres of geometry of the gathered molecules 
  // and store them in vcog
  std::vector<Vec> vcog;
  for(int i=0; i<sys().numMolecules();++i){
    Vec cog(0.0,0.0,0.0);
    int numat=0;
    for (int j=0; j<sys().mol(i).numAtoms(); ++j) {
      cog = cog + sys().mol(i).pos(j);
      ++numat;
    }
    cog = (1.0/double(numat))*cog;
    vcog.push_back(cog);
  }

  // Gather nearest image of cog of molecule 1 w.r.t. origin
  // vcog[0]=nim((0.0,0.0,0.0),vcog[0],sys().box());

  // Now gather cog's w.r.t. cog of previous molecule
  // ocog: buffer to store the overall cog of already gathered cog's
  Vec ocog=vcog[0];
  for(int i=0; i<sys().numMolecules()-1;++i){
    // crs:
    // Vec nimj=nim(vcog[i],vcog[i+1],sys().box());
    // seq:
    Vec nimj=nim(ocog/double(i+1),vcog[i+1],sys().box());
    vcog[i+1]=nimj;
    ocog+=vcog[i+1];
  }

  // Now gather the atoms of the solute molecules with
  // the newly determined cog's of the respective molecule
  for(int i=0;i<sys().numMolecules();++i){
    Molecule &mol=sys().mol(i);
    for(int j=0;j<mol.numAtoms();++j){
      mol.pos(j)=nim(vcog[i],mol.pos(j),sys().box());
    }
  }

  // Gather the solvent molecules with ocog as a reference
  ocog=ocog/double(sys().numMolecules());
  Solvent &sol=sys().sol(0);
  for(int i=0;i<sol.numPos();i+=sol.topology().numAtoms()){
    sol.pos(i)=nim(ocog,sol.pos(i),sys().box());
    for(int j=i+1;j<(i+sol.topology().numAtoms());++j){
      sol.pos(j)=nim(sol.pos(j-1),sol.pos(j),sys().box());
    }
  }
}

void RectBox::gengather(){

  if (!sys().hasBox)
    throw gromos::Exception("Gather problem",
               "System does not contain Box block! Abort!");

  if (sys().box().K().abs() == 0 || sys().box().L().abs() == 0 || sys().box().M().abs() == 0)
    throw gromos::Exception("Gather problem",
                "Box block contains element(s) of value 0.0! Abort!");

  // Reconstruct the connectivity of the submolecules
  for(int i=0; i<sys().numMolecules();++i){
    Molecule &mol=sys().mol(i);
    mol.pos(0)=nim(reference(i),mol.pos(0),sys().box());
    for(int j=1;j<mol.numAtoms();++j)
      mol.pos(j)=nim(mol.pos(j-1),mol.pos(j),sys().box());
  }

  // Determine the positions of the centres of geometry of the gathered molecules 
  // and store them in vcog
  std::vector<Vec> vcog;
  for(int i=0; i<sys().numMolecules();++i){
    Vec cog(0.0,0.0,0.0);
    int numat=0;
    for (int j=0; j<sys().mol(i).numAtoms(); ++j) {
      cog = cog + sys().mol(i).pos(j);
      ++numat;
    }
    cog = (1.0/double(numat))*cog;
    vcog.push_back(cog);
  }

  // Gather nearest image of cog of molecule 1 w.r.t. origin
  // vcog[0]=nim((0.0,0.0,0.0),vcog[0],sys().box());

  // Use vcogi to make the graph connecting the closest cog's
  std::vector<int> vcogi;
  for(int i=0; i<sys().numMolecules();++i){
    vcogi.push_back(i);
  }

  // Now gather cog's w.r.t. each other
  // ocog: buffer to store the overall cog of already gathered cog's
  Vec ocog=vcog[0];
  for(int i=0; i<sys().numMolecules()-1;++i){
    // Determine closest nim to i among remaining molecules (using vcogi)
    int bufi=vcogi[i];
    int inimcogi=vcogi[i+1];
    Vec nimcogi=nim(vcog[bufi],vcog[inimcogi],sys().box());
    int jclose=i+1;
    for(int j=i+2; j<sys().numMolecules();++j){
      int bufj=vcogi[j];
      if( (nim(vcog[bufi],vcog[bufj],sys().box())-vcog[bufi]).abs()<(nimcogi-vcog[bufi]).abs()){
        nimcogi=nim(vcog[bufi],vcog[bufj],sys().box());
        inimcogi=bufj;
        jclose=j;
      }
    }
    // Now swap inimcogi with i+1 in vcogi
    int bufci=vcogi[i+1];
    vcogi[i+1]=inimcogi;
    vcogi[jclose]=bufci;
    
    // Set vcog[i+1] either to its nim to vcog[i], or to
    // nim to overall cog of molecules[1 ... i], depending
    // on what corresponds with the closest distance
    Vec nic1=nimcogi;
    Vec nic2=nim(ocog/double(i+1),nimcogi,sys().box());
    if((nic1-vcog[bufi]).abs()<(nic2-ocog/double(i+1)).abs()){
      vcog[inimcogi]=nic1;
    }
    else{
      vcog[inimcogi]=nic2;
    }
    ocog+=vcog[inimcogi];
  }

  // Now gather the atoms of the solute molecules with
  // the newly determined cog's of the respective molecule
  // as a reference
  for(int i=0;i<sys().numMolecules();++i){
    Molecule &mol=sys().mol(i);
    for(int j=0;j<mol.numAtoms();++j){
      mol.pos(j)=nim(vcog[i],mol.pos(j),sys().box());
    }
  }

  // Gather the solvent molecules with ocog as a reference
  ocog=ocog/double(sys().numMolecules());
  Solvent &sol=sys().sol(0);
  for(int i=0;i<sol.numPos();i+=sol.topology().numAtoms()){
    sol.pos(i)=nim(ocog,sol.pos(i),sys().box());
    for(int j=i+1;j<(i+sol.topology().numAtoms());++j){
      sol.pos(j)=nim(sol.pos(j-1),sol.pos(j),sys().box());
    }
  }
}

void RectBox::bondgather(){

  if (!sys().hasBox) throw gromos::Exception("Gather problem",  
                              "System does not contain Box block! Abort!");

  if (sys().box().K().abs() == 0 || sys().box().L().abs() == 0 || sys().box().M().abs() == 0)
    throw gromos::Exception("Gather problem",  
			    "Box block contains element(s) of value 0.0! Abort!");  

  for(int i=0; i<sys().numMolecules();++i){
    Molecule &mol=sys().mol(i);
    mol.pos(0)=nim(reference(0),mol.pos(0),sys().box());
    for(int j=1;j<mol.numPos();++j){
      
      //find a previous atom to which we are bound
      BondIterator bi(mol.topology());
      int k=0;
      for(;bi;++bi)
	if(bi()[1]==j) { k = bi()[0]; break; }
      mol.pos(j)=nim(mol.pos(k),mol.pos(j),sys().box());
    }
  }
  
  // do the solvent 
  Solvent &sol=sys().sol(0);
  for(int i=0;i<sol.numPos();i+= sol.topology().numAtoms()){
    sol.pos(i)=nim(reference(0),sol.pos(i),sys().box());  
    for (int j=i+1;j < (i+sol.topology().numAtoms());++j){
      sol.pos(j)=nim(sol.pos(j-1),sol.pos(j),sys().box());
    }
  }
}

void RectBox::refgather() {
  if (!sys().hasBox) throw gromos::Exception("Gather problem",
          "System does not contain Box block! Abort!");

  const gcore::Box & box = sys().box();

  if (sys().box().M().abs() == 0 || sys().box().L().abs() == 0 || sys().box().M().abs() == 0)
    throw gromos::Exception("Gather problem",
          "Box block contains element(s) of value 0.0! Abort!");


  if (sys().numMolecules() != refSys().numMolecules())
    throw gromos::Exception("Gather problem",
            "Number of molecules in reference and frame are not the same.");
  for (int m = 0; m < sys().numMolecules(); ++m) {
    Molecule &mol = sys().mol(m);
    Molecule &refMol = refSys().mol(m);
    if (mol.numPos() != refMol.numPos()) {
      std::ostringstream msg;
      msg << "Number of positions in molecule " << m + 1 << " in reference and frame"
              " are not the same.";
      throw gromos::Exception("Gather problem", msg.str());
    }
    for (int a = 0; a < mol.numPos(); ++a) {
      // gather to the reference
      mol.pos(a) = nim(refMol.pos(a), mol.pos(a), box);
      // set the current frame as the reference for the next frame
      refMol.pos(a) = mol.pos(a);
    }
  }
  // do the solvent
  Solvent &sol = sys().sol(0);
  Solvent &refSol = sys().sol(0);
  if (sol.numPos() != refSol.numPos()) {
    throw gromos::Exception("Gather problem", "Number of solvent positions in "
            "reference and frame are not the same.");
  }
  for (int a = 0; a < sol.numPos(); ++a) {
    // gather to the reference
    refSol.pos(a) = nim(refSol.pos(a), sol.pos(a), box);
    // set the current frame as the reference for the next frame;
    refSol.pos(a) = sol.pos(a);
  }
}
