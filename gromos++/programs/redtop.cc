#include <iostream>
#include <strstream>
#include "../src/args/Arguments.h"
#include "../src/gio/InTopology.h"
#include "../src/gio/OutTopology.h"
#include "../src/gcore/BondType.h"
#include "../src/gcore/Bond.h"
#include "../src/gcore/AngleType.h"
#include "../src/gcore/Angle.h"
#include "../src/gcore/DihedralType.h"
#include "../src/gcore/Dihedral.h"
#include "../src/gcore/ImproperType.h"
#include "../src/gcore/Improper.h"
#include "../src/gcore/LJType.h"
#include "../src/gcore/AtomPair.h"
#include "../src/gcore/Exclusion.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/GromosForceField.h"
#include "../src/gcore/Solvent.h"
#include "../src/gcore/SolventTopology.h"
#include "../src/gcore/System.h"
#include "../src/utils/AtomSpecifier.h"

using namespace gcore;
using namespace gio;
using namespace args;


int main(int argc, char *argv[]){

  char *knowns[] = {"topo", "mol"};
  int nknowns = 2;
  
  string usage = argv[0];
  usage += "\n\t@topo <topology>\n";
  usage += "\t@mol <mol>[:<first atom>-<last atom>]\n";

  try{
    Arguments args(argc, argv, nknowns, knowns, usage);

    InTopology it(args["topo"]);
    System sys(it.system());
    System syo;
    ostrstream addtitle;

    
    addtitle << "\nREDTOP molecules: ";
    
    // Parse atom specifiers, loop over all molecules to consider
    for(Arguments::const_iterator iter=args.lower_bound("mol"), 
          to=args.upper_bound("mol"); iter!=to;++iter) {
      string s=iter->second.c_str();
      utils::AtomSpecifier atspec(sys,s);
      
      int atoff=atspec.atom(0);
      int atend=atspec.atom(atspec.size()-1)+1;
      int mol=atspec.mol(0);
      

      MoleculeTopology mti=sys.mol(mol).topology();
      MoleculeTopology mto;

      AtomTopology ato;
  
      // First, do the atoms
      int lastres=-1, currres=-1;
      
      for(int i=atoff; i<atend; i++){
        Exclusion tempex, tempex14;
        
        ato.setName(mti.atom(i).name());
        ato.setIac(mti.atom(i).iac());
        ato.setChargeGroup(mti.atom(i).chargeGroup());
        ato.setMass(mti.atom(i).mass());
        ato.setCharge(mti.atom(i).charge());

        int num=mti.atom(i).exclusion().size();
        int num14=mti.atom(i).exclusion14().size();
        for(int k=0; k<num; ++k){
          int ati=mti.atom(i).exclusion().atom(k);
          if(ati<atend)
            tempex.insert(ati-atoff);
        }
        for(int k=0; k<num14; ++k){
          int ati14=mti.atom(i).exclusion14().atom(k);
          if(ati14<atend)
            tempex14.insert(ati14-atoff);
        }
      
        ato.setExclusion(tempex);
        ato.setExclusion14(tempex14);

        mto.addAtom(ato);
        if(mti.resNum(i)!=lastres){
	  lastres=mti.resNum(i);
	  currres++;
          mto.setResName(currres, mti.resName(mti.resNum(i)));
	}
	
        mto.setResNum(i-atoff,currres);
      }

      // Now, for the bonds, angles etc.
      BondIterator bi(mti);
      for(;bi;++bi){
        if(bi()[0]>=atoff&&bi()[0]<atend&&
            bi()[1]>=atoff&&bi()[1]<atend){
          Bond b(bi()[0]-atoff,bi()[1]-atoff);
          b.setType(bi().type());
          mto.addBond(b);
        }
      }
  
      AngleIterator ai(mti);
      for(;ai;++ai){
        if(ai()[0]>=atoff&&ai()[0]<atend&&
           ai()[1]>=atoff&&ai()[1]<atend&&
           ai()[2]>=atoff&&ai()[2]<atend){
          Angle a(ai()[0]-atoff,ai()[1]-atoff,ai()[2]-atoff);
          a.setType(ai().type());
          mto.addAngle(a);
        }
      }
      ImproperIterator ii(mti);
      for(;ii;++ii){
        if(ii()[0]>=atoff&&ii()[0]<atend&&
           ii()[1]>=atoff&&ii()[1]<atend&&
           ii()[2]>=atoff&&ii()[2]<atend&&
           ii()[3]>=atoff&&ii()[3]<atend){
          Improper imp(ii()[0]-atoff, ii()[1]-atoff,
                       ii()[2]-atoff, ii()[3]-atoff);
          imp.setType(ii().type());
          mto.addImproper(imp);
        }
      }
      DihedralIterator di(mti);
      for(;di;++di){
        if(di()[0]>=atoff&&di()[0]<atend&&
           di()[1]>=atoff&&di()[1]<atend&&
           di()[2]>=atoff&&di()[2]<atend&&
           di()[3]>=atoff&&di()[3]<atend){
          Dihedral dih(di()[0]-atoff, di()[1]-atoff,
                       di()[2]-atoff, di()[3]-atoff);
          dih.setType(di().type());
          mto.addDihedral(dih);
        }
      }
 
      syo.addMolecule(mto);
      addtitle << mol+1 << ":" << atoff+1 << "-" << atend << " ";
    }
    
    syo.addSolvent(sys.sol(0));
    
    
    OutTopology ot(cout);
    addtitle << "\nfrom: " << args["topo"] <<  '\0';
    ot.setTitle(it.title()+addtitle.str());

    ot.write(syo,it.forceField());
    
    return 0;
  }
  catch(gromos::Exception e){
    cerr << e.what() << endl;
    return 1;
  }
}





