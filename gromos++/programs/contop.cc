#include <iostream>
#include <cassert>

#include "../src/args/Arguments.h"
#include "../src/gio/InTopology.h"
#include "../src/gio/InParameter.h"
#include "../src/gio/OutTopology.h"
#include "../src/gcore/GromosForceField.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/Bond.h"
#include "../src/gcore/Angle.h"
#include "../src/gcore/Dihedral.h"
#include "../src/gcore/Improper.h"
#include "../src/gcore/Exclusion.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gcore/Solvent.h"
#include "../src/gcore/SolventTopology.h"
#include "../src/gio/Ginstream.h"

using namespace std;
using namespace gcore;
using namespace gio;
using namespace args;

void renumber_types(System &sys, GromosForceField &gff, string renum);
void check_types(System &sys, GromosForceField &gff);

int main(int argc, char *argv[]){

  char *knowns[] = {"topo", "param", "renum"};
  int nknowns = 3;
  
  string usage = argv[0];
  usage += "\n\t@topo <topology>\n";
  usage += "\t@param <gromos parameter file/topo>\n";
  usage += "\t@renum <renumber file>\n";
  
  try{
    Arguments args(argc, argv, nknowns, knowns, usage);

    InTopology it(args["topo"]);

    System sys(it.system());
    
    OutTopology ot(cout);
    string addtitle;
    addtitle+="CONTOP parameters: "+args["param"];
    if (args.count("renum")>0) 
      addtitle += "\nrenumbering from: " + args["renum"];
    
    ot.setTitle(it.title()+addtitle);

    InParameter ip(args["param"]);
    GromosForceField gff(ip.forceField());
    gff.setFpepsi(it.forceField().fpepsi());
    gff.setHbar(it.forceField().hbar());

    // maybe some types are to be renumbered?
    if(args.count("renum")>0)
      renumber_types(sys, gff, args["renum"]);

    // check if the topology is now not referring to non-existing types
    check_types(sys, gff);
    
    ot.write(sys,gff);
    
  }
  catch(gromos::Exception e){
    cerr << e.what() << endl;
    return 1;
  }
}

void check_types(System &sys, GromosForceField &gff)
{
  int num_atom_types=gff.numAtomTypeNames();
  int num_bond_types=gff.numBondTypes();
  int num_angle_types=gff.numAngleTypes();
  int num_improper_types=gff.numImproperTypes();
  int num_dihedral_types=gff.numDihedralTypes();
  
  for(int m=0; m< sys.numMolecules();m++){
    for(int a=0; a<sys.mol(m).numAtoms(); a++){
      if(sys.mol(m).topology().atom(a).iac() >= num_atom_types){
	ostringstream os;
	os << "Atom " << m+1 << ":" << a+1 << " has a higher integer "
	   << "atom code (" << sys.mol(m).topology().atom(a).iac()+1
	   << ") than defined in the parameter file ("
	   << num_atom_types << ")";
	
	throw gromos::Exception("contop", os.str());
      }
    }
    BondIterator bi(sys.mol(m).topology());
    for(;bi;++bi){
      if(bi().type() >= num_bond_types){
	ostringstream os;
	os << "Bond " << m+1 << ":" << bi()[0] << " - " << m+1 << ":" 
	   << bi()[1] << " has a higher bondtype (" << bi().type()+1 
	   << ") than defined in the parameter file ("
	   << num_bond_types << ")";
	throw gromos::Exception("contop", os.str());
      }
    }
    AngleIterator ai(sys.mol(m).topology());
    for(;ai;++ai){
      if(ai().type() >= num_angle_types){
	ostringstream os;
	os << "Angle " << m+1 << ":" << ai()[0]+1 << " - " << m+1 << ":" 
	   << ai()[1]+1 << " - " << m+1 << ":" << ai()[2]+1
	   << " has a higher angletype (" << ai().type()+1 
	   << ") than defined in the parameter file ("
	   << num_angle_types << ")";
	throw gromos::Exception("contop", os.str());
      }
    }
    ImproperIterator ii(sys.mol(m).topology());
    for(;ii;++ii){
      if(ii().type() >= num_improper_types){
	ostringstream os;
	os << "Improper " << m+1 << ":" << ii()[0]+1 << " - " << m+1 << ":" 
	   << ii()[1]+1 << " - " << m+1 << ":" << ii()[2]+1 << " - " 
	   << m+1 << ":" << ii()[3]+1 
	   << " has a higher impropertype (" << ii().type()+1 
	   << ") than defined in the parameter file ("
	   << num_improper_types << ")";
	throw gromos::Exception("contop", os.str());
      }
    }    
    DihedralIterator di(sys.mol(m).topology());
    for(;di;++di){
      if(di().type() >= num_dihedral_types){
	ostringstream os;
	os << "Improper " << m+1 << ":" << di()[0]+1 << " - " << m+1 << ":" 
	   << di()[1]+1 << " - " << m+1 << ":" << di()[2]+1 << " - " 
	   << m+1 << ":" << di()[3]+1 
	   << " has a higher dihedraltype (" << di().type()+1 
	   << ") than defined in the parameter file ("
	   << num_improper_types << ")";
	throw gromos::Exception("contop", os.str());
      }
    }
  }
}

void renumber_types(System &sys, GromosForceField &gff, string renum)
{
  // read in the renumber file
  Ginstream gin(renum);

  int num_atom_types=gff.numAtomTypeNames();
  int num_bond_types=gff.numBondTypes();
  int num_angle_types=gff.numAngleTypes();
  int num_improper_types=gff.numImproperTypes();
  int num_dihedral_types=gff.numDihedralTypes();
    
  map<int, int> atomtypes, bondtypes, angletypes, impropertypes, 
    dihedraltypes;
  std::vector<std::string> buffer;
  std::vector<std::vector<std::string > > content;
  while(!gin.stream().eof()){
    gin.getblock(buffer);
    if(!gin.stream().eof()){
      if(buffer[buffer.size()-1].find("END")!=0)
	throw gromos::Exception("contop", "Renumber file " + gin.name() +
				" is corrupted. No END in "+buffer[0]+
				" block. Got\n"
				+ buffer[buffer.size()-1]);
      content.push_back(buffer);
    }    
  }
  // now loop over the content
  std::vector<std::vector<std::string > >::const_iterator 
    iter=content.begin();
  for( ; iter!=content.end(); ++iter){
    map<int, int> *pointer_to_a_map=NULL;
    if ((*iter)[0]=="ATOMTYPE")          pointer_to_a_map = &atomtypes;
    else if ((*iter)[0]=="BONDTYPE") 	   pointer_to_a_map = &bondtypes;
    else if ((*iter)[0]=="ANGLETYPE")    pointer_to_a_map = &angletypes;
    else if ((*iter)[0]=="IMPROPERTYPE") pointer_to_a_map = &impropertypes;
    else if ((*iter)[0]=="DIHEDRALTYPE") pointer_to_a_map = &dihedraltypes;
    else
      throw gromos::Exception("renumber", 
			      "Don't know how to handle "+(*iter)[0] + "-block");
    
    int a, b;
    
    // now loop over the contents of the block
    for(unsigned int i=1; i< (*iter).size()-1; i++){
      std::istringstream linestream((*iter)[i]);
      linestream >> a >> b;
      (*pointer_to_a_map)[a]=b;
    }
  }
  // let's fill up all types that are not used with themselves
  
  //atomtypes
  for(int i=1; i< num_atom_types; i++)
    if(atomtypes[i]==0) atomtypes[i]=i;
  for(int i=1; i< num_bond_types; i++) 
    if(bondtypes[i]==0) bondtypes[i]=i;
  for(int i=1; i< num_angle_types; i++) 
    if(angletypes[i]==0) angletypes[i]=i;
  for(int i=1; i< num_improper_types; i++) 
    if(impropertypes[i]==0) impropertypes[i]=i;
  for(int i=1; i< num_dihedral_types; i++) 
    if(dihedraltypes[i]==0) dihedraltypes[i]=i;
  
  // Now loop over all the bonds, angles and atoms in the topology to
  // replace the types

  // molecules
  for(int m=0; m<sys.numMolecules(); m++){

    MoleculeTopology mt;
    
    // atoms
    for(int a=0; a < sys.mol(m).numAtoms(); a++){
      
      sys.mol(m).topology().atom(a).setIac( 
	atomtypes[sys.mol(m).topology().atom(a).iac()+1] - 1);
      mt.addAtom(sys.mol(m).topology().atom(a));
      mt.setResNum(a,sys.mol(m).topology().resNum(a));
    }
    for(int r=0; r < sys.mol(m).topology().numRes(); r++)
      mt.setResName(r, sys.mol(m).topology().resName(r));
    
    // bonds
    BondIterator bi(sys.mol(m).topology());
    for(;bi;++bi){
      Bond b=bi();
      b.setType(bondtypes[bi().type()+1] - 1);
      mt.addBond(b);
    }
    
    // Angles
    AngleIterator ai(sys.mol(m).topology());
    for(;ai; ++ai){
      Angle a=ai();
      a.setType(angletypes[ai().type()+1] - 1);
      mt.addAngle(a);
    }
    
    // Impropers
    ImproperIterator ii(sys.mol(m).topology());
    for(;ii; ++ii){
      Improper i=ii();
      i.setType(impropertypes[ii().type()+1] - 1);
      mt.addImproper(i);
    }
    
    // Dihedrals
    DihedralIterator di(sys.mol(m).topology());
    for(;di; ++di){
      Dihedral d=di();
      d.setType(dihedraltypes[di().type()+1] - 1);
      mt.addDihedral(d);
    }
    sys.mol(m).topology() = mt;

    // after possibly renumbering atom types, or new masses
    // we have to assign hydrogen atoms new.
    sys.mol(m).topology().clearH();
    // determine H atoms based on the masses
    sys.mol(m).topology().setHmass(1.008);
    // or on the iac-code
    // sys.mol(m).topology().setHmass(17);
  }
  // And don't forget the solvent!! Okay, I did forget the solven.
  for(int s=0; s<sys.numSolvents(); s++){
    for(int a=0; a<sys.sol(s).topology().numAtoms(); a++){
      sys.sol(s).topology().atom(a).setIac(
	   atomtypes[sys.sol(s).topology().atom(a).iac()+1] - 1);
    }
  }
  
}








