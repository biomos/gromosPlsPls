/**
 * @file renumber.cc
 * re-assigns force field types in building blocks
 */

/**
 * @page contrib Contrib Program Documentation
 *
 * @anchor renumber
 * @section renumber re-assigns force field types in building blocks
 * @author @ref co
 * @date 31-10-08
 *
 * Program renumber can help to convert building blocks from one force field definition to the other. It uses a renumber file, that lists for all bond-, angle-, dihedral-, improper- and atomtypes the original type number and the new type number. Note that it only translate the types for a building block, but does not introduce new types or modify charge distributions.
 *
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@build</td><td>&lt;mtb-file&gt; </td></tr>
 * <tr><td> \@block</td><td>&lt;buildingblock name&gt; </td></tr>
 * <tr><td> \@renum</td><td>&lt;renumber file&gt; </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  renumber
    @build   mtb45a4.dat
    @block   ALA
    @renum   ren_45a4_to_53a6.dat
 @endverbatim
 *
 * <hr>
 */
#include <cassert>
#include <iostream>
#include <iomanip>
#include <map>
#include "../src/args/Arguments.h"
#include "../src/gio/InBuildingBlock.h"
#include "../src/gcore/BuildingBlock.h"
#include "../src/gcore/BbSolute.h"


#include "../src/gcore/Bond.h"
#include "../src/gcore/Angle.h"
#include "../src/gcore/Dihedral.h"
#include "../src/gcore/Improper.h"
#include "../src/gcore/Exclusion.h"
#include "../src/gcore/AtomTopology.h"

#include "../src/gio/Ginstream.h"

using namespace gcore;
using namespace gio;
using namespace args;
using namespace std;


int main(int argc, char *argv[]){
  Argument_List knowns;
  knowns << "build" << "block" << "renum";
  
  string usage = "# " + string(argv[0]);
  usage += "\n\t@build <mtb-file>\n";
  usage += "\t@block <buildingblock name>\n";
  usage += "\t@renum <renumber file>\n";

  try{
    Arguments args(argc, argv, knowns, usage);

    // read in the building block file
    InBuildingBlock ibb(args["build"]);
    BuildingBlock mtb(ibb.building());
    
    // define the selected block
    int index=mtb.findBb(args["block"]);
    if(index==0)
      throw gromos::Exception("renumber",
	      "Building block " + args["block"] + " not found" );
    BbSolute bb;
    int endgroup=0;
    
    if(index>0)
      bb=mtb.bb(index-1);
    else{
      
      bb=mtb.be(-index-1);
      endgroup=1;
    }
    
    // read in the renumber file
    Ginstream gin(args["renum"]);

    // quite ugly to define this here, but who cares
    const int max_number_of_types=100;
    
    map<int, int> atomtypes, bondtypes, angletypes, impropertypes, 
      dihedraltypes;
    std::vector<std::string> buffer;
    std::vector<std::vector<std::string > > content;
    while(!gin.stream().eof()){
      gin.getblock(buffer);
      if(!gin.stream().eof()){ 
	if(buffer[buffer.size()-1].find("END")!=0)
	  throw gromos::Exception("renumber", "Renumber file " + gin.name() +
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
      if      ((*iter)[0]=="ATOMTYPE"     || (*iter)[0]=="ATOMTYPECONV")
	pointer_to_a_map = &atomtypes;
      else if ((*iter)[0]=="BONDTYPE"     || (*iter)[0]=="BONDTYPECONV")
	pointer_to_a_map = &bondtypes;
      else if ((*iter)[0]=="ANGLETYPE"    || (*iter)[0]=="ANGLETYPECONV")
	pointer_to_a_map = &angletypes;
      else if ((*iter)[0]=="IMPROPERTYPE" || (*iter)[0]=="IMPROPERTYPECONV")
	pointer_to_a_map = &impropertypes;
      else if ((*iter)[0]=="DIHEDRALTYPE" || (*iter)[0]=="DIHEDRALTYPECONV") 
	pointer_to_a_map = &dihedraltypes;
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
    for(int i=1; i< max_number_of_types; i++) {
      if(atomtypes[i]==0)     atomtypes[i]=i;
      if(bondtypes[i]==0)     bondtypes[i]=i;
      if(angletypes[i]==0)    angletypes[i]=i;
      if(impropertypes[i]==0) impropertypes[i]=i;
      if(dihedraltypes[i]==0) dihedraltypes[i]=i;
    }

    int numBonds=0;
    int numAngles=0;
    int numImpropers=0;
    int numDihedrals=0;
    

    //Now we need the output, but there is no outBuildingBlock yet!
    cout.precision(5);
    int last_few=0;
    
    if(endgroup){
      cout << "MTBUILDBLEND" << endl;
      last_few=bb.rep();
    }
    else{
      cout << "MTBUILDBLSOLUTE" << endl;
      last_few=bb.numPexcl();
    }
    cout << "# building block (residue, nucleotide, etc.)" << endl;
    cout << "# RNME" << endl;
    cout << bb.resName() << endl;
    if(endgroup){
      cout << "# number of atoms, number of atoms to be replaced" << endl;
      cout << "# NMAT,NREP" << endl;
    }
    else{
      cout << "# number of atoms, number of preceding exclusions" << endl;
      cout << "# NMAT,NLIN" << endl;
    }
    cout << setw(5) << bb.numAtoms();
    if(endgroup)
      cout << setw(5) << bb.rep() << endl;
    else{
      cout << setw(5) << bb.numPexcl() << endl;
      cout << "# preceding exclusions" << endl;
      cout << "#ATOM                               MAE MSAE" << endl;
      for(int i=0; i< bb.numPexcl(); i++){
	cout << setw(5) << i+1-bb.numPexcl()
	     << setw(34) << bb.pexcl(i).size();
	for(int j=0; j< bb.pexcl(i).size();j++)
	  cout << setw(5) << bb.pexcl(i).atom(j)+1;
	cout << endl;
      }
    }
    
    cout << "# atoms" << endl;
    cout << "#ATOM ANM  IACM MASS        CGMICGM MAE MSAE" << endl;
    for(int i=0; i<bb.numAtoms(); i++){
      if(i== bb.numAtoms() - last_few)
	if(endgroup)
	  cout << "# replacing atoms" << endl;
	else 
	  cout << "# trailing atoms" << endl
	       << "#ATOM ANM  IACM MASS        CGMICGM" << endl;
      cout << setw(5) << i+1 << ' ';

      cout.setf(ios::left, ios::adjustfield);
      
      cout << setw(4) << bb.atom(i).name();
      cout.setf(ios::fixed, ios::adjustfield);
      cout.precision(5);
      cout.setf(ios::fixed, ios::floatfield);
      
      cout << setw(5) << atomtypes[bb.atom(i).iac()+1]
	   << setw(5) << int(bb.atom(i).mass())+1
	   << setw(11) << bb.atom(i).charge()
	   << setw(4) << bb.atom(i).chargeGroup();
      
      if(i < bb.numAtoms() - last_few){
	cout << setw(4) << bb.atom(i).exclusion().size();
	for(int j=0; j< bb.atom(i).exclusion().size(); j++){
	  cout << setw(5) << bb.atom(i).exclusion().atom(j)+1;
	  if((j+1)%6==0 && j+1 < bb.atom(i).exclusion().size()) 
	    cout << endl << setw(39) << " ";
	}
      }
      cout << endl;
    }
    cout << "# bonds" << endl;
    cout << "#  NB" << endl;
    {
      BondIterator bi(bb);
      for(; bi; ++bi) numBonds++;
    }
    cout << setw(5) << numBonds << endl;
    cout << "#  IB   JB  MCB" << endl;
    BondIterator bi(bb);
    for(;bi;++bi)
      cout << setw(5) << bi()[0]+1
	   << setw(5) << bi()[1]+1
	   << setw(5) << bondtypes[bi().type()+1] << endl;
    cout << "# bond angles" << endl;
    cout << "# NBA" << endl;
    {
      AngleIterator ai(bb);
      for(;ai;++ai) numAngles++;
    }
    cout << setw(5) << numAngles << endl;
    cout << "#  IB   JB   KB  MCB" << endl;

    AngleIterator ai(bb);
    for(;ai;++ai)
      cout << setw(5) << ai()[0]+1
	   << setw(5) << ai()[1]+1
	   << setw(5) << ai()[2]+1
	   << setw(5) << angletypes[ai().type()+1] << endl;
    cout << "# improper dihedrals" << endl;
    cout << "# NIDA" << endl;
    {
      ImproperIterator ii(bb);
      for(;ii;++ii) numImpropers++;
    }
    
    cout << setw(5) << numImpropers << endl;
    cout << "#  IB   JB   KB   LB  MCB" << endl;
    ImproperIterator ii(bb);
    for(;ii;++ii)
      cout <<  setw(5) << ii()[0]+1
	   <<  setw(5) << ii()[1]+1
	   <<  setw(5) << ii()[2]+1
	   <<  setw(5) << ii()[3]+1
	   <<  setw(5) << impropertypes[ii().type()+1] << endl;
    cout << "# dihedrals" << endl;
    cout << "# NDA" << endl;
    {
      DihedralIterator di(bb);
      for(;di; ++di) numDihedrals++;
    }
    
    cout << setw(5) << numDihedrals << endl;
    cout << "#  IB   JB   KB   LB  MCB" << endl;
    DihedralIterator di(bb);
    for(;di;++di)
      cout <<  setw(5) << di()[0]+1
	   <<  setw(5) << di()[1]+1
	   <<  setw(5) << di()[2]+1
	   <<  setw(5) << di()[3]+1
	   <<  setw(5) << dihedraltypes[di().type()+1] << endl;   
    cout << "END" << endl;
    
    

  }
  catch(gromos::Exception e){
    cerr << e.what() << endl;
    return 1;
  }
}





