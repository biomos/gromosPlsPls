/**
 * @file prepnoe.cc
 * Converts X-plor NOE data to GROMOS format
 */

/**
 * @page programs Program Documentation
 *
 * @anchor prepnoe
 * @section prepnoe Converts X-plor NOE data to GROMOS format
 * @author @ref mk co
 * @date 11-8-2006
 *
 * This program has been renamed to @ref prep_noe and will no longer be
 * maintained under this name
 *
 * <hr>
 */

#include <cassert>
#include <vector>
#include <map>
#include <iomanip>
#include <iostream>
#include <cmath>
#include <sstream>

#include <args/Arguments.h>
#include <gio/Ginstream.h>
#include <gio/InG96.h>
#include <gcore/System.h>
#include <gio/InTopology.h>
#include <gio/StringTokenizer.h>
#include <bound/Boundary.h>
#include <args/BoundaryParser.h>
#include <utils/VirtualAtom.h>
#include <utils/Neighbours.h>
#include "../src/gcore/AtomTopology.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/utils/AtomSpecifier.h"

using namespace gcore;
using namespace args;
using namespace gio;
using namespace bound;
using namespace std;
using namespace utils;

vector<VirtualAtom*> getvirtual(int atom, int type, int subtype, System &sys);

class Noeprep {
public:
  int residA;
  string atomA;
  int residB;
  string atomB;
  double dis;
  
  Noeprep(int resA, string A, int resB, string B, double d){
    residA=resA;
    atomA=A;
    residB=resB;
    atomB=B;
    dis=d;}
  ~Noeprep(){}
};

class Noelib{
public:
  string resname;
  string orgatomname;
  string gratomname;
  int NOETYPE;
  int NOESUBTYPE;
  
  Noelib(string A, string B, string C, string D){
    resname=A;
    orgatomname=B;
    gratomname=C;
    NOETYPE=atoi(D.c_str());
    NOESUBTYPE=0;
  }
  Noelib(string A, string B, string C, string D, string E){
    resname=A;
    orgatomname=B;
    gratomname=C;
    NOETYPE=atoi(D.c_str());
    NOESUBTYPE=atoi(E.c_str());
  }
  
  ~Noelib(){}
};

int main(int argc,char *argv[]){


  // Usage string

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo         <topology>\n";
  usage += "\t@title        <NOE title for output>\n";
  usage += "\t@noe          <NOE specification file>\n";
  usage += "\t@lib          <NOE specification library>\n";
  usage += "\t[@dish        <carbon-hydrogen distance; default: 0.1 nm>]\n";
  usage += "\t[@disc        <carbon-carbon distance; default: 0.153 nm>]\n";
  usage += "\t[@parsetype   <Upper bound parse type: 1, 2 or 3> ]\n";
  usage += "\t        Choices are:\n";
  usage += "\t        1: Upper bound == first number\n";
  usage += "\t        2: Upper bound == first + third number (most common, default)\n";
  usage += "\t        3: Upper bound == first - second number (commonly the lower bound)\n";
  usage += "\t[@correction  <correction file> [<correction type>] ]\n";
  usage += "\t[@action      <add> or <subtract> correction from upper bound; default: add ]\n";
  usage += "\t[@filter      <discard NOE's above a certain distance [nm]; default 10000 nm>]\n";
  usage += "\t[@factor      <conversion factor Ang to nm; default is 10>]\n";

  // known arguments...
  char *knowns[]={"topo", "title", "filter", "factor", "noe", "lib", 
		  "parsetype", "correction", "dish", "disc", "action"};
  int nknowns = 11;
    
  // prepare cout for formatted output
  cout.setf(ios::right, ios::adjustfield);
  cout.setf(ios::fixed, ios::floatfield);
  cout.precision(3);
    
  try{

    // Getting arguments and checking if everything is known.
    Arguments args(argc,argv,nknowns,knowns,usage);

    // read topology
    InTopology it(args["topo"]);
    System sys(it.system());

    //read in title
    string tit;
    {
      Arguments::const_iterator iter=args.lower_bound("title");
      for( ; iter!=args.upper_bound("title"); ++iter){
	tit+= iter->second + " ";
      }
    }

    //read in filter threshold
    double filt=10000;
    {
      Arguments::const_iterator iter=args.lower_bound("filter");
      if(iter!=args.upper_bound("filter")){
	filt=atof(iter->second.c_str());
      }
    }

    //read in conversion factor
    double conv=10.0;
    {
      Arguments::const_iterator iter=args.lower_bound("factor");
      if(iter!=args.upper_bound("factor")){
	conv=atof(iter->second.c_str());
      }
    }

    // Read in and create the NOE list
    Ginstream nf(args["noe"]);
    vector<string> buffer;
    nf.getblock(buffer);
    if(buffer[0]!="NOESPEC")
      throw gromos::Exception("main",
			      "NOESPEC file does not contain an NOESPEC block!");
    if(buffer[buffer.size()-1].find("END")!=0)
      throw gromos::Exception("prepnoe", "NOE file " + nf.name() +
			      " is corrupted. No END in NOESPEC"
			      " block. Got\n"
			      + buffer[buffer.size()-1]);

    // in noe all noes will be stored.
    vector<Noeprep> noevec;

    // a map with noe's that need to be connected for postnoe
    map<int, vector<int> > connections;
    
    int ptype=2;
    string line;
    for(unsigned int j=1; j< buffer.size()-1; j++){
      StringTokenizer tok(buffer[j]);
      vector<string> tokens = tok.tokenize();
      
      int a=atoi(tokens[0].c_str());
      int b=atoi(tokens[2].c_str());
      double d=atof(tokens[4].c_str());
      double e=atof(tokens[5].c_str());
      double f=atof(tokens[6].c_str());
      if(tokens.size()>7){
	unsigned int g=atoi(tokens[7].c_str());
	if(g!=j)
	  throw gromos::Exception("prepnoe",
				  "Numbering in NOESPEC file (7th column) is not correct");
	int h=atoi(tokens[8].c_str());
	vector<int> links(h-1);
	for(int ii=0; ii<h-1; ii++)
	  links[ii]=atoi(tokens[9+ii].c_str());
	connections[g-1]=links;
      }
      //check for parsetype
      try{
	args.check("parsetype");
	{
	  Arguments::const_iterator iter=args.lower_bound("parsetype");
	  if(iter!=args.upper_bound("parsetype")){
	    ptype=atoi(iter->second.c_str());
	  }
	}
      }
      catch(Arguments::Exception e){
	ptype=2;
      } 

      switch(ptype){
	case 1: d = d;
	  break;
	case 2: d = d+f;
	  break;
	case 3: d = d-e; 
	  break;
	default:
	  throw gromos::Exception("prepnoe", args["parsetype"] + 
				  " unknown. Known types are 1, 2 and 3");
      }     
      
      
      //put crap in vector
      //this selects only backbone-backbone NOE's
      //  if (tokens[1] == "HN" && tokens[3] == "HN"){
      noevec.push_back(Noeprep(a,tokens[1],b,tokens[3], d));
      // }
    }
    nf.close();

    // Read in and create the NOE library
    Ginstream nff(args["lib"]);
    buffer.clear();
    nff.getblock(buffer);
    
    if(buffer[0]!="NOELIB")
      throw gromos::Exception("main",
			      "NOELIB file does not contain an NOELIB block!");
    if(buffer[buffer.size()-1].find("END")!=0)
      throw gromos::Exception("prepnoe", "Library file " + nff.name() +
			      " is corrupted. No END in NOELIB"
			      " block. Got\n"
			      + buffer[buffer.size()-1]);

    vector<Noelib> noelib;
    for(unsigned int j=1; j< buffer.size()-1; j++){
      StringTokenizer tok(buffer[j]);
      vector<string> tokens = tok.tokenize();
      
      // check for number of tokens and store the data
      if(tokens.size()==4)
	noelib.push_back(Noelib(tokens[0],tokens[1],tokens[2],tokens[3]));
      else if(tokens.size()==5)
	noelib.push_back(Noelib(tokens[0],tokens[1],tokens[2],tokens[3],
				tokens[4]));
    }
    nff.close();
    
    //check for inconsistency in library
    for (int i=0; i < int (noelib.size()); ++i){
      Noelib A = noelib[i];
      for (int j=0; j < int (noelib.size()); ++j){
	Noelib B = noelib[j];
	if ((A.resname == B.resname) && 
	    (A.orgatomname == B.orgatomname) &&
	    (A.gratomname == B.gratomname) && 
	    (A.NOETYPE != B.NOETYPE || A.NOESUBTYPE != B.NOESUBTYPE)) {
	  std::stringstream sa;
	  sa << A.resname << " " 
	     << A.orgatomname << " " 
	     << A.gratomname << " " 
	     << A.NOETYPE;
	  if(A.NOESUBTYPE || B.NOESUBTYPE)
	    sa << " " << A.NOESUBTYPE;
	  sa << " !AND! "   
	     << B.resname << " " 
	     << B.orgatomname << " " 
	     << B.gratomname << " " 
	     << B.NOETYPE; 
	  if(A.NOESUBTYPE || B.NOESUBTYPE)
	    sa << " " << B.NOESUBTYPE;
	  sa << endl;
	  throw gromos::Exception("prepnoe ",  sa.str() +
				  " Inconsistent assigment of NOETYPE within library!");
	}       
      }
    }
    
    //check whether to add or subtract correction
    bool add = true;
    bool sub = false;
    if(args.count("action")>0){
      if(args["action"]=="add" || args["action"]=="ADD"){
	add= true;
	sub=false;
      }
      else if(args["action"]=="sub" || args["action"]=="SUB"){
	add=false;
	sub= true;
      }
      else
	throw gromos::Exception("prepnoe",
				"action type " + args["action"] + " not known!");
    }

    //read in the correction file if it exists
    map<int,double> pseudocorrectiondata;
    map<int,double> multiplicitydata;
    bool pseudocorrection = false;
    bool multiplicitycorrection = false;
    try{
      args.check("correction");
      Ginstream corf(args["correction"]);
      //get NOECORGROMOS block
      buffer.clear();
      corf.getblock(buffer);
      
      if(buffer[0]!="NOECORGROMOS")
	throw gromos::Exception("main",
				"NOE correction file does not contain the "
				"NOECORGROMOS block!");
      if(buffer[buffer.size()-1].find("END")!=0)
	throw gromos::Exception("prepnoe","Correction file " + corf.name() +
				" is corrupted. No END in NOECORGROMOS"
				" block. Got\n"
				+ buffer[buffer.size()-1]);
      
      for(unsigned int j=1; j< buffer.size()-1; j++){
        StringTokenizer tok(buffer[j]);
        vector<string> tokens = tok.tokenize();
	pseudocorrectiondata[atoi(tokens[0].c_str())] = atof(tokens[1].c_str());
      }
      
      
      //get MULTIPLICITY block
      buffer.clear();
      corf.getblock(buffer);
      
      if(buffer[0]!="MULTIPLICITY")
	throw gromos::Exception("main",
				"NOE correction file does not contain the"
				" MULTIPLICITY block!");
      if(buffer[buffer.size()-1].find("END")!=0)
	throw gromos::Exception("prepnoe","Correction file " + corf.name() +
				" is corrupted. No END in MULTIPLICITY"
				" block. Got\n"
				+ buffer[buffer.size()-1]);

      for(unsigned int j=1; j< buffer.size()-1; j++){
        StringTokenizer tok(buffer[j]);
        vector<string> tokens = tok.tokenize();
	multiplicitydata[atoi(tokens[0].c_str())] = atof(tokens[1].c_str());
      }
      corf.close();
      
      //determine the correction type
      Arguments::const_iterator it=args.lower_bound("correction");
      ++it;  

      string ctype;

      if(it == args.upper_bound("correction")) {
       pseudocorrection = true;
       multiplicitycorrection = true;
      }
      else ctype =  it->second;
   
      if (ctype == "pseudo") pseudocorrection = true;
      else if (ctype == "multiplicity") multiplicitycorrection = true;
      else if (ctype != "") 
	throw gromos::Exception("prepnoe",
				"Correction type " + ctype + " not known!" +
				"Use 'pseudo' or 'multiplicity' to apply " +
				"only one type of correction");

    }
    catch(Arguments::Exception e){
      cout << "# No correction file used!" << endl; 
    }

    //try for disc and dish
    double dish = 0.1;
    if(args.count("dish")>0) dish=atof(args["dish"].c_str());
    double disc = 0.153;
    if(args.count("disc")>0) disc=atof(args["disc"].c_str());
    
    //cout the title and that kind of stuff
    cout << "TITLE" << endl;
    cout << "NOE specification file for: " << tit << endl;
    cout << "END" << endl;
    cout << "DISRESSPEC" << endl;
    cout << "# DISH: carbon-hydrogen distance" << endl;
    cout << "# DISC: carbon-carbon distance" << endl;
    cout << "# DISH,DISC" << endl;
    cout << dish << " " << disc << endl;
    
    //open the filter file...
    ofstream filterfile; filterfile.open("noe.filter");
    filterfile << "TITLE" << endl;
    filterfile << "NOE filter file for: " << tit << endl;
    filterfile << "END" << endl;
    filterfile << "NOEFILTER" << endl;
    filterfile << "# noe"
	       << setw(4) << "mol"
	       << setw(11) << "residue"
	       << setw(5) << "atom"
	       << setw(5) << "atom"
	       << setw(4) << "mol"
	       << setw(11) << "residue"
	       << setw(5) << "atom"
	       << setw(5) << "atom"
	       << setw(6) << "r0"
	       << " filter noe" << endl;
    
    //here goes the crappy code
    string resnameA, resnameB;
    int molA=0, molB=0, resnumA=0, resnumB=0, mol=0, atA=0, atB=0, atNumA=0, atNumB=0, count=0;
    int atNOE = 0, totalnoecount = 1;
    double filterbound = 0;
    
    for (int i=0; i < int (noevec.size()); ++i) {
      
      vector<int> links;
      if(connections.count(i)) links=connections[i];
      
      atNOE = i;
      ostringstream atname;
      Noeprep NOE = noevec[i];
      //first find the residue-name corresponding to
      //your input-number in the topology
      mol=0;
      atNumA=0;
      while(NOE.residA > (atNumA+=sys.mol(mol).topology().numRes())){
        ++mol;
        if(mol > sys.numMolecules())
          throw gromos::Exception("prepnoe" ,
				  "Residue number too high in input line:\n");
      }
      atA=NOE.residA;
      atA-=atNumA-sys.mol(mol).topology().numRes();
      resnumA = (atA-1);
      molA = mol;
      resnameA = (sys.mol(mol).topology().resName(atA-1));
      
      mol=0;
      atNumB=0;
      while(NOE.residB > (atNumB+=sys.mol(mol).topology().numRes())){
        ++mol;
        if(mol > sys.numMolecules())
          throw gromos::Exception("prepnoe", + "Residue number too high in input line:\n");
      }
      atB=NOE.residB;
      atB-=atNumB-sys.mol(mol).topology().numRes();
      resnumB = (atB-1);
      molB = mol;
      resnameB = (sys.mol(mol).topology().resName(atB-1)); 
      
      //then map the correct gromos-topology atomname
      //to your input-atomname based on the residue-name
      //and the input-atomname
      bool foundA = false;
      vector<VirtualAtom*> vatomA;
      vector<VirtualAtom*> vatomB;
      
      int p=0;
      for (int k=0;k< int (noelib.size());++k){
	Noelib NOELIB = noelib[k];
	if (NOELIB.resname==resnameA && NOELIB.orgatomname==NOE.atomA){
	  //back to topology to get the atom number
	  for(int f=0;f<sys.mol(molA).numAtoms() && foundA==false;++f){
	    if (sys.mol(molA).topology().atom(f).name() == NOELIB.gratomname &&
		sys.mol(molA).topology().resNum(f) == resnumA) {
	      int addA=0; 
	      foundA = true;
	      p=k;
	      
	      for(int i=0;i < molA;++i)	addA+=sys.mol(i).numAtoms();
	      
	      if (NOE.dis/conv > filt)	cout << "#";
              vatomA = getvirtual(f+addA, NOELIB.NOETYPE, NOELIB.NOESUBTYPE, sys);
	    }
	  }
	}
      }
      
      if (!foundA) {
	std::stringstream ss;
	string a;
	ss << NOE.residA;
	ss >> a; 
	string b = NOE.atomA;
	string c = " ";
	string d = a+c+b;
	throw gromos::Exception("prepnoe ",  d +
				" Noe specification not found in library!");
      }   
      Noelib NA = noelib[p];
      
      
      bool foundB = false; 
      for (int z=0;z< int (noelib.size());++z){  
	Noelib NOELIBB = noelib[z];    
	if (NOELIBB.resname==resnameB && NOELIBB.orgatomname==NOE.atomB){
	  //back to topology to get the atom number
	  for(int g=0;g<sys.mol(molB).numAtoms() && foundB==false;++g){
	    if (sys.mol(molB).topology().atom(g).name() == NOELIBB.gratomname &&
		sys.mol(molB).topology().resNum(g) == resnumB){
	      int addB=0;
	      foundB=true;
	      count+=1;
	      for(int i=0;i < molB;++i)	addB+=sys.mol(i).numAtoms();
              vatomB = getvirtual(g+addB, NOELIBB.NOETYPE, NOELIBB.NOESUBTYPE, sys);
	      
              atname << setw(3) << molA+1 << " "
                           << setw(5) << resnumA+1 << " "
			   << setw(4) << NA.resname << " " 
			   << setw(4) << NA.gratomname << " " 
			   << setw(4) << NA.orgatomname << " "
			   << setw(3) << molB+1 << " " 
                           << setw(5) << (resnumB+1) << " " 
			   << setw(4) << NOELIBB.resname << " " 
			   << setw(4) << NOELIBB.gratomname << " " 
                           << setw(4) << NOELIBB.orgatomname << " ";
	      
	    }
	  }
	}
      }
      
      if (!foundB) {
	std::stringstream s;
	string aa;
	s << NOE.residB;
	s >> aa; 
	string bb = NOE.atomB;
	string cc = " ";
	string dd = aa+cc+bb;
	throw gromos::Exception("prepnoe ",  dd +
				" Noe specification not found in library!");
      }
      
      //spit out disresblock...
      int atomsA[4], atomsB[4];
      int creatednoe = 0;
      
      for (int va=0; va < (int) vatomA.size(); ++va) {
	int offsetA = 1;
	VirtualAtom VA(*vatomA[va]);
	
	int mol = VA.conf().mol(0);
	for(int l=0;l<mol;++l) offsetA += sys.mol(l).numAtoms();
	
	for (int aa=0; aa < 4; ++aa) {
	  int att;
	  if(VA.conf().size() > aa)
	    att = VA.conf().atom(aa);
	  else
	    att=-1;
	  
	  atomsA[aa] = att;
	}
	
	for (int vb=0; vb < (int) vatomB.size(); ++vb) {
          int offsetB = 1;
	  VirtualAtom VB(*vatomB[vb]);
	  int mol = VB.conf().mol(0);                
	  for(int l=0;l<mol;++l) offsetB += sys.mol(l).numAtoms();
	  
	  for (int bb=0; bb < 4; ++bb) {
	    int att;
	    if(VB.conf().size()>bb)
	      att = VB.conf().atom(bb);
	    else 
	      att = -1;
	    

	    atomsB[bb] = att;
	  }
	  
	  ostringstream ss;
	  ss.setf(ios::right, ios::adjustfield);
	  ss.setf(ios::fixed, ios::floatfield);
	  ss.precision(3);
	  
	  
	  for (int kk=0; kk < 4; ++kk) {
	    if(atomsA[kk] == -1) ss << setw(5) << 0;
	    else ss << setw(5) << atomsA[kk] + offsetA;
	  }
	  ss << setw(3) << VA.type();
	  
	  for (int kk=0; kk < 4; ++kk) {
	    if(atomsB[kk] == -1) ss << setw(5) << 0;
	    else ss << setw(5) << atomsB[kk] + offsetB;
	  }
	  ss << setw(3) << VB.type();
	  
	  
	  
	  double bound = NOE.dis/conv;
	  //set up corrections
	  vector<int> type; 
	  type.push_back(VA.type());
	  type.push_back(VB.type());
	  
	  // first do the multiplicity correction and then the 
	  // pseudo atom correction
	  
	  double mult=1;
	  double cor=0;
	  
	  //check for multiplicity-correction              
	  if (multiplicitycorrection) {
	    for (int i=0; i < (int) type.size(); ++i) {
	      std::map<int,double>::iterator k = multiplicitydata.find(type[i]);
	      if (k != multiplicitydata.end()) mult *= k->second; 
	    }                
	  } //end if (multiplicitycorrection) 
	  
	  
	  //check for gromos-correction
	  if (pseudocorrection) {
	    for (int i=0; i < (int) type.size(); ++i) {
	      std::map<int,double>::iterator k = pseudocorrectiondata.find(type[i]);
	      if (k != pseudocorrectiondata.end()) cor += k->second;
	    }	       
	  } //end if (pseudocorrection) 
	  if(add) bound = bound*mult + cor;
	  else if(sub) bound = (bound - cor) / mult;
	  
	  // in the filter file I also want the corrected bound
	  filterbound = bound;
	  
	  cout << ss.str() << "   "   << bound << " 1.0000" << endl;
	  if (va > 0 || vb > 0) 	cout << "#automatically generated NOE distance to monitor!" << endl;
	  
	  ++creatednoe;
	  
	}
      }
      cout << "# " << count << atname.str() << endl;
      
      //smack out the filterfile...
      for (int ii=0; ii < creatednoe; ++ii) {
	filterfile << setw(5) << totalnoecount << " ";
	filterfile << atname.str();
	filterfile << setw(5) << filterbound << " ";
	int offset = totalnoecount-i+creatednoe-ii-2;
	filterfile << " " << creatednoe + links.size() << " ";             	       
	for (int iii=0; iii < creatednoe; ++iii) {
	  //                if (ii != iii) filterfile << " " << atNOE+1+iii;
	  if(ii!=iii) filterfile << " " << totalnoecount - ii + iii;
	}
	for(unsigned int iii=0; iii < links.size(); ++iii){
	  int iiii=links.size()-1-iii;
	  filterfile << " " << offset+links[iiii];
	}
	
	filterfile << endl;
	++totalnoecount;
      }
      
      
      
    } //end for (int i=0; i < noevec.size()) ++i) ...
    
    cout << "END" << endl;
    filterfile << "END" << endl;
    filterfile.close();           
    
  }
  
  catch(gromos::Exception e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}

vector<VirtualAtom*> getvirtual(int at, int type, int subtype, System &sys) {

  int mol=0, atNum=0;

  // parse into mol and atom rather than high atom nr.
   while(at >= (atNum+=sys.mol(mol).numAtoms())){
    ++mol;
    if(mol >= sys.numMolecules())
     throw gromos::Exception("prepnoe ", +"Atom number too high in input atom number:\n"+at);
   }
    at-=atNum-sys.mol(mol).numAtoms();

    vector<VirtualAtom*> vat;
    
    if (type == 4 ){
      if(subtype==0){
	
	//we need to automatically generate the r, l atoms...
	vat.push_back(new VirtualAtom(sys,mol,at, VirtualAtom::virtual_type(type)));
	vat.push_back(new VirtualAtom(sys,mol,at, VirtualAtom::virtual_type(type),1));
      }
      else if(subtype==1){
	vat.push_back(new VirtualAtom(sys,mol,at,VirtualAtom::virtual_type(type)));
      }
      else if(subtype==2){
	vat.push_back(new VirtualAtom(sys,mol,at,VirtualAtom::virtual_type(type),1));
      }
    }
    else vat.push_back(new VirtualAtom(sys,mol,at, VirtualAtom::virtual_type(type)));     
    
    return vat;

}
