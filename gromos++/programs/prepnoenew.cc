// noeprep.cc; 
#include <cassert>

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

#include <vector>
#include <map>
#include <iomanip>
#include <iostream>
#include <cmath>
#include <sstream>

using namespace gcore;
using namespace args;
using namespace gio;
using namespace bound;
using namespace std;
using namespace utils;

vector<VirtualAtom*> getvirtual(int atom, int type, System &sys);


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

  Noelib(string A, string B, string C, string D){
    resname=A;
    orgatomname=B;
    gratomname=C;
    NOETYPE=atoi(D.c_str());
  }
  ~Noelib(){}
};

int main(int argc,char *argv[]){


  // Usage string

  string usage = argv[0];
  usage += "\n\t@topo <topology>\n";
  usage += "\t@title <NOE title for output>\n";
  usage += "\t@filter <discard NOE's above a certain distance [nm]; default 10000 nm, so you will omit none>\n";
  usage += "\t@factor <conversion factor Ang <-> nm; default is 10>\n";
  usage += "\t@noe <NOE specification file>\n";
  usage += "\t@lib <NOE specification library>\n";
  usage += "\t[@dish <carbon-hydrogen distance; default: 0.1 nm>\n";
  usage += "\t[@disc <carbon-carbon distance; default: 0.153 nm>\n";
  usage += "\t[@parsetype [Upper bound parse type <1 2 3>] ]\n";
  usage += "            Choices are:\n";
  usage += "             1: Upper bound == first number\n";
  usage += "             2: Upper bound == first + third number (most common, default)\n";
  usage += "             3: Upper bound == first - second number (commonly the lower bound)\n";
  usage += "\t[ @correction [correction file] ]\n";
  usage += "\t[ @addorsubtractcorrection [add or substract correction from upper bound; default: add] ]\n";

  // defining all sorts of constants
  // known arguments...
  char *knowns[]={"topo", "title", "filter", "factor", "noe", "lib", "parsetype", "correction", "dish", "disc", "addorsubtractcorrection"};
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
      if(iter!=args.upper_bound("title")){
	tit=(iter->second.c_str());
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

    // in noe all noes will be stored.
    vector<Noeprep> noevec;
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

    // in noe all noes will be stored.
    vector<Noelib> noelib;
    for(unsigned int j=1; j< buffer.size()-1; j++){
      StringTokenizer tok(buffer[j]);
      vector<string> tokens = tok.tokenize();
      
      //put crap in vector
      noelib.push_back(Noelib(tokens[0],tokens[1],tokens[2], tokens[3]));
    }
    nff.close();

    //check for inconsistency in library
    for (int i=0; i < int (noelib.size()); ++i){
      Noelib A = noelib[i];
      for (int j=0; j < int (noelib.size()); ++j){
	Noelib B = noelib[j];
	if ((A.resname == B.resname) && (A.orgatomname == B.orgatomname) 
	    && (A.gratomname == B.gratomname) && (A.NOETYPE != B.NOETYPE)) {
	  std::stringstream sa;
	  sa << A.resname << " " << A.orgatomname << " " << A.gratomname << " " << A.NOETYPE
	     << " !AND! "   
	     << B.resname << " " << B.orgatomname << " " << B.gratomname << " " << B.NOETYPE << endl << "\0";
	  string ex;
	  ex = sa.str();
   
	  throw gromos::Exception("prepnoe ",  ex +
				  " Inconsistent assigment of NOETYPE within library!");
	}       
      }
    }

    //check whether to add or substract correction
    bool add = false;
    bool sub = false;

    try{
      args.check("addorsubtractcorrection");      
      //get NOECORGROMOS block
       Arguments::const_iterator iter=args.lower_bound("addorsubtractcorrection");
       // ++it;  

       string ctype;
      if(iter!=args.upper_bound("addorsubtractcorrection")){
     
	
       ctype = iter->second;
   
      if (ctype == "add") add = true;
      else if (ctype == "sub") sub = true;
      else if (ctype != "") throw gromos::Exception("prepnoe",
		     "addorsubtractcorrection type " + ctype + " not known!");

      }
    }
    catch(Arguments::Exception e){
      add = true; 
    }
      

    //read in the correction file if it exists
    map<int,double> gromoscorrectiondata;
    map<int,double> multiplicitydata;
    bool gromoscorrection = false;
    bool multiplicitycorrection = false;
    try{
      args.check("correction");      
      Ginstream corf(args["correction"]);
      //get NOECORGROMOS block
      buffer.clear();
      corf.getblock(buffer);
      
      if(buffer[0]!="NOECORGROMOS")
	throw gromos::Exception("main",
		     "NOE correction file does not contain the NOECORGROMOS block!");
      for(unsigned int j=1; j< buffer.size()-1; j++){
        StringTokenizer tok(buffer[j]);
        vector<string> tokens = tok.tokenize();
	gromoscorrectiondata[atoi(tokens[0].c_str())] = atof(tokens[1].c_str());
      }
      

      //get MULTIPLICITY block
      buffer.clear();
      corf.getblock(buffer);
      
      if(buffer[0]!="MULTIPLICITY")
	throw gromos::Exception("main",
		     "NOE correction file does not contain the MULTIPLICITY block!");
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
       gromoscorrection = true;
       multiplicitycorrection = true;
      }
      else ctype =  it->second;
   
      if (ctype == "gromos") gromoscorrection = true;
      else if (ctype == "multiplicity") multiplicitycorrection = true;
      else if (ctype != "") throw gromos::Exception("prepnoe",
		     "Correction type " + ctype + " not known!");

    }
    catch(Arguments::Exception e){
      cout << "# No correction file used!" << endl; 
    }

    //try for disc and dish
    double dish = 0.1;

    try{
	args.check("dish");
	{
	  Arguments::const_iterator iter=args.lower_bound("dish");
	  if(iter!=args.upper_bound("dish")){
	    dish=atof(iter->second.c_str());
	  }
	}
      }
      catch(Arguments::Exception e){
	dish= 0.1;
      }

    double disc = 0.153;

    try{
	args.check("disc");
	{
	  Arguments::const_iterator iter=args.lower_bound("disc");
	  if(iter!=args.upper_bound("disc")){
	    disc=atof(iter->second.c_str());
	  }
	}
      }
      catch(Arguments::Exception e){
	disc= 0.153;
      }
    
    //cout the title and that kind of stuff
      cout << "TITLE" << endl;
      cout << "NOE specification file for " << tit << endl;
      cout << "END" << endl;
      cout << "DISRESSPEC" << endl;
      cout << "# DISH: carbon-hydrogen distance" << endl;
      cout << "# DISC: carbon-carbon distance" << endl;
      cout << "# DISH,DISC" << endl;
      cout << dish << " " << disc << endl;

      //open the filter file...
      ofstream filterfile; filterfile.open("noe.filter");
      filterfile << "TITLE" << endl;
      filterfile << "socially inept, nose-picking (not that i'm not), ultra-sound-whistling rtrds" << endl;
      filterfile << "END" << endl;
      filterfile << "NOEFILTER" << endl;
      filterfile << "# mol residue atom atom mol residue atom atom     r0  filter    noe" << endl;


    //here goes the crappy code
    string resnameA, resnameB;
    int molA=0, molB=0, resnumA=0, resnumB=0, mol=0, atA=0, atB=0, atNumA=0, atNumB=0, count=0;
    int atNOE = 0, totalnoecount = 1;
    double filterbound = 0;
      
    for (int i=0; i < int (noevec.size()); ++i) {
      atNOE = i;
      ostringstream atname;
      ostringstream filteratname;
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
              vatomA = getvirtual(f+addA, NOELIB.NOETYPE, sys);
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
              vatomB = getvirtual(g+addB, NOELIBB.NOETYPE, sys);

	      atname << "# " << count << " " 
		     << (resnumA+1) << NA.resname << " " 
                     << NA.gratomname << " " <<  NA.orgatomname << " " 
                     <<  " "  << (resnumB+1) << NOELIBB.resname << " " 
                     << NOELIBB.gratomname << " " <<  NOELIBB.orgatomname;
 
              filteratname << mol+1 << " "
                           << (resnumA+1) << " " << NA.resname << " " << NA.gratomname << " " <<  NA.orgatomname << " "
			   << molB+1 << " " 
                           << (resnumB+1) << " " << NOELIBB.resname << " " << NOELIBB.gratomname 
                           << " " <<  NOELIBB.orgatomname << " ";

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
              int offsetA = 1, offsetB = 1;
              int creatednoe = 0;

	      for (int va=0; va < (int) vatomA.size(); ++va) {
	       VirtualAtom VA(*vatomA[va]);
               for (int aa=0; aa < 4; ++aa) {
                 int att = VA.operator[](aa);
                 int mol = VA.mol();
                 for(int l=0;l<mol;++l) offsetA += sys.mol(l).numAtoms();
                 atomsA[aa] = att;
		}
               for (int vb=0; vb < (int) vatomB.size(); ++vb) {
		VirtualAtom VB(*vatomB[vb]);
                for (int bb=0; bb < 4; ++bb) {
                 int att = VB.operator[](bb);
                 int mol = VB.mol();                
                 for(int l=0;l<mol;++l) offsetB += sys.mol(l).numAtoms();
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
	      filterbound = bound;
	      //set up corrections
              vector<int> type; 
              type.push_back(VA.type());
              type.push_back(VB.type());
              //check for gromos-correction
              if (gromoscorrection) {
               double cor = 0;
               for (int i=0; i < (int) type.size(); ++i) {
                std::map<int,double>::iterator k = gromoscorrectiondata.find(type[i]);
                if (k != gromoscorrectiondata.end()) cor += k->second;
	       }	       

                if (add) bound += cor;
                else if (sub) bound -= cor;
	      } //end if (gromoscorrection) {

	      //check for multiplicity-correction              
              if (multiplicitycorrection) {
	        double mult = 1;
                
                for (int i=0; i < (int) type.size(); ++i) {
		 std::map<int,double>::iterator k = multiplicitydata.find(type[i]);
                 if (k != multiplicitydata.end()) mult *= k->second; 
		}                

                if (add) bound *= mult;
                else if (sub)  bound /= mult;
	      } //end if (multiplicitycorrection) {
  

	         cout << ss.str() << "   "   << bound << " 1.0000" << endl;
                 if (va > 0 || vb > 0) 	cout << "#automatically generated NOE distance to monitor!" << endl;

		 ++creatednoe;

	       }
	      }
	      cout << atname.str() << endl;

	      //smack out the filterfile...
              
              for (int ii=0; ii < creatednoe; ++ii) {
	       filterfile << totalnoecount << " ";
               filterfile << filteratname.str();
               filterfile << filterbound << " ";
               filterfile << " " << creatednoe << " ";             	       
               for (int iii=0; iii < creatednoe; ++iii) {
		 //                if (ii != iii) filterfile << " " << atNOE+1+iii;
		 if(ii!=iii) filterfile << " " << totalnoecount - ii + iii;
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

vector<VirtualAtom*> getvirtual(int at, int type, System &sys) {


  int mol=0, atNum=0;

  // parse into mol and atom rather than high atom nr.
   while(at >= (atNum+=sys.mol(mol).numAtoms())){
    ++mol;
    if(mol >= sys.numMolecules())
     throw gromos::Exception("prepnoe ", +"Atom number too high in input atom number:\n"+at);
   }
    at-=atNum-sys.mol(mol).numAtoms();

    vector<VirtualAtom*> vat;

    if (type == 4) {
     //we need to automatically generate the r, l atoms...
     vat.push_back(new VirtualAtom(sys,mol,at, type));
     vat.push_back(new VirtualAtom(sys,mol,at, type,1));
    }
    
    else vat.push_back(new VirtualAtom(sys,mol,at, type));     

    return vat;

}
