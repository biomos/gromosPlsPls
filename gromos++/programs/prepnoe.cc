// noeprep.cc; 

#include <args/Arguments.h>
#include <gio/Ginstream.h>
#include <gio/InG96.h>
#include <gcore/System.h>
#include <gio/InTopology.h>
#include <bound/Boundary.h>
#include <args/BoundaryParser.h>
#include "../src/gcore/AtomTopology.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/MoleculeTopology.h"



#include <vector>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <cmath>
#include <sstream>

using namespace gcore;
using namespace args;
using namespace gio;
using namespace bound;
using namespace std;


void Tokenize(const string& str,
                      vector<string>& tokens,
                      const string& delimiters = " ")
{
    // Skip delimiters at beginning.
    string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
    string::size_type pos     = str.find_first_of(delimiters, lastPos);

    while (string::npos != pos || string::npos != lastPos)
    {
        // Found a token, add it to the vector.
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        // Skip delimiters.  Note the "not_of"
        lastPos = str.find_first_not_of(delimiters, pos);
        // Find next "non-delimiter"
        pos = str.find_first_of(delimiters, lastPos);
    }
}


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
  string NOETYPE;

  Noelib(string A, string B, string C, string D){
    resname=A;
    orgatomname=B;
    gratomname=C;
    NOETYPE=D;
  }
  ~Noelib(){}
};

int main(int argc,char *argv[]){


  // Usage string

  string usage = argv[0];
  usage += "\n\t@topo <topology>\n";
  usage += "\t@title <NOE title for output>\n";
  usage += "\t@filter <discard NOE's above a certain distance [nm]; default 0.6 nm>\n";
  usage += "\t@factor <conversion factor Ang <-> nm; default is 10>\n";
  usage += "\t@noe <NOE specification file>\n";
  usage += "\t@lib <NOE specification library>\n";

  // defining all sorts of constants


  // known arguments...
  char *knowns[]={"topo", "title", "filter", "factor", "noe", "lib"};
  int nknowns = 6;
    
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
       double filt=0.6;
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

    if(!nf.check("NOESPEC"))
      throw gromos::Exception("main","NOESPEC file does not contain an NOESPEC block!");

    // in noe all noes will be stored.
    vector<Noeprep> noevec;

    string line;
    while(nf.getline(line),line!="END"){
     vector<string> tokens;
     Tokenize(line, tokens);
     int a=atoi(tokens[0].c_str());
     int b=atoi(tokens[2].c_str());
     double d=atof(tokens[4].c_str());
     //put crap in vector
     //this selects only backbone-backbone NOE's
     //  if (tokens[1] == "HN" && tokens[3] == "HN"){
     noevec.push_back(Noeprep(a,tokens[1],b,tokens[3], d));
     // }
	 }
  nf.close();

    // Read in and create the NOE library
    Ginstream nff(args["lib"]);

       if(!nff.check("NOELIB"))
      throw gromos::Exception("main","NOELIB file does not contain an NOELIB block!");

    // in noe all noes will be stored.
    vector<Noelib> noelib;

 while(nff.getline(line),line!="END"){
     vector<string> tokens;
     Tokenize(line, tokens);
     //put crap in vector
     noelib.push_back(Noelib(tokens[0],tokens[1],tokens[2], tokens[3]));
    }
  nff.close();

     
    nf.close();

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

    //cout the title and that kind of stuff
    cout << "TITLE" << endl;
    cout << "NOE specification file for " << tit << endl;
    cout << "END" << endl;
    cout << "NOE" << endl;



    //here goes the crappy code
    string resnameA, resnameB;
   int molA=0, molB=0, resnumA=0, resnumB=0, mol=0, atA=0, atB=0, atNumA=0, atNumB=0, count=0;
   
  for (int i=0; i < int (noevec.size()); ++i){
   Noeprep NOE = noevec[i];
   //first find the residue-name corresponding to
   //your input-number in the topology

 
   mol=0;
   atNumA=0;
   while(NOE.residA > (atNumA+=sys.mol(mol).topology().numRes())){
        ++mol;
        if(mol > sys.numMolecules())
          throw gromos::Exception("prepnoe" , + "Residue number too high in input line:\n");
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


    int checkA=0;
   //then map the correct gromos-topology atomname
   //to your input-atomname based on the residue-name
   //and the input-atomname
    int p=0;
   for (int k=0;k< int (noelib.size());++k){
    Noelib NOELIB = noelib[k];
     if (NOELIB.resname==resnameA && NOELIB.orgatomname==NOE.atomA){
       //back to topology to get the atom number
	 for(int f=0;f<sys.mol(molA).numAtoms() && checkA==0;++f){
	   if (sys.mol(molA).topology().atom(f).name() == NOELIB.gratomname &&
               sys.mol(molA).topology().resNum(f) == resnumA){
             int addA=0; 
             checkA=1;
             p=k;
	     for(int i=0;i < molA;++i){
               addA+=sys.mol(i).numAtoms();
	     }
             if (NOE.dis/conv > filt){
	       cout << "#";
	     }
             cout << f+addA+1 << NOELIB.NOETYPE << ' ';
	   }
	 }
     }
   }
          if (checkA == 0) {
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
          

	  int checkB=0; 
   for (int z=0;z< int (noelib.size());++z){  
    Noelib NOELIBB = noelib[z];    
     if (NOELIBB.resname==resnameB && NOELIBB.orgatomname==NOE.atomB){
       //back to topology to get the atom number
	 for(int g=0;g<sys.mol(molB).numAtoms() && checkB==0;++g){
	   if (sys.mol(molB).topology().atom(g).name() == NOELIBB.gratomname &&
               sys.mol(molB).topology().resNum(g) == resnumB){
             int addB=0;
             checkB=1;
             count+=1;
	     for(int i=0;i < molB;++i){
               addB+=sys.mol(i).numAtoms();
	     }
             cout << g+addB+1 << NOELIBB.NOETYPE << ' ' <<  NOE.dis/conv 
                  <<   " #"   << count << ' ' << resnumA+1 << NA.resname
                  << ' '      << NA.gratomname << ' ' <<  NA.orgatomname
                  << "  "     << resnumB+1 << NOELIBB.resname
                  << ' '      << NOELIBB.gratomname << ' ' <<  NOELIBB.orgatomname << endl;
	   }
	 }
     }
   }
  
  
 if (checkB == 0) {
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
  
  }

  
  cout << "END" << endl;            
         
  }
  
  catch(gromos::Exception e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}
