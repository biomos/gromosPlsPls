// mkscript.cc

#include <cassert>

#include "../src/args/Arguments.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/Solvent.h"
#include "../src/gcore/SolventTopology.h"
#include "../src/gio/InTopology.h"
#include "../src/gio/Ginstream.h"
#include <fstream>
#include <sstream>
#include <string>
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <unistd.h>

using namespace std;
using namespace gcore;
using namespace gio;
using namespace args;

#include "mkscript.h"

void printWarning(int &numw, int &nume, string s);
void printError(int &numw, int &nume, string s);
void printInput(string ifile, string ofile, int nstlim, double t, double dt);
void readLibrary(string file,  vector<filename> &names,
		 vector<filename> &misc, 
		 vector<string> &linknames, vector<int> &linkadditions, 
		 string system, string q, string submitcommand, double t, 
		 double dt, int &w, int &e, int ns);


int main(int argc, char **argv){

  char *knowns[] = {"sys", "script", "bin", "dir", "queue", 
		    "files", "template", "XX", "cmd"};
  int nknowns = 9;

  string usage = argv[0];
  usage += "\n\t@sys  <system name>\n";
  usage += "\t@script <script number> <number of scripts>\n";
  usage += "\t@bin    <gromos96 binary to use>\n";
  usage += "\t[@dir    <where should the files be>]\n";
  usage += "\t[@queue  <which queue?>]\n";
  usage += "\t@files\n";
  usage += "\t\ttopo        <topology>\n";
  usage += "\t\tinput       <input file>\n";
  usage += "\t\tcoord       <initial coordinates>\n";
  usage += "\t\t[refpos     <reference positions>]\n";
  usage += "\t\t[posresspec <position restraints specifications>]\n";
  usage += "\t\t[disres     <distance restraints>]\n";
  usage += "\t\t[dihres     <dihedral restraints>]\n";
  usage += "\t\t[jvalue     <j-value restraints>]\n";
  usage += "\t\t[ledih      <local elevation dihedrals>]\n";
  usage += "\t\t[pttopo     <perturbation topology>]\n";
  usage += "\t[@template   <template filenames>]\n";
  usage += "\t[@XX         gromosXX script]\n";
  usage += "\t[@cmd        <last command>]\n";

  try{
    Arguments args(argc, argv, nknowns, knowns, usage);

    // set the number of warnings and the number of errors
    int numWarnings=0;
    int numErrors=0;

    // first get some input parameters
    int scriptNumber=1, numScripts=1;
    string simuldir,q, submitcommand;
    {
      
      Arguments::const_iterator iter=args.lower_bound("script");
      if(iter!=args.upper_bound("script")){
	  scriptNumber=atoi(iter->second.c_str());
	  ++iter;
      }
      if(iter!=args.upper_bound("script"))
	  numScripts=atoi(iter->second.c_str());
      iter=args.lower_bound("dir");
      if(iter!=args.upper_bound("dir"))
	 simuldir=iter->second;
      else
	  simuldir="`pwd`";
      q="\"put your favourite queue name here\"";
      if(args.count("queue")>0) q=args["queue"];
      if(q=="penguin")
	submitcommand="psub -s "+q+" 2 ";
      else
	// if(q=="oxen" || q=="moose" || q=="ccpc")
	submitcommand="ssub -s "+q+" ";
    }
    string systemname=args["sys"];
    
    // parse the files
    int l_coord = 0,l_topo=0, l_input=0, l_refpos=0, l_posresspec=0;
    int l_disres = 0, l_dihres=0, l_jvalue=0, l_ledih=0, l_pttopo=0;
    string s_coord, s_topo, s_input, s_refpos, s_posresspec;
    string s_disres, s_dihres, s_jvalue, s_ledih, s_pttopo;
    for(Arguments::const_iterator iter=args.lower_bound("files"),
	    to=args.upper_bound("files"); iter!=to; ++iter){
      switch(FILETYPE[iter->second]){
	case coordfile  : ++iter; s_coord=iter->second;  l_coord=1;  break;
	case inputfile  : ++iter; s_input=iter->second;  l_input=1;  break;
	case topofile   : ++iter; s_topo=iter->second;   l_topo=1;   break;
	case refposfile : ++iter; s_refpos=iter->second; l_refpos=1; break;
	case posresspecfile :++iter; s_posresspec=iter->second; 
	                                             l_posresspec=1; break;
	case disresfile : ++iter; s_disres=iter->second; l_disres=1; break;
	case dihresfile : ++iter; s_dihres=iter->second; l_dihres=1; break;
	case jvaluefile : ++iter; s_jvalue=iter->second; l_jvalue=1; break;
	case ledihfile  : ++iter; s_ledih=iter->second;  l_ledih=1;  break;
        case pttopofile : ++iter; s_pttopo=iter->second; l_pttopo=1; break;
        case outputfile : ++iter; printWarning(numWarnings, numErrors, 
					  iter->second+" not used"); break;
	case outtrxfile : ++iter; printWarning(numWarnings, numErrors, 
					  iter->second+" not used"); break;
	case outtrvfile : ++iter; printWarning(numWarnings, numErrors, 
					  iter->second+" not used"); break; 
	case outtrefile : ++iter; printWarning(numWarnings, numErrors, 
					  iter->second+" not used"); break;
	case outtrgfile : ++iter; printWarning(numWarnings, numErrors, 
					  iter->second+" not used"); break;
	case scriptfile : ++iter; printWarning(numWarnings, numErrors, 
					  iter->second+" not used"); break;
	case unknownfile: printError(numWarnings, numErrors,
	      "Don't know how to handle file "+iter->second);
      }
    }

    // check which outformat we want (gromos96 or gromosXX)
    bool gromosXX = false;
    if (args.count("XX") != -1)
      gromosXX = true;
    
    // read topology
    if(!l_topo){
      throw gromos::Exception("mkscript", "You have to specify a topology\n"+usage);
    }
    InTopology it(s_topo);
    System sys=it.system();
    
    // read the input file
    if(!l_input){
      throw gromos::Exception("mkscript", "You have to specify an input file\n"+usage);
    }
    Ginstream imd(s_input);
    
    input gin;
    imd >> gin;
    
    imd.close();
 
    // create names for automated file names
    vector<filename> filenames;
    vector<filename> misc;
    vector<int>      linkadditions;
    vector<string>   linknames;
    
    for(int i=0; i<numFiletypes; i++){
      filename newname(systemname, gin.step.t, gin.step.nstlim*gin.step.dt, 
		       scriptNumber, q);
      filenames.push_back(newname);
    }
    for(int i=0; i<2; i++){
      filename newname(systemname, gin.step.t, gin.step.nstlim*gin.step.dt, 
		       scriptNumber, q);
      misc.push_back(newname);
    }
    
    // set the standard templates
    filenames[FILETYPE["script"]].setTemplate("jmd%system%_%number%.sh");
    filenames[FILETYPE["input"]].setTemplate("imd%system%_%number%.dat");
    filenames[FILETYPE["topo"]].setTemplate("%system%mta1.dat");
    filenames[FILETYPE["refpos"]].setTemplate("%system%px_%number%.dat");
    filenames[FILETYPE["posresspec"]].setTemplate("%system%pr_%number%.dat");
    filenames[FILETYPE["disres"]].setTemplate("%system%drts_%number%.dat");
    filenames[FILETYPE["pttopo"]].setTemplate("%system%pt1.dat");
    filenames[FILETYPE["dihres"]].setTemplate("%system%arts_%number%.dat");
    filenames[FILETYPE["jvalue"]].setTemplate("%system%jlts_%number%.dat");
    filenames[FILETYPE["ledih"]].setTemplate("%system%lets_%number%.dat");
    filenames[FILETYPE["coord"]].setTemplate("o%system%sxmd_%number%.dat");
    filenames[FILETYPE["output"]].setTemplate("omd%system%_%number%.out");
    filenames[FILETYPE["outtrx"]].setTemplate("o%system%trmd_%number%.dat");
    filenames[FILETYPE["outtrv"]].setTemplate("o%system%tvmd_%number%.dat");
    filenames[FILETYPE["outtre"]].setTemplate("o%system%temd_%number%.dat");
    filenames[FILETYPE["outtrg"]].setTemplate("o%system%tgmd_%number%.dat");

    misc[0].setTemplate("/scrloc/${NAME}_%system%_%number%");
    misc[1].setTemplate(submitcommand+filenames[FILETYPE["script"]].temp());
    
    // read in the library
    if (args.count("template")>=0){
      int really_do_it=1;
      string libraryfile;
      if(args.count("template")==0){
	if(getenv("MKSCRIPT_TEMPLATE")){
	  libraryfile=getenv("MKSCRIPT_TEMPLATE");
	}
	else{
	  ostringstream os;
	  os << "Trying to read template file, but MKSCRIPT_TEMPLATE is not set\n"
	     << "Either specify a filename or set this environment variable\n"
	     << "Using defaults now\n";
	  printWarning(numWarnings, numErrors, os.str());
	  really_do_it=0;
	}
		
      }
      
      if(args.count("template")>0)
	libraryfile=args["template"];
      if(really_do_it)
	// And here is a gromos-like function call!
	readLibrary(libraryfile, filenames, misc, 
		    linknames, linkadditions, 
		    systemname, q, submitcommand, gin.step.t, 
		    gin.step.nstlim*gin.step.dt, numWarnings, numErrors,
		    scriptNumber);     
    }

    // overwrite last command if given as argument
    if (args.count("cmd") > 0){
      ostringstream os;
      Arguments::const_iterator iter = args.lower_bound("cmd"),
	to = args.upper_bound("cmd");
      for(; iter != to; ++iter){
	std::string s = iter->second;
	if(s.find("\\n")!=string::npos) 
	  s.replace(s.find("\\n"), 2, "\n");
	else s+=" ";
	
	// os << iter->second << " ";
	os << s;
	
      }
      
      misc[1].setTemplate(os.str());
    }
    
    // read what is in the coordinate file
    fileInfo crd;
    if(!l_coord){
      // try to open it from the template
      ifstream fin(filenames[FILETYPE["coord"]].name(-1).c_str());
      if(fin){
	ostringstream os;
	os << "No coordinate file specified, but I found "
	   << filenames[FILETYPE["coord"]].name(-1)
	   << " which I will use\n";
	printWarning(numWarnings, numErrors, os.str());
	s_coord=filenames[FILETYPE["coord"]].name(-1);
	l_coord=1;
      }
      else{
	ostringstream os;
	os << "No coordinate file is specified, some checks are not performed\n";
	os << "Assuming it does not exist yet, I will use " 
	   << filenames[FILETYPE["coord"]].name(-1) 
	   << " in the script\n";
	printWarning(numWarnings, numErrors, os.str());
      }
     
    }
    if(l_coord){
      Ginstream icrd(s_coord);
      icrd >> crd;
      icrd.close();
    }

    // calculate some standard numbers
    int numSoluteAtoms=0;
    for(int i=0; i< sys.numMolecules(); i++) 
	numSoluteAtoms+=sys.mol(i).topology().numAtoms();
    int numSolventAtoms=sys.sol(0).topology().numAtoms();
    int numTotalAtoms=gin.system.npm * numSoluteAtoms +
	              gin.system.nsm * numSolventAtoms;
    
    // carry out a thousand tests:

    // Does the binary exist?
    {
      ifstream fin(args["bin"].c_str());
      if(!fin)
	printWarning(numWarnings, numErrors, "Specified binary not found! "
		     +args["bin"]+"\n");
      else fin.close();
    }
    
    //SYSTEM block
    if(gin.system.found){
      if(l_coord){
        for(unsigned int i=0; i<crd.blocks.size(); i++){
          if(crd.blocks[i]=="POSITION"||
	     crd.blocks[i]=="VELOCITY"||
	     crd.blocks[i]=="REFPOSITION"||
	     crd.blocks[i]=="REDPOSITION"){
	    if(numTotalAtoms!=crd.blockslength[i]){
	      ostringstream os;
	      os << "From topology and SYSTEM block, I calculate "
	         << numTotalAtoms << " atoms in the system" << endl;
	      os << "But coordinate file has " << crd.blockslength[i]
	         << " lines in " << crd.blocks[i] << " block." << endl;
	      os << "Maybe NSM should be " 
	         << (crd.blockslength[i] - numSoluteAtoms)/numSolventAtoms 
	         << " ?" << endl;
	      printError(numWarnings, numErrors,os.str());
  	    }
          }
        }
      }
    }
    else printError(numWarnings, numErrors, 
		    "Could not find SYSTEM block\n");

 
    // START
    if(gin.start.found){
      if(l_coord){
        int velblock=0;
        for(unsigned int i=0; i< crd.blocks.size(); i++)
	  if(crd.blocks[i]=="VELOCITY") velblock=1;
        if(velblock&&gin.start.ntx==1){
	  ostringstream os;
	  os << "NTX = 1 in START block, which means that you don't want to "
             << "read in the\n";
	  os << "velocity block I found in the coordinate file.\n";
	  printWarning(numWarnings, numErrors, os.str());
        }
        if(gin.start.ntx>1&&!velblock){
	  ostringstream os;
          os << "NTX = " << gin.start.ntx << " in START block means that you "
	     << "want to read in the velocities from the coordinate file.\n";
	  os << "But coordinate file does not have a VELOCITY block.\n";
          printError(numWarnings, numErrors, os.str());
        }
      }
    }
    else printError(numWarnings, numErrors, 
		    "Could not find START block\n");

    // STEP
    if(!gin.step.found)
	printError(numWarnings, numErrors,
		   "Could not find STEP block\n");

    // BOUNDARY
    double box[3];
    if(gin.boundary.found){
      int boxblock=0;
      if(gin.boundary.nrdbox && l_coord){
        for(unsigned int i=0; i< crd.blocks.size(); i++)
	  if(crd.blocks[i]=="BOX" || 
	     crd.blocks[i]=="TRICLINICBOX") boxblock=1;
        if(!boxblock){
	  string s="NRDBOX = 1 in BOUNDARY block, ";
          s+="but no BOX block in coordinate file.\n"; 
	  printError(numWarnings, numErrors, s); 
        }
        else 
	  for(int i=0;i<3;i++)
	    box[i]=crd.box[0];
      }
      if(gin.boundary.nrdbox==0)
	for(int i=0;i<3;i++)
	  box[i]=gin.boundary.box[i];
      if(gin.boundary.ntb < 0 && (gin.boundary.nrdbox==0 || boxblock)){
	if(box[0]!=box[1] || box[0]!=box[2] || box[1]!=box[2]){
	  ostringstream os;
	  os << "NTB = " << gin.boundary.ntb << " in BOUNDARY block "
	     << "means truncated octahedron.\n";
	  os << "But boxdimensions (from ";
	  if(boxblock) os << "coordinate file ";
	  else os << "input file ";
	  os << ") are not the same:\n";
	  os << box[0] << "\t" << box[1] << "\t" << box[2] << endl;
          printError(numWarnings, numErrors, os.str());
	}
      }
    }
    else printError(numWarnings,  numErrors, 
		    "Could not find BOUNDARY block\n");

    // SUBMOLECULES
    if(gin.submolecules.found){
      int na=0;
      int error=0;
      if(int(gin.submolecules.nsp.size())!=sys.numMolecules()) error=1;
      else
	for(int i=0;i<sys.numMolecules();i++){
	  na+=sys.mol(i).topology().numAtoms();
	  if(na!=gin.submolecules.nsp[i]) error=1;
	}
      if(error){
	ostringstream os;
        na=0;
	os << "SUBMOLECULES block does not match topology; should be\n";
        os << "SUBMOLECULES\n";
        os << "#     NSPM  NSP(1.. NSPM)\n";
	os << setw(10) << sys.numMolecules();
        int countsbm=0;
	
	for(int i=0; i<sys.numMolecules(); i++){
	  na+=sys.mol(i).topology().numAtoms();
	  os << setw(6) << na;
          countsbm++;
	  if(countsbm%8==0) os << endl;
	  
	}
        os << "\nEND\n";
        printError(numWarnings, numErrors, os.str());
      }
    }
    else printError(numWarnings, numErrors, 
		    "Could not find SUBMOLECULES block\n");

    // TCOUPLE
    if(gin.tcouple.found){
      if(gin.system.nsm==0 && gin.tcouple.ntt[2]!=0){
        string s="There are no solvent ";
        s+="molecules specified in the SYSTEM block, but\n";
        s+="you do couple its temperature in the TCOUPLE ";
        s+="block (third line)\n";
        printError(numWarnings, numErrors, s);
      }
    }

    // PCOUPLE
    if(gin.pcouple.found){
      if(gin.pcouple.ntp!=0 && abs(gin.boundary.ntb)!=2){
        ostringstream os;
        int ntb;
        if(gin.boundary.ntb<0) ntb = -2;
        else ntb=2;
        os << "NTP = " << gin.pcouple.ntp << " in PCOUPLE block, but "
	   << "NTB = " << gin.boundary.ntb << " in BOUNDARY block, which "
	   << "means that we do not calculate the virial\n";
        os << "Set NTB = " << ntb << " to calculate the virial\n";
        printError(numWarnings, numErrors, os.str());
      }
      if(gin.pcouple.ntp==2 && gin.boundary.ntb <0){
        ostringstream os;
        os << "NTP = " << gin.pcouple.ntp << " in PCOUPLE block specifies "  
  	   << "anisotropic pressure scaling.\n";
        os << "But NTB = " << gin.boundary.ntb << " in BOUNDARY block means "
	   << "truncated octahedral boundary conditions\n";
        printError(numWarnings, numErrors, os.str());
      }
      if(gin.tcouple.found){
        double tautmax=0;
        for(int i=0;i<3;i++)
	  if(gin.tcouple.ntt[i]!=0 && gin.tcouple.taut[i]>tautmax)
	    tautmax=gin.tcouple.taut[i];
        if(tautmax>=gin.pcouple.taup){
	  ostringstream os;
	  os << "TAUP in PCOUPLE block (" << gin.pcouple.taup << ") is not "
	     << "larger than TAUT in TCOUPLE block (" << tautmax << ").\n";
	  printWarning(numWarnings, numErrors, os.str());
        }
      }
    }

    //CENTROFMASS
    //PRINT
    if(!gin.print.found) printError(numWarnings, numErrors,
				    "Could not find PRINT block\n");

    //WRITE

    //SHAKE
    if(((!gin.shake.found || gin.shake.ntc==1) && gin.step.dt > 0.0005) ||
       (( gin.shake.found && gin.shake.ntc==2) && gin.step.dt > 0.001) ||
       (( gin.shake.found && gin.shake.ntc==3) && gin.step.dt > 0.002)) {
      ostringstream os;
      string comment;
      double suggest=0.0005;
      if(!gin.shake.found){
	  comment = "no SHAKE block: no shake on solute"; suggest=0.0005; 
          gin.shake.ntc=1;}
      else if(gin.shake.ntc==1){
	  comment = "no shake on solute"; suggest=0.0005;}
      else if(gin.shake.ntc==2){
	  comment = "shake bonds with H"; suggest=0.001;}
      else if(gin.shake.ntc==3){
	  comment = "shake all bonds";    suggest=0.002;}
      os << "DT in STEP block is set to " << gin.step.dt << ", which is "
	 << "considered to be too large\n";
      os << "if NTC = " << gin.shake.ntc << " in SHAKE block.\n";
      os << "For NTC = " << gin.shake.ntc << " (" << comment << ") rather "
         << "use DT = " << suggest << ".\n";
      printWarning(numWarnings, numErrors, os.str());
    }

    // FORCE
    if(gin.force.found){
      string comment;
      if(gin.shake.ntc==2) comment = "shake bonds with H"; 
      if(gin.shake.ntc==3) comment = "shake all bonds";

      if((gin.shake.ntc>1  && gin.force.ntf[0]) ||
         (gin.shake.ntc==3 && gin.force.ntf[1])){
        ostringstream os;
        os << "NTC = " << gin.shake.ntc << " in SHAKE block ("
	   << comment << ")\n";
        os << "so there is no need to calculate the forces due to these "
	   << "bonds, as indicated in the FORCE block.\n";
        printWarning(numWarnings, numErrors, os.str());
      }
      if((gin.shake.ntc==1&&(gin.force.ntf[0]==0||gin.force.ntf[1]==0))||
         (gin.shake.ntc==2&&gin.force.ntf[1]==0)){
	ostringstream os;
        os << "NTC = " << gin.shake.ntc << " in SHAKE block (";
	os << comment << ")\n";
        os << "But you do not want to calculate the forces on non-shaken "
	   << "bonds?\n";
	os << "ntf[0] = " << gin.force.ntf[0] << " ntf[1] = "
	   << gin.force.ntf[1] << endl;
        printWarning(numWarnings, numErrors, os.str());
      }
      if(gin.force.nre[gin.force.nre.size()-1]!=numTotalAtoms){
        ostringstream os;
        os << "The last energy group in the FORCE block ("
  	   << gin.force.nre[gin.force.nre.size()-1] << ") should be equal "
	   << "to the total number of atoms in the system: "
	   << numTotalAtoms << ".\n";
        printError(numWarnings, numErrors, os.str());
      }
    }
    else printError(numWarnings, numErrors, 
		    "Could not find FORCE block\n");
    //PLIST
    if(gin.plist.found){
      if(gin.plist.rcutp>gin.plist.rcutl){
        ostringstream os;
        os << "In PLIST block RCUTP = " <<gin.plist.rcutp <<  " and RCUTL = "
	   << gin.plist.rcutl << endl;
        os << "RCUTP should be less than or equal to RCUTL\n";
        printError(numWarnings, numErrors, os.str());
      }
      double minbox=1e6;
      int trunc=0;
      if(gin.boundary.ntb!=0){
	
	if(gin.boundary.nrdbox==0 || (gin.boundary.nrdbox==1 && l_coord))
	  for(int i=0;i<3; i++) if(box[i]<minbox) minbox=box[i];
	if(gin.boundary.ntb<0) {
	  trunc=1;
	  // check formula
	  minbox*=0.5*1.732051;
	}
	if(minbox<2*gin.plist.rcutl){
	  ostringstream os;
	  os << "RCUTL in PLIST block is " << gin.plist.rcutl << endl;
	  os << "for a ";
	  if(trunc) os << "truncated octahedral ";
	  else os << "rectangular ";
          os << "box with these dimensions\n";
          for(int i=0; i<3; i++) os << "\t" << box[i];
          if(gin.boundary.nrdbox==1) os << "\t(from coordinate file)\n";
	  else os << "\t(from input file)\n";
          os << "this is too long\n";
          printWarning(numWarnings, numErrors, os.str());
	}
      }
    }
    else printError(numWarnings, numErrors, 
		    "Could not find PLIST block\n");

    //LONGRANGE
    if(gin.longrange.found){
      if((gin.longrange.epsrf!=1.0) && (gin.plist.rcutl!=fabs(gin.longrange.rcrf))){
        ostringstream os;
        os << "We usually expect RCRF in the LONGRANGE block to be equal to "
	   << "RCUTL in the PLIST block\n";
        os << "You specified RCUTL = " << gin.plist.rcutl << " and RCRF = "
	   << gin.longrange.rcrf << ".\n";
        printWarning(numWarnings, numErrors, os.str());
      }
    }
    else printError(numWarnings, numErrors, 
		    "Could not find LONGRANGE block\n");
    
    //POSREST
    if(gin.posrest.found&&gin.posrest.ntr!=0){
      int refblock=0;
      int numref=0;
      if(gin.posrest.nrdrx==1&&l_coord){
        for(unsigned int i=0 ;i<crd.blocks.size(); i++){
	  if(crd.blocks[i]=="REFPOSITION") {
	    refblock=1; numref=crd.blockslength[i];
	  }
	}
        if(!refblock){
	  string s="NRDRX=1 in POSREST block, but no REFPOSITION block ";
	  s+="in coordinate file\n";
	  printError(numWarnings, numErrors, s);
	}
      }
      else{
	//first find out if the file exists
	if(!l_refpos && !gin.posrest.nrdrx){
          if(!l_refpos){
	    ifstream fin(filenames[refposfile].name(0).c_str());
	    if(fin) {
	      ostringstream os;
	      os << "No refpos-file specified, but I found "
		 << filenames[refposfile].name(0)
		 << " which I will use\n";
	      l_refpos=1;
	      s_refpos=filenames[refposfile].name(0);
	      printWarning(numWarnings, numErrors, os.str());
	      fin.close();
	    }
	    else{
	      ostringstream os;
	      os << "NRDRX = " << gin.posrest.nrdrx << " in POSREST block, but no "
		 << "refpos-file specified\n";
	      printError(numWarnings, numErrors, os.str());
	    }
	  }
	  
	}
	else if(l_refpos){
	  Ginstream irfp(s_refpos);
          fileInfo rfp;
          irfp >> rfp;
          irfp.close();
          for(unsigned int i=0; i<rfp.blocks.size(); i++){
	    if(rfp.blocks[i]=="REFPOSITION") {
	      refblock=1; numref=rfp.blockslength[i];
	    }
	  }
          if(!refblock){
	    ostringstream os;
	    os << "No REFPOSITION block in refpos file (" << s_refpos 
	       << ")\n";
            printError(numWarnings, numErrors, os.str());
	  }
	}
      }
      if(refblock&&numref!=numTotalAtoms){
	ostringstream os;
	os << "Number of atoms in REFPOSITION block in ";
        if(gin.posrest.nrdrx) os << s_coord;
	else os << s_refpos;
	os << " (" << numref << ")\n does not match total number of atoms ("
	   << numTotalAtoms << ")\n";
        printError(numWarnings, numErrors, os.str());
      }
      if(!l_posresspec){
	ifstream fin(filenames[posresspecfile].name(0).c_str());
	if(fin){
	  l_posresspec=1;
	  s_posresspec=filenames[posresspecfile].name(0);
	  ostringstream os;
	  os << "No posresspec-file specified, but I found"
	     << filenames[posresspecfile].name(0)
	     << " which I will use\n";
	  printWarning(numWarnings, numErrors, os.str());
	  fin.close();
	}
	else{
	  ostringstream os;
	  os << "NTR = " << gin.posrest.ntr << " in POSREST block, but no "
	     << "posresspec-file specified\n";
	  printError(numWarnings, numErrors, os.str());
	}
      }
      if(l_posresspec){
	Ginstream iprs(s_posresspec);
	fileInfo prs;
	iprs >> prs;
	iprs.close();
        int l_prs=0;
	for(unsigned int i=0; i< prs.blocks.size(); i++)
	  if(prs.blocks[i]=="POSRESSPEC") l_prs=1;
        if(!l_prs){
	  string s="No POSRESSPEC block in posresspec file (";
	  s+=s_posresspec;
	  s+=")\n";
          printError(numWarnings, numErrors, s);
	}
      }
    }
    else if(l_refpos||l_posresspec){
      ostringstream os;
      if(l_refpos) os << "reference positions file ";
      if(l_refpos && l_posresspec) os << "and ";
      if(l_posresspec) os << "position restraints specification file ";
      os << "specified\n";
      os << "But no position res/constraining according to input file\n";
      printWarning(numWarnings, numErrors, os.str());
    }

    //PERTURB
    if(gin.perturb.found&&gin.perturb.ntg!=0){
      // check if perturbation topology is present
      if(!l_pttopo){
	ifstream fin(filenames[pttopofile].name(0).c_str());
	if(fin){
	  l_pttopo=1;
	  s_pttopo=filenames[pttopofile].name(0);
	  ostringstream os;
	  os << "No perturbation topology specified, but I found "
	     << filenames[pttopofile].name(0)
	     << " which I will use\n";
	  printWarning(numWarnings, numErrors, os.str());
	  fin.close();
	}
	else{
	  ostringstream os;
	  os << "NTG = " << gin.perturb.ntg << " in PERTURB block, but "
	     << "no perturbation topology specified\n";
	  printError(numWarnings, numErrors, os.str());
	}
      }

      double rlamfin=gin.perturb.rlam + gin.perturb.dlamt * gin.step.dt * 
	gin.step.nstlim;
      if(rlamfin>1.0){
	ostringstream os;
	os << "Using RLAM = " << gin.perturb.rlam << " and DLAMT = "
	   << gin.perturb.dlamt << " in the PERTURB block and NSTLIM = "
	   << gin.step.nstlim << " in the STEP block\n";
	os << "will lead to a final lambda value of " <<rlamfin << endl;
	printWarning(numWarnings, numErrors, os.str());
      }
    }
    else if(l_pttopo){
      string s="Perturbation topology specified, but no perturbation ";
      s+="according to the input file\n";
      printWarning(numWarnings, numErrors, s);
    }
      
    // Now, write the stupid script
    if(numErrors!=0){
      if(numErrors>1) cout << "\n\nTHERE WERE " << numErrors << " ERRORS\n";
      else cout << "\n\nTHERE WAS 1 ERROR\n";
      cout << "No script will be written\n";
      exit(1);
    }

    vector<string> scripts(numScripts+2);
    for(int i =0; i<numScripts; i++){
      ofstream fout(filenames[FILETYPE["script"]].name(i).c_str());
      scripts[i+2]=filenames[FILETYPE["script"]].name(i);
      fout.setf(ios_base::left, ios_base::adjustfield);
      fout << "#!/bin/sh" << endl;
      fout << "\n# first we set some variables\n";
      fout << "NAME=`whoami`\n";
      fout << "PROGRAM=" << args["bin"] << endl;
      fout << "SIMULDIR=" <<  simuldir << endl;
      fout << "\n# create temporary directory\n";
      fout << "WORKDIR=" << misc[0].name(i) << endl;
      fout << "mkdir -p ${WORKDIR}\n";
      fout << "cd       ${WORKDIR}\n";
      fout << "\n# set the input files\n";
      fout << "TOPO=${SIMULDIR}/" << s_topo << endl;
      fout << "IUNIT=${SIMULDIR}/";
      if(i==0) fout << s_input << endl;
      else  {
	fout << filenames[FILETYPE["input"]].name(i) << endl;
        printInput(s_input, filenames[FILETYPE["input"]].name(i), 
		   gin.step.nstlim, gin.step.t+i*gin.step.nstlim*gin.step.dt, 
		   gin.step.dt);
      }
      fout << "INPUTCRD=${SIMULDIR}/";
      if(i==0&&l_coord) fout << s_coord << endl;
      else              fout << filenames[FILETYPE["coord"]].name(i-1) << endl;
      if(l_refpos) fout << "REFPOS=${SIMULDIR}/" << s_refpos << endl;
      if(l_posresspec) fout << "POSRESSPEC=${SIMULDIR}/" 
			    << s_posresspec << endl;
      if(l_disres) fout << "DISRES=${SIMULDIR}/" << s_disres << endl;
      if(l_dihres) fout << "DIHRES=${SIMULDIR}/" << s_dihres << endl;
      if(l_jvalue) fout << "JVALUE=${SIMULDIR}/" << s_jvalue << endl;
      if(l_ledih)  fout << "LEDIH=${SIMULDIR}/"  << s_ledih  << endl;
      if(l_pttopo) fout << "PTTOPO=${SIMULDIR}/" << s_pttopo << endl;
      // any additional links?
      for(unsigned int k=0; k<linkadditions.size(); k++)
	if(linkadditions[k]<0)
	  fout << linknames[k] <<  "=${SIMULDIR}/" 
	       << filenames[numFiletypes+k].name(i)
	       << endl;
	  
      fout << "\n#set the output files\n";
      fout << "OUNIT=" << filenames[FILETYPE["output"]].name(i) << endl;
      fout << "OUTPUTCRD=" << filenames[FILETYPE["coord"]].name(i) << endl;
      if(gin.write.ntwx) fout << "OUTPUTTRX=" 
			      << filenames[FILETYPE["outtrx"]].name(i) 
			      << endl;
      if(gin.write.ntwv) fout << "OUTPUTTRV=" 
			      << filenames[FILETYPE["outtrv"]].name(i)
			      << endl;
      if(gin.write.ntwe) fout << "OUTPUTTRE=" 
			      << filenames[FILETYPE["outtre"]].name(i) 
			      << endl;
      if(gin.write.ntwg) fout << "OUTPUTTRG=" 
			      << filenames[FILETYPE["outtrg"]].name(i) 
			      << endl;
      // any additional links?
      for(unsigned int k=0; k<linkadditions.size(); k++)
	if(linkadditions[k]>0)
	  fout << linknames[k] << "="
	       << filenames[numFiletypes+k].name(i) << endl;

      if (!gromosXX){
	fout << "\n# link the files\n";
	fout << "rm -f fort.*\n";
	fout << setw(25) << "ln -s ${TOPO}" << " fort.20\n";
	fout << setw(25) << "ln -s ${INPUTCRD}" << " fort.21\n";
	if(l_refpos)     fout << setw(25) 
			      << "ln -s ${REFPOS}" << " fort.22\n";
	if(l_posresspec) fout << setw(25) 
			      << "ln -s ${POSRESSPEC}" << " fort.23\n";
	if(l_disres)     fout << setw(25) 
			      << "ln -s ${DISRES}" << " fort.24\n";
	if(l_dihres)     fout << setw(25) 
			      << "ln -s ${DIHRES}" << " fort.25\n";
	if(l_jvalue)     fout << setw(25) 
			      << "ln -s ${JVALUE}" << " fort.26\n";
	if(l_ledih)      fout << setw(25)
			      << "ln -s ${LEDIH}" << " fort.27\n";
	if(l_pttopo)     fout << setw(25)
			      << "ln -s ${PTTOPO}" << " fort.30\n";

	fout << setw(25) << "ln -s ${OUTPUTCRD}" << " fort.11\n";
	if(gin.write.ntwx) fout << setw(25)
				<< "ln -s ${OUTPUTTRX}" << " fort.12\n";
	if(gin.write.ntwv) fout << setw(25) 
				<< "ln -s ${OUTPUTTRV}" << " fort.13\n";
	if(gin.write.ntwe) fout << setw(25)
				<< "ln -s ${OUTPUTTRE}" << " fort.15\n";
	if(gin.write.ntwg) fout << setw(25)
				<< "ln -s ${OUTPUTTRG}" << " fort.16\n";
	// any additional links
	for(unsigned int k=0; k<linkadditions.size(); k++){
	  string s("ln -s ${" + linknames[k] +"}");
	  fout << setw(25) << s << " fort." << abs(linkadditions[k]) 
	       << endl;
	}
	
	fout << "\n# run the program\n\n";
	fout << "${PROGRAM} < ${IUNIT} > ${OUNIT}\n";
      }
      else{

	if (q=="penguin")
		fout << "export LD_LIBRARY_PATH=/orc/markus/machines/isengard/programs/lib\n\n";	

	fout << "\n\n${PROGRAM}";
	fout << " \\\n\t" << setw(12) << "@topo" << " ${TOPO}";
	fout << " \\\n\t" << setw(12) << "@conf" << " ${INPUTCRD}";
	fout << " \\\n\t" << setw(12) << "@input" << " ${IUNIT}";
	if (l_pttopo)       fout << " \\\n\t" 
				 << setw(12) << "@pert" << " ${PTTOPO}";
	if (l_jvalue)       fout << " \\\n\t"
				 << setw(12) << "@jval" << " ${JVALUE}";
	fout << " \\\n\t" << setw(12) << "@fin" << " ${OUTPUTCRD}";
	if (gin.write.ntwx) fout << " \\\n\t" << setw(12) << "@trj"
				 <<" ${OUTPUTTRX}";
	if (gin.write.ntwv) fout << " \\\n\t" << setw(12) << "@trv"
				 << " ${OUTPUTTRV}";
	if (gin.write.ntwe) fout << " \\\n\t" << setw(12) << "@tre"
				 << " ${OUTPUTTRE}";
	if (gin.write.ntwg) fout << " \\\n\t" << setw(12) << "@trg"
				 << " ${OUTPUTTRG}";
	// any additional links
        for(unsigned int k=0; k<linkadditions.size(); k++)
	  fout << " \\\n\t@" << setw(11) <<  linknames[k] 
	       << " ${" << linknames[k]<< "}";
	
	fout << "\\\n\t" << setw(12) << ">" << " ${OUNIT}\n\n";

      }
      
      fout << "uname -a >> ${OUNIT}\n";
      
      if(gin.write.ntwx||gin.write.ntwv||gin.write.ntwe||gin.write.ntwg)
	fout << "\n# compress some files\n";
      if(gin.write.ntwx) fout << "gzip ${OUTPUTTRX}\n";
      if(gin.write.ntwv) fout << "gzip ${OUTPUTTRV}\n";
      if(gin.write.ntwe) fout << "gzip ${OUTPUTTRE}\n";
      if(gin.write.ntwg) fout << "gzip ${OUTPUTTRG}\n";

      fout << "\n# copy the files back\n";
      fout << "OK=1\n";
      fout << setw(25) << "cp ${OUNIT}" << " ${SIMULDIR} || OK=0\n";
      fout << setw(25) << "cp ${OUTPUTCRD}" << " ${SIMULDIR} || OK=0\n";
      if(gin.write.ntwx) fout << setw(25) << "cp ${OUTPUTTRX}.gz" 
			      << " ${SIMULDIR} || OK=0\n";
      if(gin.write.ntwv) fout << setw(25) << "cp ${OUTPUTTRV}.gz" 
			      << " ${SIMULDIR} || OK=0\n";
      if(gin.write.ntwe) fout << setw(25) << "cp ${OUTPUTTRE}.gz" 
			      << " ${SIMULDIR} || OK=0\n";
      if(gin.write.ntwg) fout << setw(25) << "cp ${OUTPUTTRG}.gz" 
			      << " ${SIMULDIR} || OK=0\n";
      // any additional links
      for(unsigned int k=0; k<linkadditions.size(); k++){
	if(linkadditions[k]>0){
	  string s("cp ${"+ linknames[k] + "}");
	  fout << setw(25) << s << " ${SIMULDIR} || OK=0\n";
	}
      }
      
      fout << "\n# clean up after us\n";
      fout << "if `test ${OK} -eq 0`; then\n";
      fout << "  uname -a > mess;\n";
      fout << "  echo 'cp failed for " << systemname << ", run " 
	   << i+scriptNumber << "' >> mess;\n";
      fout << "  Mail -s \"ERROR\" ${NAME} < mess;\n";
      fout << "  cd ${SIMULDIR};\n";
      fout << "else\n";
      fout << "  cd ${SIMULDIR};\n";
      fout << "  rm ${WORKDIR}/*;\n";
      fout << "  rmdir ${WORKDIR};\n";
      fout << "fi\n";

      fout << "\n# perform last command (usually submit next job)\n";
      fout << misc[1].name(i+1) << endl;
      fout.close();
    }


    // Now, some ugly code to change the permissions of the scripts.
    scripts[0]="chmod";
    scripts[1]="u+x";
    char *something[numScripts+3];
    char trialanderror[numScripts+2][100];
    
    for(int i=0; i<numScripts+2; i++){
      strncpy(trialanderror[i], scripts[i].c_str(), 100);
      something[i]=trialanderror[i];
    }
    something[numScripts+2]=NULL;
    execv("/usr/local/bin/chmod", something);
  }
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}

void printWarning(int &numw, int &nume, string s){
    numw++;
    cout << numw+nume << ". WARNING ("<< numw <<")\n";
    cout << s;
    cout << endl;
}
void printError(int &numw, int &nume, string s){
    nume++;
    cout << numw+nume << ". ERROR ("<< nume << ")\n";
    cout << s;
    cout << endl;
}
void printInput(string ifile, string ofile, int nstlim, double t, double dt){
    Ginstream fin(ifile);
    ofstream fout(ofile.c_str());
    string s;
    fout << "TITLE\n";
    fout << fin.title();
    fout << "\nEND\n";
    while(!fin.stream().eof()){
	fin.getline(s);
        if(s=="STEP"){
	    fout << "STEP\n";
	    fout << nstlim << "\t" << t << "\t" << dt << endl;
            fout << "END\n";
            while(s!="END") fin.getline(s);
	}else
	    fout << s << endl;
    }
}

void readLibrary(string file, vector<filename> &names,
		 vector<filename> &misc, 
		 vector<string> &linknames, vector<int> &linkadditions, 
		 string system, string q, string submitcommand, double t, 
		 double dt, int &w, int &e, int ns)
{
  // Open the file
  Ginstream templates(file);
  int found_filenames=0;
  string sdum, temp, first;
  templates.getline(first);
  
  while(!templates.stream().eof()){
    vector<string> buffer;
    templates.getblock(buffer);
    
    if(buffer.size() && first=="FILENAMES"){
      for(unsigned int j=0; j<buffer.size()-1; j++){
	found_filenames=1;
	istringstream iss(buffer[j]);
	iss >> sdum >> temp;
	switch(FILETYPE[sdum]){
	  case inputfile:      names[inputfile].setTemplate(temp);      break;
	  case topofile:       names[topofile].setTemplate(temp);       break;
	  case coordfile:      names[coordfile].setTemplate(temp);      break;
	  case refposfile:     names[refposfile].setTemplate(temp);     break;
	  case posresspecfile: names[posresspecfile].setTemplate(temp); break;
	  case disresfile:     names[disresfile].setTemplate(temp);     break;
	  case pttopofile:     names[pttopofile].setTemplate(temp);     break;
	  case dihresfile:     names[dihresfile].setTemplate(temp);     break;
	  case jvaluefile:     names[jvaluefile].setTemplate(temp);     break;
	  case ledihfile:      names[ledihfile].setTemplate(temp);      break;
	  case outputfile:     names[outputfile].setTemplate(temp);     break;
	  case outtrxfile:     names[outtrxfile].setTemplate(temp);     break;
	  case outtrvfile:     names[outtrvfile].setTemplate(temp);     break;
	  case outtrefile:     names[outtrefile].setTemplate(temp);     break;
	  case outtrgfile:     names[outtrgfile].setTemplate(temp);     break;
	  case scriptfile:     names[scriptfile].setTemplate(temp);     break;
	  case unknownfile:
	    printWarning(w,e, "Don't know how to handle template for "+sdum
			 +". Ingoring");
	}
      }
    }
    if(buffer.size() && first=="MISCELLANEOUS"){
      int l_lastcommand=0;
      for(unsigned int j=0; j<buffer.size()-1; j++){
	istringstream iss(buffer[j]);
	iss >> sdum;
	if(sdum=="workdir") {
	  iss >> temp;
	  misc[0].setTemplate(temp);
	}
	if(sdum=="lastcommand") {
	  l_lastcommand=1;
	  ostringstream os;
	  while(!iss.eof()){
	    iss >> sdum;
	    os << sdum << " ";
	  }
	  misc[1].setTemplate(os.str());
	}
      }
      // re-set the standard lastcommand template, in case the script template
      // has changed
      if(!l_lastcommand)
	misc[1].setTemplate(submitcommand+names[FILETYPE["script"]].temp());
    
    }
    if(buffer.size() && first=="LINKADDITION"){
      for(unsigned int j=0; j<buffer.size()-1; j++){
	istringstream iss(buffer[j]);
	int k;
	string varname;
	
	iss >> sdum >> varname >> temp >> k;
	filename newlink(system, t, dt, ns, q);
	newlink.setTemplate(temp);
	names.push_back(newlink);
	if(sdum=="input") k*=-1;
	linkadditions.push_back(k);
	linknames.push_back(varname);
	
      }
    }
    templates.getline(first);
  }
}
  



