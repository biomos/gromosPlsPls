// mkpert.cc

#include <cassert>
#include <fstream>
#include <sstream>
#include <string>
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <unistd.h>
#include <sys/syscall.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "../src/args/Arguments.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/Solvent.h"
#include "../src/gcore/SolventTopology.h"
#include "../src/gio/InTopology.h"
#include "../src/gio/Ginstream.h"

using namespace std;
using namespace gcore;
using namespace gio;
using namespace args;

#include "mkscript.h"

void printWarning(int &numw, int &nume, string s);
void printError(int &numw, int &nume, string s);
void printInput(string ofile, input gin);
void readLibrary(string file,  vector<filename> &names,
		 vector<filename> &misc, 
		 vector<string> &linknames, vector<int> &linkadditions, 
		 string system, string q, string submitcommand, double t, 
		 double dt, int &w, int &e, int ns);
void readJobinfo(string file, map<int, jobinfo> &ji);
void setParam(input &gin, jobinfo const &job);

int main(int argc, char **argv){

  char *knowns[] = {"sys", "script", "bin", "dir", "queue", 
		    "files", "template", "XX", "cmd", "joblist"} ;
  int nknowns = 10;

  string usage = argv[0];
  usage += "\n\t@sys  <system name>\n";
  usage += "\t@script        <first script> <number of scripts>\n";
  usage += "\t@bin           <gromos96 binary to use>\n";
  usage += "\t[@dir          <where should the files be>]\n";
  usage += "\t[@queue        <which queue?>]\n";
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
  usage += "\t@joblist       <joblist file>\n";
  usage += "\t[@cmd          <last command>\n";
  
  try{
    
    Arguments args(argc, argv, nknowns, knowns, usage);

    // set the number of warnings and the number of errors
    int numWarnings=0;
    int numErrors=0;

    // first get some input parameters
    int scriptNumber=1, numScripts=1;
    string simuldir,q, submitcommand;
    {
      Arguments::const_iterator iter=args.lower_bound("dir");
      if(iter!=args.upper_bound("dir")){
	 simuldir=iter->second;
	 if(simuldir[0]!='/')
	   throw gromos::Exception("mkscript", 
				   "Specified directory should be an "
				   "absolute path");
	 if(chdir(simuldir.c_str())!=0)
	   throw gromos::Exception("mkscript",
				   "Specified directory does not exist\n"
				   "Enter root password:");
      }
      else
	  simuldir="`pwd`";
      q="\"put your favourite queue name here\"";
      if(args.count("queue")>0) q=args["queue"];
      if(q=="penguin")
	submitcommand="psub -s "+q+" 2 ";
      else
	// if(q=="oxen" || q=="moose" || q=="ccpc")
	submitcommand="ssub -s "+q+" ";
      iter=args.lower_bound("script");
      if(iter!=args.upper_bound("script")){
	scriptNumber=atoi(iter->second.c_str());
	++iter;
      }
      if(iter!=args.upper_bound("script"))
	numScripts=atoi(iter->second.c_str());
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
 

    // read the jobinfo file
    if(args.count("joblist")>0 && args.count("script")>0)
      throw gromos::Exception("mkscript", "You can only specify @script OR "
			      "@joblist");
      
    map<int, jobinfo> joblist;
    if(args.count("joblist")>0){
      scriptNumber=0;
      readJobinfo(args["joblist"], joblist);
      map<int, jobinfo>::iterator iter=joblist.begin(), to=joblist.end();
      {
	
	ostringstream os;
	os << gin.step.t;
	iter->second.param["T"]=os.str();
      }
      {
	ostringstream os;
	os << gin.step.t+gin.step.nstlim*gin.step.dt;
	iter->second.param["ENDTIME"]=os.str();
      }
      for(++iter; iter!=to; ++iter){
	if(joblist.find(iter->second.prev_id)!=joblist.end() && 
	   iter->first != iter->second.prev_id){
	  iter->second.param["T"]="-1";
	  iter->second.param["ENDTIME"]="-1";
	}
      }
    }
    else{
      for(int i=0; i<numScripts; i++){
	jobinfo job;
	{
	  ostringstream os;
	  os << gin.step.t+i*gin.step.nstlim*gin.step.dt; 
	  job.param["T"]=os.str();
	}
	{
	  ostringstream os;
	  os << gin.step.t+(i+1)*gin.step.nstlim*gin.step.dt;
	  job.param["ENDTIME"]=os.str();
	}
	job.dir=".";
	job.prev_id=i+scriptNumber-1;
	joblist[i+scriptNumber]=job;
      }
    }
    
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
	s_coord= filenames[FILETYPE["coord"]].name(-1);
	
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

    map<int,jobinfo>::iterator iter=joblist.begin(), to=joblist.end();
    for(; iter!=to; ++iter){
     
      l_coord = l_coord && iter==joblist.begin();
      
      // update the input parameters
      setParam(gin, iter->second);
      {
	ostringstream os;
	os << gin.step.dt*gin.step.nstlim;
	iter->second.param["DELTAT"]=os.str();
      }
      if(iter!=joblist.begin()){
	double time=atof(joblist[iter->second.prev_id].param["ENDTIME"].c_str());
	double endtime=time+gin.step.nstlim*gin.step.dt;
	ostringstream os;
	os << endtime;
	iter->second.param["T"]=joblist[iter->second.prev_id].param["ENDTIME"];
	iter->second.param["ENDTIME"]=os.str();
	gin.step.t=time;
	
      }
      
      for(unsigned int i=0; i<filenames.size(); i++){
	filenames[i].setInfo(systemname, gin.step.t, gin.step.dt*gin.step.nstlim, 
			     iter->first, q);
      }
      for(unsigned int i=0; i<misc.size(); i++){
	misc[i].setInfo(systemname, gin.step.t, gin.step.dt*gin.step.nstlim, 
			iter->first, q);
      }
      
      // Do we go through all the checks?
      if(iter==joblist.begin() || iter->second.param.size() != 3){
	cout << "Performing checks for script " << iter->first << endl;
	cout << "--------------------------------------------" << endl;
	cout << endl;
	
	
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
	if((gin.perturb.found && gin.perturb.ntg != 0) ||
	   (gin.perturb03.found && gin.perturb03.ntg != 0)){

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
	      if (gin.perturb.found)
		os << "NTG = " << gin.perturb.ntg << " in PERTURB block, but "
		   << "no perturbation topology specified\n";
	      else
		os << "NTG = " << gin.perturb03.ntg << " in PERTURB block, but "
		   << "no perturbation topology specified\n";

	      printError(numWarnings, numErrors, os.str());
	    }
	  }

	  double rlamfin;
	  if (gin.perturb.found)
	    rlamfin = gin.perturb.rlam + gin.perturb.dlamt * gin.step.dt * 
	      gin.step.nstlim;
	  else
	    rlamfin = gin.perturb03.rlam + gin.perturb03.dlamt * gin.step.dt * 
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
      }
      
      // Now, write the stupid script
      if(numErrors!=0){
	if(numErrors>1) cout << "\n\nTHERE WERE " << numErrors << " ERRORS\n";
	else cout << "\n\nTHERE WAS 1 ERROR\n";
	cout << "No script will be written\n";
	exit(1);
      }
      if(numErrors==0 && numWarnings==0)
	cout << "OK" << endl << endl;
       

      //first check whether we should be in a different directory
      string subdir=simuldir;
      if(iter->second.dir!="."){
	subdir+="/"+iter->second.dir;
      }
      mkdir(subdir.c_str(), 00755);
      chdir(subdir.c_str());
      cout << "Writing script: " << filenames[FILETYPE["script"]].name(0) << endl;
      cout << "--------------------------------------------" << endl;
      cout << endl;
      
      ofstream fout(filenames[FILETYPE["script"]].name(0).c_str());
      fout.setf(ios_base::left, ios_base::adjustfield);
      fout << "#!/bin/sh" << endl;
      fout << "\n# first we set some variables\n";
      fout << "NAME=`whoami`\n";
      fout << "PROGRAM=" << args["bin"] << endl;
      fout << "SIMULDIR=" <<  simuldir << endl;
      fout << "\n# create temporary directory\n";
      fout << "WORKDIR=" << misc[0].name(0) << endl;
      fout << "mkdir -p ${WORKDIR}\n";
      fout << "cd       ${WORKDIR}\n";
      fout << "\n# set the input files\n";
      fout << "TOPO=${SIMULDIR}/" << s_topo << endl;
      fout << "IUNIT=${SIMULDIR}/";
      if(iter->second.dir!=".") fout << iter->second.dir << "/";
      fout << filenames[FILETYPE["input"]].name(0) << endl;

      if(iter!=joblist.begin() || iter->second.dir!=".")
	printInput(filenames[FILETYPE["input"]].name(0), gin);
	
      fout << "INPUTCRD=${SIMULDIR}/";
      if(iter==joblist.begin()){
	fout << s_coord << endl;
      }
      else{
	jobinfo prevjob=joblist[iter->second.prev_id];
	if(prevjob.dir!=".") fout << prevjob.dir << "/";
	filenames[FILETYPE["coord"]].setInfo(systemname,
					     atof(prevjob.param["T"].c_str()),
					     atof(prevjob.param["DELTAT"].c_str()),
					     iter->second.prev_id,
					     q);
	fout << filenames[FILETYPE["coord"]].name(0) << endl;
	filenames[FILETYPE["coord"]].setInfo(systemname,
					     atof(iter->second.param["T"].c_str()),
					     atof(iter->second.param["DELTAT"].c_str()),
					     iter->first,
					     q);
      }
      
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
	       << filenames[numFiletypes+k].name(0)
	       << endl;
      
      fout << "\n#set the output files\n";
      fout << "OUNIT=" << filenames[FILETYPE["output"]].name(0) << endl;
      fout << "OUTPUTCRD=" << filenames[FILETYPE["coord"]].name(0) << endl;
      if(gin.write.ntwx) fout << "OUTPUTTRX=" 
			      << filenames[FILETYPE["outtrx"]].name(0) 
			      << endl;
      if(gin.write.ntwv) fout << "OUTPUTTRV=" 
			      << filenames[FILETYPE["outtrv"]].name(0)
				<< endl;
      if(gin.write.ntwe) fout << "OUTPUTTRE=" 
			      << filenames[FILETYPE["outtre"]].name(0) 
			      << endl;
      if(gin.write.ntwg) fout << "OUTPUTTRG=" 
			      << filenames[FILETYPE["outtrg"]].name(0) 
			      << endl;
      // any additional links?
      for(unsigned int k=0; k<linkadditions.size(); k++)
	if(linkadditions[k]>0)
	  fout << linknames[k] << "="
	       << filenames[numFiletypes+k].name(0) << endl;
      
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
				 << setw(12) << "@pttopo" << " ${PTTOPO}";
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
      fout << setw(25) << "cp ${OUNIT}" << " ${SIMULDIR}";
      if(iter->second.dir!=".") fout << "/" << iter->second.dir;
      fout << " || OK=0\n";
      fout << setw(25) << "cp ${OUTPUTCRD}" << " ${SIMULDIR}";
      if(iter->second.dir!=".") fout<< "/" << iter->second.dir;
      fout << " || OK=0\n";
      if(gin.write.ntwx){
	fout << setw(25) << "cp ${OUTPUTTRX}.gz" << " ${SIMULDIR}";
	if(iter->second.dir!=".") fout << "/" << iter->second.dir;
	fout << " || OK=0\n";
      }
      if(gin.write.ntwv){
	fout << setw(25) << "cp ${OUTPUTTRV}.gz" << " ${SIMULDIR}";
	if(iter->second.dir!=".") fout << "/" << iter->second.dir;
	fout << " || OK=0\n";
      }
      if(gin.write.ntwe) {
	fout << setw(25) << "cp ${OUTPUTTRE}.gz" << " ${SIMULDIR}";
	if(iter->second.dir!=".") fout << "/" << iter->second.dir;
	fout << " || OK=0\n";
      }
      if(gin.write.ntwg){
	fout << setw(25) << "cp ${OUTPUTTRG}.gz" << " ${SIMULDIR}";
	if(iter->second.dir!=".") fout << "/" << iter->second.dir;
	fout << " || OK=0\n";
      }
      // any additional links
      for(unsigned int k=0; k<linkadditions.size(); k++){
	if(linkadditions[k]>0){
	  string s("cp ${"+ linknames[k] + "}");
	  fout << setw(25) << s << " ${SIMULDIR}";
	  if(iter->second.dir!=".") fout << "/" << iter->second.dir;
	  fout << " || OK=0\n";
	}
      }
      
      fout << "\n# clean up after us\n";
      fout << "if `test ${OK} -eq 0`; then\n";
      fout << "  uname -a > mess;\n";
      fout << "  echo 'cp failed for " << systemname << ", run " 
	   << iter->first << "' >> mess;\n";
      fout << "  Mail -s \"ERROR\" ${NAME} < mess;\n";
      fout << "  cd ${SIMULDIR};\n";
      fout << "else\n";
      fout << "  cd ${SIMULDIR};\n";
      fout << "  rm ${WORKDIR}/*;\n";
      fout << "  rmdir ${WORKDIR};\n";
      fout << "fi\n";
      
      fout << "\n# perform last command (usually submit next job)\n";
      // which job do we have to submit
      map<int, jobinfo>::const_iterator it=iter;
      while (it!=to){
	if(it->second.prev_id==iter->first) {
	  setParam(gin, it->second);
	  misc[1].setInfo(systemname, 
			  atof(iter->second.param["ENDTIME"].c_str()), 
			  gin.step.dt*gin.step.nstlim, it->first, q);
	  if(it->first != iter->first){
	    
	    fout << "cd ${SIMULDIR}";
	    if(it->second.dir!=".") fout << "/" << it->second.dir;
	    fout << "\n";
	    fout << misc[1].name(0) << endl;
	  }
	  
	}
	++it;
      }
      
      fout.close();
      chmod(filenames[FILETYPE["script"]].name(0).c_str(), 00755);
    }
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
void printInput(string ofile, input gin){
    ofstream fout(ofile.c_str());
    const time_t t=time(0);
    fout << "TITLE\n";
    fout << "Automatically generated input file\n";
    fout << getenv("USER") << " " << ctime(&t);
    fout << "END\n";
    fout << gin;
}
void readJobinfo(string file, map<int, jobinfo> &ji)
{
  Ginstream gin(file);
  vector<string> buffer;
  gin.getblock(buffer);
  if(buffer[0]!="JOBSCRIPTS")
    throw gromos::Exception("mkscript", "Reading of jobscript file failed. "
			    "No JOBSCRIPTS block");
  istringstream iss(buffer[1]);
  vector<string> head;
  string b;
  while((iss>>b) !=0) head.push_back(b);
  if(head[0]!="job_id" || head.back() !="run_after" 
     || head[head.size()-2]!="subdir")
    throw gromos::Exception("mkscript", "Reading of jobscript file failed.\n"
			    "First line syntax:\n"
			    "job_id PARAM PARAM ... subdir run_after");
  if(buffer.back().find("END")!=0)
	throw gromos::Exception("mkscript", "Jobscript file " 
				+ gin.name() +
				" is corrupted. No END in JOBSCRIPTS"
				" block. Got\n"
				+ buffer.back());
  int id=0;
  for(unsigned int i=2; i<buffer.size()-1; i++){
    vector<string> tmp(head.size());
    iss.clear();
    iss.str(buffer[i]);
    
    for(unsigned int j=0; j<head.size(); j++) iss >> tmp[j];
    jobinfo job;
    id=atoi(tmp[0].c_str());
    for(unsigned int j=1; j<head.size()-2; j++)
      job.param[head[j]]=tmp[j];
    job.dir=tmp[head.size()-2];
    job.prev_id=atoi(tmp.back().c_str());
    ji[id]=job;
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
      if(buffer[buffer.size()-1].find("END")!=0)
	throw gromos::Exception("mkscript", "Template file " 
				+ templates.name() +
				" is corrupted. No END in "+first+
				" block. Got\n"
				+ buffer[buffer.size()-1]);
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
      if(buffer[buffer.size()-1].find("END")!=0)
	throw gromos::Exception("mkscript", "Template file " + 
				templates.name() +
				" is corrupted. No END in "+first+
				" block. Got\n"
				+ buffer[buffer.size()-1]);
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
      if(buffer[buffer.size()-1].find("END")!=0)
	throw gromos::Exception("mkscript", "Template file " + 
				templates.name() +
				" is corrupted. No END in "+first+
				" block. Got\n"
				+ buffer[buffer.size()-1]);

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

void setParam(input &gin, jobinfo const &job)
{
  map<string, string>::const_iterator iter=job.param.begin(), 
    to=job.param.end();
  for(; iter!=to; ++iter){
    if(iter->first=="NPM")
      gin.system.npm=atoi(iter->second.c_str());
    else if(iter->first=="NSM")
      gin.system.nsm=atoi(iter->second.c_str());
    else if(iter->first=="NTEM")
      gin.minimise.ntem=atoi(iter->second.c_str());
    else if(iter->first=="NMIN")
      gin.minimise.nmin=atoi(iter->second.c_str());
    else if(iter->first=="NTX")
      gin.start.ntx=atoi(iter->second.c_str());
    else if(iter->first=="INIT")
      gin.start.init=atoi(iter->second.c_str());
    else if(iter->first=="IG")
      gin.start.ig=atoi(iter->second.c_str());
    else if(iter->first=="TEMPI")
      gin.start.tempi=atof(iter->second.c_str());
    else if(iter->first=="HEAT")
      gin.start.heat=atof(iter->second.c_str());
    else if(iter->first=="NTXO")
      gin.start.ntx0=atoi(iter->second.c_str());
    else if(iter->first=="BOLTZ")
      gin.start.boltz=atof(iter->second.c_str());
    else if(iter->first=="NSTLIM")
      gin.step.nstlim=atoi(iter->second.c_str());
    else if(iter->first=="T")
      gin.step.t=atof(iter->second.c_str());
    else if(iter->first=="DT")
      gin.step.dt=atof(iter->second.c_str());
    else if(iter->first=="NTB")
      gin.boundary.ntb=atoi(iter->second.c_str());
    else if(iter->first=="NRDBOX")
      gin.boundary.nrdbox=atoi(iter->second.c_str());
    else if(iter->first=="NTT[1]")
      gin.tcouple.ntt[0]=atoi(iter->second.c_str());
    else if(iter->first=="NTT[2]")
      gin.tcouple.ntt[1]=atoi(iter->second.c_str());
    else if(iter->first=="NTT[3]")
      gin.tcouple.ntt[2]=atoi(iter->second.c_str());
    else if(iter->first=="TEMP0[1]")
      gin.tcouple.temp0[0]=atof(iter->second.c_str());
    else if(iter->first=="TEMP0[2]")
      gin.tcouple.temp0[1]=atof(iter->second.c_str());
    else if(iter->first=="TEMP0[3]")
      gin.tcouple.temp0[2]=atof(iter->second.c_str());
    else if(iter->first=="TAUT[1]")
      gin.tcouple.taut[0]=atof(iter->second.c_str());
    else if(iter->first=="TAUT[2]")
      gin.tcouple.taut[1]=atof(iter->second.c_str());
    else if(iter->first=="TAUT[3]")
      gin.tcouple.taut[2]=atof(iter->second.c_str());
    else if(iter->first=="NTP")
      gin.pcouple.ntp=atoi(iter->second.c_str());
    else if(iter->first=="PRES0")
      gin.pcouple.pres0=atof(iter->second.c_str());
    else if(iter->first=="COMP")
      gin.pcouple.comp=atof(iter->second.c_str());
    else if(iter->first=="TAUP")
      gin.pcouple.taup=atof(iter->second.c_str());
    else if(iter->first=="NDFMIN")
      gin.centreofmass.ndfmin=atoi(iter->second.c_str());
    else if(iter->first=="NTCM")
      gin.centreofmass.ntcm=atoi(iter->second.c_str());
    else if(iter->first=="NSCM")
      gin.centreofmass.nscm=atoi(iter->second.c_str());
    else if(iter->first=="NTPR")
      gin.print.ntpr=atoi(iter->second.c_str());
    else if(iter->first=="NTPL")
      gin.print.ntpl=atoi(iter->second.c_str());
    else if(iter->first=="NTPP")
      gin.print.ntpp=atoi(iter->second.c_str());
    else if(iter->first=="NTWX")
      gin.write.ntwx=atoi(iter->second.c_str());
    else if(iter->first=="NTWSE")
      gin.write.ntwse=atoi(iter->second.c_str());
    else if(iter->first=="NTWV")
      gin.write.ntwv=atoi(iter->second.c_str());
    else if(iter->first=="NTWE")
      gin.write.ntwe=atoi(iter->second.c_str());
    else if(iter->first=="NTWG")
      gin.write.ntwg=atoi(iter->second.c_str());
    else if(iter->first=="NTPW")
      gin.write.ntpw=atoi(iter->second.c_str());
    else if(iter->first=="NTC")
      gin.shake.ntc=atoi(iter->second.c_str());
    else if(iter->first=="TOL")
      gin.shake.tol=atof(iter->second.c_str());
    else if(iter->first=="NTNB")
      gin.plist.ntnb=atoi(iter->second.c_str());
    else if(iter->first=="NSNB")
      gin.plist.nsnb=atoi(iter->second.c_str());
    else if(iter->first=="RCUTP")
      gin.plist.rcutp=atof(iter->second.c_str());
    else if(iter->first=="RCUTL")
      gin.plist.rcutl=atof(iter->second.c_str());
    else if(iter->first=="EPSRF")
      gin.longrange.epsrf=atof(iter->second.c_str());
    else if(iter->first=="APPAK")
      gin.longrange.appak=atof(iter->second.c_str());
    else if(iter->first=="RCRF")
      gin.longrange.rcrf=atoi(iter->second.c_str());
    else if(iter->first=="NTR")
      gin.posrest.ntr=atoi(iter->second.c_str());
    else if(iter->first=="CHO")
      gin.posrest.cho=atof(iter->second.c_str());
    else if(iter->first=="NRDRX")
      gin.posrest.nrdrx=atoi(iter->second.c_str());
    else if(iter->first=="NTG"){
      if (gin.perturb.found)
	gin.perturb.ntg=atoi(iter->second.c_str());
      else
	gin.perturb03.ntg=atoi(iter->second.c_str());
    }
    else if(iter->first=="NRDGL")
      gin.perturb.nrdgl=atoi(iter->second.c_str());
    else if(iter->first=="RLAM"){
      if (gin.perturb.found)
	gin.perturb.rlam=atof(iter->second.c_str());
      else
	gin.perturb03.rlam=atof(iter->second.c_str());
    }
    else if(iter->first=="DLAMT"){
      if (gin.perturb.found)
	gin.perturb.dlamt=atof(iter->second.c_str());
      else
	gin.perturb03.dlamt=atof(iter->second.c_str());
    }
    else if(iter->first=="RMU")
      gin.perturb.rmu=atof(iter->second.c_str());
    else if(iter->first=="DMUT")
      gin.perturb.dmut=atof(iter->second.c_str());
    else if(iter->first=="ALPHLJ")
      gin.perturb.alphlj=atof(iter->second.c_str());
    else if(iter->first=="ALPHC")
      gin.perturb.alphc=atof(iter->second.c_str());
    else if(iter->first=="NLAM"){
      if (gin.perturb.found)
	gin.perturb.nlam=atoi(iter->second.c_str());
      else
	gin.perturb03.nlam=atoi(iter->second.c_str());
    }
    else if(iter->first=="MMU")
      gin.perturb.mmu=atoi(iter->second.c_str());
    else if(iter->first=="ENDTIME" || iter->first=="DELTAT"){}
    else
      throw gromos::Exception("mkscript", "Cannot automatically change "
			      +iter->first +" in input file");
  }
}

