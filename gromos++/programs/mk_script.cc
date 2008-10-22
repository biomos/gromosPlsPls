/**
 * @file mk_script.cc
 * generate scripts to run MD simulations
 */
/**
 * @page programs Program Documentation
 *
 * @anchor mk_script
 * @section mk_script generate scripts to run MD simulations
 * @author @ref co
 * @date 10-6-07
 *
 * A molecular dynamics simulation is usually performed by executing a small
 * script that combines all the necessary files and redirects the output to the
 * appropriate places. Especially when simulations are performed on a queue,
 * such scripts become indispensible. In many simulation projects the user is 
 * preparing similar input files and scripts over and over again. Program
 * mk_script can write input files for promd and md++ and the scripts to run
 * the simulation. The format (promd or md++ (or g96)) of the generated md input
 * files will be idenitcal to the format of the md input file given via the 
 * files->input flag. Program mk_script cannot convert program-specific md input
 * blocks into similar input blocks that are specific for another program!
 * Program mk_script can either create a series of similar job scripts that create
 * a sequential list of simulations (keyword \@script) or it can create a more
 * complex set of simulations that are performing a specific task (start-up,
 * perturbation; keyword \@joblist).
 *
 * GROMOS does not require specific filenames for files of specific types.
 * However, most users find it useful to retain some order in their filenames.
 * mk_script has a standard way of constructing filenames that depends on the
 * script-number and the system-code. The user can specify another set of rules
 * to create filenames through the mk_script-library file. In this file, also
 * machine dependent modifications to the scripts that are to be written can be
 * specified. A standard location of the mk_script-library file can be
 * specified through the environment variable MK_SCRIPT_TEMPLATE.
 *
 * By default, mk_script generates job scripts that launches simulation jobs
 * that make use of the Fortran version of the md code (promd). If the c++ program
 * md++ is to be used, the user has to give the @XX (@cpp??) flag.
 *
 * In addition, mk_script performs a number of tests on the input files given
 * to prevent the user from submitting a simulation that fails within the first
 * seconds. In this respect, it distinguishes warnings and errors. A warning is
 * an inconsistency in the input that may lead to an erronous simulation, but 
 * could also be intended. An error is an inconsistency that will for sure lead
 * the simulation to program to crash. Note that this is not a failsafe list of
 * errors, but that it is merely a summary of mistakes that have been made in
 * the past. The tests include:
 * <ol>
 * <li> if the specified binary cannot be found (warning)
 * <li> if no SYSTEM block was given in the input file (error)
 * <li> if the total number of atoms from the topology and the number of
 *      solvent atoms that is specified in the input file does not match the
 *      number of atoms in the coordinate file (error)
 * <li> if no START block was given in the input file (warning)
 * <li> if a VELOCITY block was present in the coordinate file, but NTX in the
 *      START block was set to 1 (warning)
 * <li> if no VELOCITY block was present in the coordinate file, although NTX
 *      in the START block specifies that there should (error)
 * <li> if no STEP block was given in the input file (error)
 * <li> if no BOUNDARY block was given in the input file (error) 
 * <li> if no BOX or TRICLINICBOX block was found in the input coordinate,
 *      although NRDBOX was set to 1 in the BOUNDARY block (error)
 * <li> if the BOX shape was set to truncated octahedral, but the X-, Y- and 
 *      Z-dimensions in the BOX block were not equal (error)
 * <li> if no SUBMOLECULES block was given in the input file (error)
 * <li> if the information in the SUBMOLECULES block does not match the
 *      molecules that are determined from the topology by analysing the bonds
 *      (warning)
 * <li> if the TCOUPLE block in the input file refers to solvent, but there are
 *      no solvent molecules specified (error)
 * <li> if the pressure scaling is requested in the PCOUPLE block, but the
 *      virial is not calculated according to the BOUNDARY block (error)
 * <li> if anisotropic pressure scaling is requested in the PCOUPLE block for a
 *      truncated octahedral box shape (error)
 * <li> if TAUP in the PCOUPLE block is not larger than TAUT in the TCOUPLE
 *      block (warning)
 * <li> if no PRINT block was given in the input file (error)
 * <li> if DT in the STEP block is larger than the suggested values according 
 *      to the SHAKE block (no SHAKE, NTC=1: DT=0.0005; SHAKE bonds with H, 
 *      NTC=2: DT=0.001; SHAKE all bonds, NTC=3: DT=0.002) (warning)
 * <li> if no FORCE block was given in the input file (error)
 * <li> if the FORCE for bonds that are SHAKEn is calculated (warning)
 * <li> if the FORCE for bonds that are not SHAKEn is not calculated (warning)
 * <li> if the last atom of the last energy group is not equal to the total 
 *      number of atoms in the system (error)
 * <li> if no PLIST or PLIST03 block was given in the input file (error)
 * <li> if the short-range cutoff RCUTP is larger than the long-range cutoff 
 *      RCUTL (error)
 * <li> if the shortest edge-to-edge distance of the periodic box is shorter
 *      than twice the long-range cutoff (warning)
 * <li> if no LONGRANGE block was given in the input file (error)
 * <li> if the reaction field cutoff distance RCRF is not equal to the 
 *      long-range cutoff RCUTL (warning)
 * <li> if no REFPOSITION block was found in the coordinate file, even though 
 *      NRDRX was set to 1 in the POSREST block (error)
 * <li> if a reference position file was required, no reference position file 
 *      was specified, but a file with the appropriate name was found instead
 *      (warning)
 * <li> if a reference position file was required but could not be found at all
 *      (error)
 * <li> if the reference position file does not contain a REFPOSITION block
 *      (error)
 * <li> if the total number of atoms in the REFPOSITION block does not match
 *      the total number of atoms
 * <li> if a position restraints specification file was required, no file was
 *      specified but a file with the appropriate name was found instead
 *      (warning)
 * <li> if a position restraints specification file was required but could not 
 *      be found at all (error)
 * <li> if the position restraints specification file does not contain a 
 *      POSRESSPEC (promd) or POSRES (md++) block
 * <li> if a position restraints specification file or a reference position 
 *      file was specified, but position restraints were turned off in the 
 *      input file (warning)
 * <li> if a perturbation topology was required, no file was specified, but a
 *      file with the appropriate name was found instead (warning)
 * <li> if a perturbation topology was required but could not be found at all
 *      (error)
 * <li> if the combination of RLAM, DLAMT in the PERTURB block and the number 
 *      of steps from the STEP block will lead to a lambda value larger than 1
 *      (warning)
 * <li> if a perturbation topologies was specified but no perturbation was 
 *      requested (warning)
 * </ol>
 * In case no errors remain, all requested input files and scripts will be 
 * written to disc. In case of errors, the scripts will not be written, unless 
 * the user forces this by the use of the flag \@force. Scripts for several
 * special cases can be written, such as scripts for REMD simulations, or
 * scripts that run two independent simulations on a single 2 CPU node.
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@sys</td><td>&lt;system name&gt; </td></tr>
 * <tr><td> \@bin</td><td>&lt;gromos96 binary to use&gt; </td></tr>
 * <tr><td> \@dir</td><td>&lt;where should the files be&gt; </td></tr>
 * <tr><td> [\@script</td><td>&lt;first script&gt; &lt;number of scripts&gt;] </td></tr>
 * <tr><td> [\@joblist</td><td>&lt;joblist file&gt;] </td></tr>
 * <tr><td> \@files</td><td></td></tr>
 * <tr><td> [\@template</td><td>&lt;template filenames&gt;] </td></tr>
 * <tr><td> [\@queue</td><td>&lt;which queue?&gt;] </td></tr>
 * <tr><td> [\@XX</td><td>gromosXX script] </td></tr>
 * <tr><td> [\@remd</td><td>&lt;master / slave hostname port&gt; (replica exchange MD)] </td></tr>
 * <tr><td> [\@dual</td><td>&lt;job nr offset&gt; (run two jobs simultaneously)] </td></tr>
 * <tr><td> [\@cmd</td><td>&lt;overwrite last command&gt;] </td></tr>
 * <tr><td> [\@force</td><td>(write script regardless of errors)] </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  mk_script
    @sys         ex
    @bin         /usr/local/gromos/gromos96/bin/promd.64
    @dir         /home/user
    @joblist     joblist.startup
    @files
       topo      ex.top
       coord     exref.coo
       input     imd.dat
    @template    mk_script.lib
    @queue       igcpc
 @endverbatim
 *
 * <hr>
 */

#include <cassert>
#include <fstream>
#include <sstream>
#include <string>
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <algorithm>
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

#include "mk_script.h"

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

  Argument_List knowns;
  knowns << "sys" << "script" << "bin" << "dir" << "queue" << "remd" << "dual"
         << "files" << "template" << "XX" << "cmd" << "joblist"
         << "force";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@sys  <system name>\n";
  usage += "\t@bin           <gromos96 binary to use>\n";
  usage += "\t@dir           <where should the files be>\n";
  usage += "\t[@script       <first script> <number of scripts>]\n";
  usage += "\t[@joblist      <joblist file>]\n";
  usage += "\t@files\n";
  usage += "\t\ttopo         <molecular topology file>\n";
  usage += "\t\tinput        <input file>\n";
  usage += "\t\tcoord        <initial coordinates>\n";
  usage += "\t\t[refpos      <reference positions>]\n";
  usage += "\t\t[posresspec  <position restraints specifications>]\n";
  usage += "\t\t[disres      <distance restraints>]\n";
  usage += "\t\t[dihres      <dihedral restraints>]\n";
  usage += "\t\t[jvalue      <j-value restraints>]\n";
  usage += "\t\t[ledih       <local elevation dihedrals>]\n";
  usage += "\t\t[pttopo      <perturbation topology>]\n";
  usage += "\t[@template     <template filenames>]\n";
  usage += "\t[@queue        <which queue?>]\n"; 
  usage += "\t[@XX           md++ script]\n";
  usage += "\t[@remd         <master / slave hostname port> (replica exchange MD)]\n";
  usage += "\t[@dual         <job nr offset> (run two jobs simultaneously)]\n";
  usage += "\t[@cmd          <overwrite last command>]\n";
  usage += "\t[@force        (write script regardless of errors)]\n";

  try {

    Arguments args(argc, argv, knowns, usage);

    // set the number of warnings and the number of errors
    int numWarnings = 0;
    int numErrors = 0;

    // first get some input parameters
    int scriptNumber = 1, numScripts = 1;
    string simuldir, q, submitcommand;
    {
      Arguments::const_iterator iter = args.lower_bound("dir");
      if (iter != args.upper_bound("dir")) {
        simuldir = iter->second;
        if (simuldir[0] != '/')
          throw gromos::Exception("mk_script",
                "Specified directory should be an "
                "absolute path");
        if (chdir(simuldir.c_str()) != 0)
          throw gromos::Exception("mk_script",
                "Specified directory does not exist\n"
                "Enter root password:");
      } else
        simuldir = "`pwd`";
      q = "\"put your favourite queue name here\"";
      if (args.count("queue") > 0) q = args["queue"];
      if (q == "penguin")
        submitcommand = "psub -s " + q + " 2 ";
      else
        // if(q=="oxen" || q=="moose" || q=="ccpc")
        submitcommand = "ssub -s " + q + " ";
      iter = args.lower_bound("script");
      if (iter != args.upper_bound("script")) {
        scriptNumber = atoi(iter->second.c_str());
        ++iter;
      }
      if (iter != args.upper_bound("script"))
        numScripts = atoi(iter->second.c_str());

      if (numScripts < 0)
        throw Arguments::Exception("Can't deal with negativ number of scripts in @script argument");
    }
    string systemname = args["sys"];

    bool do_remd = false;
    std::string hostname = "";
    int port = -1;

    if (args.count("remd") >= 0) {
      if (args.count("remd") == 0)
        throw Arguments::Exception("remd: expected slave / master");

      do_remd = true;

      Arguments::const_iterator iter = args.lower_bound("remd");
      if (iter != args.upper_bound("remd")) {
        if (iter->second != "slave")
          throw Arguments::Exception("remd: mk_script only for slave");
        ++iter;
      }
      if (iter != args.upper_bound("remd")) {
        hostname = iter->second;
        std::cout << "\ttrying to connect to host " << hostname << "\n";
        ++iter;
      }
      if (iter != args.upper_bound("remd")) {
        std::istringstream is(iter->second);
        if (!(is >> port))
          throw Arguments::Exception("could not read port");
        std::cout << "\ton port " << port << "\n";
      }
    }

    bool dual = false;
    int dual_offset = 1000;
    if (args.count("dual") >= 0)
      dual = true;
    if (args.count("dual") > 0) {
      istringstream is(args["dual"]);
      if (!(is >> dual_offset))
        throw gromos::Exception("mk_script",
              "dual offset wrong");
    }

    // parse the files
    int l_coord = 0, l_topo = 0, l_input = 0, l_refpos = 0, l_posresspec = 0;
    int l_disres = 0, l_dihres = 0, l_jvalue = 0, l_ledih = 0, l_pttopo = 0;
    string s_coord, s_topo, s_input, s_refpos, s_posresspec;
    string s_disres, s_dihres, s_jvalue, s_ledih, s_pttopo;
    for (Arguments::const_iterator iter = args.lower_bound("files"),
            to = args.upper_bound("files"); iter != to; ++iter) {
      switch (FILETYPE[iter->second]) {
        case coordfile: ++iter;
          s_coord = iter->second;
          l_coord = 1;
          break;
        case inputfile: ++iter;
          s_input = iter->second;
          l_input = 1;
          break;
        case topofile: ++iter;
          s_topo = iter->second;
          l_topo = 1;
          break;
        case refposfile: ++iter;
          s_refpos = iter->second;
          l_refpos = 1;
          break;
        case posresspecfile:++iter;
          s_posresspec = iter->second;
          l_posresspec = 1;
          break;
        case disresfile: ++iter;
          s_disres = iter->second;
          l_disres = 1;
          break;
        case dihresfile: ++iter;
          s_dihres = iter->second;
          l_dihres = 1;
          break;
        case jvaluefile: ++iter;
          s_jvalue = iter->second;
          l_jvalue = 1;
          break;
        case ledihfile: ++iter;
          s_ledih = iter->second;
          l_ledih = 1;
          break;
        case pttopofile: ++iter;
          s_pttopo = iter->second;
          l_pttopo = 1;
          break;
        case outputfile: ++iter;
          printWarning(numWarnings, numErrors,
                  iter->second + " not used");
          break;
        case outtrxfile: ++iter;
          printWarning(numWarnings, numErrors,
                  iter->second + " not used");
          break;
        case outtrvfile: ++iter;
          printWarning(numWarnings, numErrors,
                  iter->second + " not used");
          break;
        case outtrefile: ++iter;
          printWarning(numWarnings, numErrors,
                  iter->second + " not used");
          break;
        case outtrgfile: ++iter;
          printWarning(numWarnings, numErrors,
                  iter->second + " not used");
          break;
        case outbaefile: ++iter;
          printWarning(numWarnings, numErrors,
                  iter->second + " not used");
          break;
        case outbagfile: ++iter;
          printWarning(numWarnings, numErrors,
                  iter->second + " not used");
          break;
        case scriptfile: ++iter;
          printWarning(numWarnings, numErrors,
                  iter->second + " not used");
          break;
        case unknownfile: printError(numWarnings, numErrors,
                  "Don't know how to handle file " + iter->second);
      }
    }

    // check which outformat we want (gromos96 or gromosXX)
    // in g08 only relevant for job scripts!
    bool gromosXX = false;
    if (args.count("XX") != -1)
      gromosXX = true;


    // read topology
    if (!l_topo) {
      throw gromos::Exception("mk_script", "You have to specify a topology\n" + usage);
    }
    InTopology it(s_topo);
    System sys = it.system();

    // read the input file
    if (!l_input) {
      throw gromos::Exception("mk_script", "You have to specify an input file\n" + usage);
    }
    Ginstream imd(s_input);

    input gin;
    imd >> gin;

    imd.close();


    // read the jobinfo file
    if (args.count("joblist") > 0 && args.count("script") > 0)
      throw gromos::Exception("mk_script", "You can only specify @script OR "
            "@joblist");

    map<int, jobinfo> joblist;
    if (args.count("joblist") > 0) {
      scriptNumber = 0;
      readJobinfo(args["joblist"], joblist);
      map<int, jobinfo>::iterator iter = joblist.begin(), to = joblist.end();
      {

        ostringstream os;
        os << gin.step.t;
        iter->second.param["T"] = os.str();
      }
      {
        ostringstream os;
        os << gin.step.t + gin.step.nstlim * gin.step.dt;
        iter->second.param["ENDTIME"] = os.str();
      }
      for (++iter; iter != to; ++iter) {
        if (joblist.find(iter->second.prev_id) != joblist.end() &&
            iter->first != iter->second.prev_id) {
          iter->second.param["T"] = "-1";
          iter->second.param["ENDTIME"] = "-1";
        }
      }
    } else { // no joblist
      for (int i = 0; i < numScripts; i++) {
        jobinfo job;
        {
          ostringstream os;
          os << gin.step.t + i * gin.step.nstlim * gin.step.dt;
          job.param["T"] = os.str();
        }
        {
          ostringstream os;
          os << gin.step.t + (i + 1) * gin.step.nstlim * gin.step.dt;
          job.param["ENDTIME"] = os.str();
        }
        job.dir = ".";
        job.prev_id = i + scriptNumber - 1;
        joblist[i + scriptNumber] = job;
      }
    }

    // create names for automated file names
    vector<filename> filenames;
    vector<filename> misc;
    vector<int> linkadditions;
    vector<string> linknames;
    // and a second set for "dual" scripts
    vector<filename> filenames2;
    vector<filename> misc2;
    vector<int> linkadditions2;
    vector<string> linknames2;

    for (int i = 0; i < numFiletypes; i++) {
      filename newname(systemname, gin.step.t, gin.step.nstlim * gin.step.dt,
              scriptNumber, q);
      filenames.push_back(newname);

      if (dual) {
        filename newname(systemname, gin.step.t, gin.step.nstlim * gin.step.dt,
                scriptNumber + dual_offset, q);
        filenames2.push_back(newname);
      }
    }
    // workdir lastcommand firstcommand mpi command
    for (int i = 0; i < 4; i++) {
      filename newname(systemname, gin.step.t, gin.step.nstlim * gin.step.dt,
              scriptNumber, q);
      misc.push_back(newname);
      if (dual) {
        filename newname(systemname, gin.step.t, gin.step.nstlim * gin.step.dt,
                scriptNumber + dual_offset, q);
        misc2.push_back(newname);
      }
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
    filenames[FILETYPE["outbae"]].setTemplate("o%system%baemd_%number%.dat");
    filenames[FILETYPE["outbag"]].setTemplate("o%system%bagmd_%number%.dat");

    misc[0].setTemplate("/scrloc/${NAME}_%system%_%number%");
    misc[1].setTemplate(submitcommand + filenames[FILETYPE["script"]].temp());
    misc[2].setTemplate("");
    misc[3].setTemplate("");

    if (dual) {
      filenames2[FILETYPE["script"]].setTemplate("jmd%system%_%number%.sh");
      filenames2[FILETYPE["input"]].setTemplate("imd%system%_%number%.dat");
      filenames2[FILETYPE["topo"]].setTemplate("%system%mta1.dat");
      filenames2[FILETYPE["refpos"]].setTemplate("%system%px_%number%.dat");
      filenames2[FILETYPE["posresspec"]].setTemplate("%system%pr_%number%.dat");
      filenames2[FILETYPE["disres"]].setTemplate("%system%drts_%number%.dat");
      filenames2[FILETYPE["pttopo"]].setTemplate("%system%pt1.dat");
      filenames2[FILETYPE["dihres"]].setTemplate("%system%arts_%number%.dat");
      filenames2[FILETYPE["jvalue"]].setTemplate("%system%jlts_%number%.dat");
      filenames2[FILETYPE["ledih"]].setTemplate("%system%lets_%number%.dat");
      filenames2[FILETYPE["coord"]].setTemplate("o%system%sxmd_%number%.dat");
      filenames2[FILETYPE["output"]].setTemplate("omd%system%_%number%.out");
      filenames2[FILETYPE["outtrx"]].setTemplate("o%system%trmd_%number%.dat");
      filenames2[FILETYPE["outtrv"]].setTemplate("o%system%tvmd_%number%.dat");
      filenames2[FILETYPE["outtre"]].setTemplate("o%system%temd_%number%.dat");
      filenames2[FILETYPE["outtrg"]].setTemplate("o%system%tgmd_%number%.dat");
      filenames2[FILETYPE["outbae"]].setTemplate("o%system%baemd_%number%.dat");
      filenames2[FILETYPE["outbag"]].setTemplate("o%system%bagmd_%number%.dat");

      misc2[0].setTemplate("/scrloc/${NAME}_%system%_%number%");
      misc2[1].setTemplate(submitcommand + filenames[FILETYPE["script"]].temp());
      misc2[2].setTemplate("");
      misc2[3].setTemplate("");
    }

    // read in the library
    if (args.count("template") >= 0) {
      int really_do_it = 1;
      string libraryfile;
      if (args.count("template") == 0) {
        if (getenv("MK_SCRIPT_TEMPLATE")) {
          libraryfile = getenv("MK_SCRIPT_TEMPLATE");
        } else {
          ostringstream os;
          os << "Trying to read template file, but MK_SCRIPT_TEMPLATE is not set\n"
                  << "Either specify a filename or set this environment variable\n"
                  << "Using defaults now\n";
          printWarning(numWarnings, numErrors, os.str());
          really_do_it = 0;
        }

      }

      if (args.count("template") > 0)
        libraryfile = args["template"];
      if (really_do_it) {
        // And here is a gromos-like function call!
        readLibrary(libraryfile, filenames, misc,
                linknames, linkadditions,
                systemname, q, submitcommand, gin.step.t,
                gin.step.nstlim * gin.step.dt, numWarnings, numErrors,
                scriptNumber);
        // yes, i can make it even worse!
        if (dual)
          readLibrary(libraryfile, filenames2, misc2,
                linknames2, linkadditions2,
                systemname, q, submitcommand, gin.step.t,
                gin.step.nstlim * gin.step.dt, numWarnings, numErrors,
                scriptNumber + dual_offset);
      }
    }

    // overwrite last command if given as argument
    if (args.count("cmd") > 0) {
      ostringstream os;
      Arguments::const_iterator iter = args.lower_bound("cmd"),
              to = args.upper_bound("cmd");
      for (; iter != to; ++iter) {
        std::string s = iter->second;
        if (s.find("\\n") != string::npos)
          s.replace(s.find("\\n"), 2, "\n");
        else s += " ";

        // os << iter->second << " ";
        os << s;

      }

      misc[1].setTemplate(os.str());
      if (dual)
        misc2[1].setTemplate(os.str());
    }

    // read what is in the coordinate file
    fileInfo crd;
    if (!l_coord) {
      // try to open it from the template
      ifstream fin(filenames[FILETYPE["coord"]].name(-1).c_str());
      if (fin) {
        ostringstream os;
        os << "No coordinate file specified, but I found "
                << filenames[FILETYPE["coord"]].name(-1)
                << " which I will use\n";
        printWarning(numWarnings, numErrors, os.str());
        s_coord = filenames[FILETYPE["coord"]].name(-1);
        l_coord = 1;
      } else {
        ostringstream os;
        os << "No coordinate file is specified, some checks are not performed\n";
        os << "Assuming it does not exist yet, I will use "
                << filenames[FILETYPE["coord"]].name(-1)
                << " in the script\n";
        s_coord = filenames[FILETYPE["coord"]].name(-1);

        printWarning(numWarnings, numErrors, os.str());
      }

    }
    if (l_coord) {
      Ginstream icrd(s_coord);
      icrd >> crd;
      icrd.close();
    }

    // calculate some standard numbers
    int numSoluteAtoms = 0;
    for (int i = 0; i < sys.numMolecules(); i++)
      numSoluteAtoms += sys.mol(i).topology().numAtoms();
    int numSolventAtoms = sys.sol(0).topology().numAtoms();
    int numTotalAtoms = gin.system.npm * numSoluteAtoms +
            gin.system.nsm * numSolventAtoms;

    // carry out a thousand tests:

    // Does the binary exist?
    {
      ifstream fin(args["bin"].c_str());
      if (!fin)
        printWarning(numWarnings, numErrors, "Specified binary not found! "
              + args["bin"] + "\n");
      else fin.close();
    }

    map<int, jobinfo>::iterator iter = joblist.begin(), to = joblist.end();
    for (; iter != to; ++iter) {

      //make sure we start in the right directory
      chdir(simuldir.c_str());

      l_coord = l_coord && iter == joblist.begin();

      // update the input parameters
      setParam(gin, iter->second);
      {
        ostringstream os;
        os << gin.step.dt * gin.step.nstlim;
        iter->second.param["DELTAT"] = os.str();
      }
      if (iter != joblist.begin()) {
        double time = atof(joblist[iter->second.prev_id].param["ENDTIME"].c_str());
        double endtime = time + gin.step.nstlim * gin.step.dt;
        ostringstream os;
        os << endtime;
        iter->second.param["T"] = joblist[iter->second.prev_id].param["ENDTIME"];
        iter->second.param["ENDTIME"] = os.str();
        gin.step.t = time;

      }

      for (unsigned int i = 0; i < filenames.size(); i++) {
        filenames[i].setInfo(systemname, gin.step.t, gin.step.dt * gin.step.nstlim,
                iter->first, q);
      }
      for (unsigned int i = 0; i < misc.size(); i++) {
        misc[i].setInfo(systemname, gin.step.t, gin.step.dt * gin.step.nstlim,
                iter->first, q);
      }

      if (dual) {
        if (filenames2.size() != filenames.size())
          cerr << "filenames size does not match!" << std::endl;

        for (unsigned int i = 0; i < filenames2.size(); i++) {
          filenames2[i].setInfo(systemname, gin.step.t, gin.step.dt * gin.step.nstlim,
                  iter->first + dual_offset, q);
        }
        for (unsigned int i = 0; i < misc2.size(); i++) {
          misc2[i].setInfo(systemname, gin.step.t, gin.step.dt * gin.step.nstlim,
                  iter->first + dual_offset, q);
        }
      }

      // Do we go through all the checks?
      if (iter == joblist.begin() || iter->second.param.size() != 3) {
        cout << "Performing checks for script " << iter->first << endl;
        cout << "--------------------------------------------" << endl;
        cout << endl;

        // Ignore gromos96 specific blocks if gromos08 input is to be written
        // and vice versa (start)
        if (args::Arguments::inG96 == false) {
          if (gin.minimise.found) {
            if (gin.energymin.found) {
              printWarning(numWarnings, numErrors,
                      "Ignored GROMOS96 specific block MINIMISE\n");
              gin.minimise.found = 0;
            } else {
              printError(numWarnings, numErrors,
                      "You want to perform a GROMOS08 run and you have specified the GROMOS96 specific\n"
                      "block MINIMISE. Maybe you want to specify the ENERGYMIN block in stead?\n");
            }
          }
          if (gin.stochastic.found) {
            if (gin.stochdyn.found) {
              printWarning(numWarnings, numErrors,
                      "Ignored GROMOS96 specific block STOCHASTIC\n");
              gin.stochastic.found = 0;
            } else {
              printError(numWarnings, numErrors,
                      "You want to perform a GROMOS08 run and you have specified the GROMOS96 specific\n"
                      "block STOCHASTIC. Maybe you want to specify the STOCHDYN block in stead?\n");
            }
          }
          if (gin.start.found) {
            if (gin.initialise.found) {
              printWarning(numWarnings, numErrors,
                      "Ignored GROMOS96 specific block START\n");
              gin.start.found = 0;
            } else {
              printError(numWarnings, numErrors,
                      "You want to perform a GROMOS08 run and you have specified the GROMOS96 specific\n"
                      "block START. Maybe you want to specify the INITIALISE block in stead?\n");
            }
          }
          if (gin.boundary.found) {
            if (gin.boundcond.found) {
              printWarning(numWarnings, numErrors,
                      "Ignored GROMOS96 specific block BOUNDARY\n");
              gin.boundary.found = 0;
            } else {
              printError(numWarnings, numErrors,
                      "You want to perform a GROMOS08 run and you have specified the GROMOS96 specific\n"
                      "block BOUNDARY. Maybe you want to specify the BOUNDCOND block in stead?\n");
            }
          }
          if (gin.submolecules.found) {
            printWarning(numWarnings, numErrors,
                    "Ignored GROMOS96 specific block SUBMOLECULES\n");
            gin.start.found = 0;
          }
          if (gin.tcouple.found) {
            printWarning(numWarnings, numErrors,
                    "Ignored GROMOS96 specific block TCOUPLE\n");
            gin.tcouple.found = 0;
            if (!gin.thermostat.found && !gromosXX) {
              printError(numWarnings, numErrors,
                      "You want to perform a PROMD run and you have specified the GROMOS96 specific\n"
                      "block TCOUPLE. Maybe you want to specify the THERMOSTAT block in stead?\n");
            }
            if (!gin.multibath.found && gromosXX) {
              printError(numWarnings, numErrors,
                      "You want to perform a MD++ run and you have specified the GROMOS96 specific\n"
                      "block TCOUPLE. Maybe you want to specify the MULTIBATH block in stead?\n");
            }
          }
          if (gin.pcouple.found) {
            printWarning(numWarnings, numErrors,
                    "Ignored GROMOS96 specific block PCOUPLE\n");
            gin.pcouple.found = 0;
            if (!gin.barostat.found && !gromosXX) {
              printError(numWarnings, numErrors,
                      "You want to perform a PROMD run and you have specified the GROMOS96 specific\n"
                      "block PCOUPLE. Maybe you want to specify the BAROSTAT block in stead?\n");
            }
            if (!gin.pressurescale.found && gromosXX) {
              printError(numWarnings, numErrors,
                      "You want to perform a MD++ run and you have specified the GROMOS96 specific\n"
                      "block PCOUPLE. Maybe you want to specify the PRESSURESCALE block in stead?\n");
            }
          }
          if (gin.centreofmass.found) {
            printWarning(numWarnings, numErrors,
                    "Ignored GROMOS96 specific block CENTREOFMASS\n");
            gin.centreofmass.found = 0;
            if (!gin.overalltransrot.found && !gromosXX) {
              printError(numWarnings, numErrors,
                      "You want to perform a PROMD run and you have specified the GROMOS96 specific\n"
                      "block CENTREOFMASS. Maybe you want to specify the OVERALLTRANSROT block in stead?\n");
            }
            if (!gin.comtransrot.found && gromosXX) {
              printError(numWarnings, numErrors,
                      "You want to perform a MD++ run and you have specified the GROMOS96 specific\n"
                      "block CENTREOFMASS. Maybe you want to specify the COMTRANSROT block in stead?\n");
            }
          }
          if (gin.print.found) {
            if (gin.printout.found) {
              printWarning(numWarnings, numErrors,
                      "Ignored GROMOS96 specific block PRINT\n");
              gin.print.found = 0;
            } else {
              printError(numWarnings, numErrors,
                      "You want to perform a GROMOS08 run and you have specified the GROMOS96 specific\n"
                      "block PRINT. Maybe you want to specify the PRINTOUT block in stead?\n");
            }
          }
          if (gin.write.found) {
            if (gin.writetraj.found) {
              printWarning(numWarnings, numErrors,
                      "Ignored GROMOS96 specific block WRITE\n");
              gin.write.found = 0;
            } else {
              printError(numWarnings, numErrors,
                      "You want to perform a GROMOS08 run and you have specified the GROMOS96 specific\n"
                      "block WRITE. Maybe you want to specify the WRITETRAJ block in stead?\n");
            }
          }
          if (gin.shake.found) {
            printWarning(numWarnings, numErrors,
                    "Ignored GROMOS96 specific block SHAKE\n");
            gin.shake.found = 0;
            if (!gin.geomconstraint.found && !gromosXX) {
              printError(numWarnings, numErrors,
                      "You want to perform a PROMD run and you have specified the GROMOS96 specific\n"
                      "block SHAKE. Maybe you want to specify the GEOMCONSTRAINT block in stead?\n");
            }
            if (!gin.constraint.found && gromosXX) {
              printError(numWarnings, numErrors,
                      "You want to perform a MD++ run and you have specified the GROMOS96 specific\n"
                      "block SHAKE. Maybe you want to specify the CONSTRAINT block in stead?\n");
            }
          }
          if (gin.plist03.found) {
            printWarning(numWarnings, numErrors,
                    "Ignored old md++ specific block PLIST03\n");
            gin.plist03.found = 0;
          }
          if (gin.plist.found) {
            printWarning(numWarnings, numErrors,
                    "Ignored GROMOS96 specific block PLIST\n");
            gin.plist.found = 0;
            if (!gin.neighbourlist.found && !gromosXX) {
              printError(numWarnings, numErrors,
                      "You want to perform a PROMD run and you have specified the GROMOS96 specific\n"
                      "block PLIST. Maybe you want to specify the NEIGHBOURLIST block in stead?\n");
            }
            if (!gin.pairlist.found && gromosXX) {
              printError(numWarnings, numErrors,
                      "You want to perform a MD++ run and you have specified the GROMOS96 specific\n"
                      "block PLIST. Maybe you want to specify the PAIRLIST block in stead?\n");
            }
          }
          if (gin.posrest.found) {
            if (gin.positionres.found) {
              printWarning(numWarnings, numErrors,
                      "Ignored GROMOS96 specific block POSREST\n");
              gin.posrest.found = 0;
            } else {
              printError(numWarnings, numErrors,
                      "You want to perform a GROMOS08 run and you have specified the GROMOS96 specific\n"
                      "block POSREST. Maybe you want to specify the POSITIONRES block in stead?\n");
            }
          }
          if (gin.distrest.found) {
            if (gin.distanceres.found) {
              printWarning(numWarnings, numErrors,
                      "Ignored GROMOS96 specific block DISTREST\n");
              gin.distrest.found = 0;
            } else {
              printError(numWarnings, numErrors,
                      "You want to perform a GROMOS08 run and you have specified the GROMOS96 specific\n"
                      "block DISTREST. Maybe you want to specify the DISTANCERES block in stead?\n");
            }
          }
          if (gin.diherest.found) {
            if (gin.dihedralres.found) {
              printWarning(numWarnings, numErrors,
                      "Ignored GROMOS96 specific block DIHEREST\n");
              gin.diherest.found = 0;
            } else {
              printError(numWarnings, numErrors,
                      "You want to perform a GROMOS08 run and you have specified the GROMOS96 specific\n"
                      "block DIHEREST. Maybe you want to specify the DIHEDRALRES block in stead?\n");
            }
          }
          if (gin.jval.found) {
            if (gin.jvalueres.found) {
              printWarning(numWarnings, numErrors,
                      "Ignored GROMOS96 specific block J-VAL\n");
              gin.jval.found = 0;
            } else {
              printError(numWarnings, numErrors,
                      "You want to perform a GROMOS08 run and you have specified the GROMOS96 specific\n"
                      "block J-VAL. Maybe you want to specify the JVALUERES block in stead?\n");
            }
          }
          if (gin.localelevation.found) {
            printWarning(numWarnings, numErrors,
                    "Ignored GROMOS96 specific block LOCALELEVATION\n");
            gin.localelevation.found = 0;
            if (!gin.localelev.found && !gromosXX) {
              printError(numWarnings, numErrors,
                      "You want to perform a PROMD run and you have specified the GROMOS96 specific\n"
                      "block LOCALELEVATION. Maybe you want to specify the LOCALELEV block in stead?\n");
            }
          }
      if(gin.perturb.found){
        if(gin.perturbation.found){
          printWarning(numWarnings, numErrors,
            "Ignored GROMOS96 specific block PERTURB\n");
          gin.perturb.found=0;
            } else {
              printError(numWarnings, numErrors,
                      "You want to perform a GROMOS08 run and you have specified the GROMOS96 specific\n"
                      "block PERTURB. Maybe you want to specify the PERTURBATION block in stead?\n");
            }
          }
          if (gin.perturb03.found) {
            printWarning(numWarnings, numErrors,
                    "Ignored old md++ specific block PERTURB03\n");
            gin.perturb03.found = 0;
          }
          if (gin.fourdim.found) {
            printWarning(numWarnings, numErrors,
                    "Ignored GROMOS96 specific block FOURDIM\n");
            gin.fourdim.found = 0;
          }
        } else { // Now we have to ignore GROMOS08 specific blocks!
          if (gin.energymin.found) {
            if (gin.minimise.found) {
              printWarning(numWarnings, numErrors,
                      "Ignored GROMOS08 specific block ENERGYMIN\n");
              gin.energymin.found = 0;
            } else {
              printError(numWarnings, numErrors,
                      "You want to perform a GROMOS96 run and you have specified the GROMOS08 specific\n"
                      "block ENERGYMIN. Maybe you want to specify the MINIMISE block in stead?\n");
            }
          }
          if (gin.stochdyn.found) {
            if (gin.stochastic.found) {
              printWarning(numWarnings, numErrors,
                      "Ignored GROMOS08 specific block STOCHDYN\n");
              gin.stochdyn.found = 0;
            } else {
              printError(numWarnings, numErrors,
                      "You want to perform a GROMOS96 run and you have specified the GROMOS08 specific\n"
                      "block STOCHDYN. Maybe you want to specify the STOCHASTIC block in stead?\n");
            }
          }
          if (gin.initialise.found) {
            if (gin.start.found) {
              printWarning(numWarnings, numErrors,
                      "Ignored GROMOS08 specific block INITIALISE\n");
              gin.initialise.found = 0;
            } else {
              printError(numWarnings, numErrors,
                      "You want to perform a GROMOS96 run and you have specified the GROMOS08 specific\n"
                      "block INITIALISE. Maybe you want to specify the START block in stead?\n");
            }
          }
          if (gin.readtraj.found) {
            printWarning(numWarnings, numErrors,
                    "Ignored GROMOS08 specific block READTRAJ\n");
            gin.energymin.found = 0;
          }
          if (gin.consistencycheck.found) {
            printWarning(numWarnings, numErrors,
                    "Ignored GROMOS08 specific block CONSISTENCYCHECK\n");
            gin.consistencycheck.found = 0;
          }
          if (gin.boundcond.found) {
            if (gin.boundary.found) {
              printWarning(numWarnings, numErrors,
                      "Ignored GROMOS08 specific block BOUNDCOND\n");
              gin.boundcond.found = 0;
            } else {
              printError(numWarnings, numErrors,
                      "You want to perform a GROMOS96 run and you have specified the GROMOS08 specific\n"
                      "block BOUNDCOND. Maybe you want to specify the BOUNDARY block in stead?\n");
            }
          }
          if (gin.multicell.found) {
            printWarning(numWarnings, numErrors,
                    "Ignored GROMOS08 specific block MULTICELL\n");
            gin.multicell.found = 0;
          }
          if (gin.thermostat.found) {
            if (gin.tcouple.found) {
              printWarning(numWarnings, numErrors,
                      "Ignored GROMOS08 specific block THERMOSTAT\n");
              gin.thermostat.found = 0;
            } else {
              printError(numWarnings, numErrors,
                      "You want to perform a GROMOS96 run and you have specified the GROMOS08 specific\n"
                      "block THERMOSTAT. Maybe you want to specify the TCOUPLE block in stead?\n");
            }
          }
          if (gin.multibath.found) {
            if (gin.tcouple.found) {
              printWarning(numWarnings, numErrors,
                      "Ignored GROMOS08 specific block MULTIBATH\n");
              gin.multibath.found = 0;
            } else {
              printError(numWarnings, numErrors,
                      "You want to perform a GROMOS96 run and you have specified the GROMOS08 specific\n"
                      "block MULTIBATH. Maybe you want to specify the TCOUPLE block in stead?\n");
            }
          }
          if (gin.barostat.found) {
            if (gin.pcouple.found) {
              printWarning(numWarnings, numErrors,
                      "Ignored GROMOS08 specific block BAROSTAT\n");
              gin.barostat.found = 0;
            } else {
              printError(numWarnings, numErrors,
                      "You want to perform a GROMOS96 run and you have specified the GROMOS08 specific\n"
                      "block BAROSTAT. Maybe you want to specify the PCOUPLE block in stead?\n");
            }
          }
          if (gin.virial.found) {
            printWarning(numWarnings, numErrors,
                    "Ignored GROMOS08 specific block VIRIAL\n");
            gin.virial.found = 0;
          }
          if (gin.pressurescale.found) {
            if (gin.pcouple.found) {
              printWarning(numWarnings, numErrors,
                      "Ignored GROMOS08 specific block PRESSURESCALE\n");
              gin.pressurescale.found = 0;
            } else {
              printError(numWarnings, numErrors,
                      "You want to perform a GROMOS96 run and you have specified the GROMOS08 specific\n"
                      "block PRESSURESCALE. Maybe you want to specify the PCOUPLE block in stead?\n");
            }
          }
          if (gin.overalltransrot.found) {
            if (gin.centreofmass.found) {
              printWarning(numWarnings, numErrors,
                      "Ignored GROMOS08 specific block OVERALLTRANSROT\n");
              gin.overalltransrot.found = 0;
            } else {
              printError(numWarnings, numErrors,
                      "You want to perform a GROMOS96 run and you have specified the GROMOS08 specific\n"
                      "block OVERALLTRANSROT. Maybe you want to specify the CENTREOFMASS block in stead?\n");
            }
          }
          if (gin.comtransrot.found) {
            if (gin.centreofmass.found) {
              printWarning(numWarnings, numErrors,
                      "Ignored GROMOS08 specific block COMTRANSROT\n");
              gin.comtransrot.found = 0;
            } else {
              printError(numWarnings, numErrors,
                      "You want to perform a GROMOS96 run and you have specified the GROMOS08 specific\n"
                      "block COMTRANSROT. Maybe you want to specify the CENTREOFMASS block in stead?\n");
            }
          }
          if (gin.printout.found) {
            if (gin.print.found) {
              printWarning(numWarnings, numErrors,
                      "Ignored GROMOS08 specific block PRINTOUT\n");
              gin.printout.found = 0;
            } else {
              printError(numWarnings, numErrors,
                      "You want to perform a GROMOS96 run and you have specified the GROMOS08 specific\n"
                      "block PRINTOUT. Maybe you want to specify the PRINT block in stead?\n");
            }
          }
          if (gin.writetraj.found) {
            if (gin.write.found) {
              printWarning(numWarnings, numErrors,
                      "Ignored GROMOS08 specific block WRITETRAJ\n");
              gin.writetraj.found = 0;
            } else {
              printError(numWarnings, numErrors,
                      "You want to perform a GROMOS96 run and you have specified the GROMOS08 specific\n"
                      "block WRITETRAJ. Maybe you want to specify the WRITE block in stead?\n");
            }
          }
          if (gin.ewarn.found) {
            printWarning(numWarnings, numErrors,
                    "Ignored GROMOS08 specific block EWARN\n");
            gin.ewarn.found = 0;
          }
          if (gin.debug.found) {
            printWarning(numWarnings, numErrors,
                    "Ignored GROMOS08 specific block DEBUG\n");
            gin.debug.found = 0;
          }
          if (gin.geomconstraint.found) {
            if (gin.shake.found) {
              printWarning(numWarnings, numErrors,
                      "Ignored GROMOS08 specific block GEOMCONSTRAINTS\n");
              gin.geomconstraint.found = 0;
            } else {
              printError(numWarnings, numErrors,
                      "You want to perform a GROMOS96 run and you have specified the GROMOS08 specific\n"
                      "block GEOMCONSTRAINTS. Maybe you want to specify the SHAKE block in stead?\n");
            }
          }
          if (gin.constraint.found) {
            if (gin.shake.found) {
              printWarning(numWarnings, numErrors,
                      "Ignored GROMOS08 specific block CONSTRAINT\n");
              gin.constraint.found = 0;
            } else {
              printError(numWarnings, numErrors,
                      "You want to perform a GROMOS96 run and you have specified the GROMOS08 specific\n"
                      "block CONSTRAINT. Maybe you want to specify the SHAKE block in stead?\n");
            }
          }
          if (gin.covalentform.found) {
            printWarning(numWarnings, numErrors,
                    "Ignored GROMOS08 specific block COVALENTFORM\n");
            gin.covalentform.found = 0;
          }
          if (gin.neighbourlist.found) {
            if (gin.plist.found) {
              printWarning(numWarnings, numErrors,
                      "Ignored GROMOS08 specific block NEIGHBOURLIST\n");
              gin.neighbourlist.found = 0;
            } else {
              printError(numWarnings, numErrors,
                      "You want to perform a GROMOS96 run and you have specified the GROMOS08 specific\n"
                      "block NEIGHBOURLIST. Maybe you want to specify the PLIST block in stead?\n");
            }
          }
          if (gin.pairlist.found) {
            if (gin.plist.found) {
              printWarning(numWarnings, numErrors,
                      "Ignored GROMOS08 specific block PAIRLIST\n");
              gin.pairlist.found = 0;
            } else {
              printError(numWarnings, numErrors,
                      "You want to perform a GROMOS96 run and you have specified the GROMOS08 specific\n"
                      "block PAIRLIST. Maybe you want to specify the PLIST block in stead?\n");
            }
          }
          if (gin.nonbonded.found) {
            if (gin.longrange.found) {
              printWarning(numWarnings, numErrors,
                      "Ignored GROMOS08 specific block NONBONDED\n");
              gin.nonbonded.found = 0;
            } else {
              printError(numWarnings, numErrors,
                      "You want to perform a GROMOS96 run and you have specified the GROMOS08 specific\n"
                      "block NONBONDED. Maybe you want to specify the LONGRANGE block in stead?\n");
            }
          }
          if (gin.cgrain.found) {
            printWarning(numWarnings, numErrors,
                    "Ignored GROMOS08 specific block CGRAIN\n");
            gin.cgrain.found = 0;
          }
          if (gin.positionres.found) {
            if (gin.posrest.found) {
              printWarning(numWarnings, numErrors,
                      "Ignored GROMOS08 specific block POSITIONRES\n");
              gin.positionres.found = 0;
            } else {
              printError(numWarnings, numErrors,
                      "You want to perform a GROMOS96 run and you have specified the GROMOS08 specific\n"
                      "block POSITIONRES. Maybe you want to specify the POSREST block in stead?\n");
            }
          }
          if (gin.distanceres.found) {
            if (gin.distrest.found) {
              printWarning(numWarnings, numErrors,
                      "Ignored GROMOS08 specific block DISTANCERES\n");
              gin.distanceres.found = 0;
            } else {
              printError(numWarnings, numErrors,
                      "You want to perform a GROMOS96 run and you have specified the GROMOS08 specific\n"
                      "block DISTANCERES. Maybe you want to specify the DISTREST block in stead?\n");
            }
          }
          if (gin.dihedralres.found) {
            if (gin.diherest.found) {
              printWarning(numWarnings, numErrors,
                      "Ignored GROMOS08 specific block DIHEDRALRES\n");
              gin.dihedralres.found = 0;
            } else {
              printError(numWarnings, numErrors,
                      "You want to perform a GROMOS96 run and you have specified the GROMOS08 specific\n"
                      "block DIHEDRALRES. Maybe you want to specify the DIHEREST block in stead?\n");
            }
          }
          if (gin.jvalueres.found) {
            if (gin.jval.found) {
              printWarning(numWarnings, numErrors,
                      "Ignored GROMOS08 specific block JVALUERES\n");
              gin.jvalueres.found = 0;
            } else {
              printError(numWarnings, numErrors,
                      "You want to perform a GROMOS96 run and you have specified the GROMOS08 specific\n"
                      "block JVALUERES. Maybe you want to specify the J-VAL block in stead?\n");
            }
          }
          if (gin.localelev.found) {
            if (gin.localelevation.found) {
              printWarning(numWarnings, numErrors,
                      "Ignored GROMOS08 specific block LOCALELEV\n");
              gin.localelev.found = 0;
            } else {
              printError(numWarnings, numErrors,
                      "You want to perform a GROMOS96 run and you have specified the GROMOS08 specific\n"
                      "block LOCALELEV. Maybe you want to specify the LOCALELEVATION block in stead?\n");
            }
          }
          if (gin.rottrans.found) {
            printWarning(numWarnings, numErrors,
                    "Ignored GROMOS08 specific block ROTTRANS\n");
            gin.rottrans.found = 0;
          }
          if (gin.perturbation.found) {
            if (gin.perturb.found) {
              printWarning(numWarnings, numErrors,
                      "Ignored GROMOS08 specific block PERTURBATION\n");
              gin.perturbation.found = 0;
            } else {
              printError(numWarnings, numErrors,
                      "You want to perform a GROMOS96 run and you have specified the GROMOS08 specific\n"
                      "block PERTURBATION. Maybe you want to specify the PERTURB block in stead?\n");
            }
          }
          if (gin.lambdas.found) {
            printWarning(numWarnings, numErrors,
                    "Ignored GROMOS08 specific block LAMBDAS\n");
            gin.lambdas.found = 0;
          }
          if (gin.umbrella.found) {
            printWarning(numWarnings, numErrors,
                    "Ignored GROMOS08 specific block UMBRELLA\n");
            gin.umbrella.found = 0;
          }
          if (gin.perscale.found) {
            printWarning(numWarnings, numErrors,
                    "Ignored GROMOS08 specific block PERSCALE\n");
            gin.perscale.found = 0;
          }
          if (gin.replica.found) {
            printWarning(numWarnings, numErrors,
                    "Ignored GROMOS08 specific block REPLICA\n");
            gin.replica.found = 0;
          }
          if (gin.innerloop.found) {
            printWarning(numWarnings, numErrors,
                    "Ignored GROMOS08 specific block INNERLOOP\n");
            gin.innerloop.found = 0;
          }
          if (gin.gromos96compat.found) {
            printWarning(numWarnings, numErrors,
                    "Ignored GROMOS08 specific block ENERGYMIN\n");
            gin.energymin.found = 0;
          }
          if (gin.integrate.found) {
            printWarning(numWarnings, numErrors,
                    "Ignored GROMOS08 specific block INTEGRATE\n");
            gin.integrate.found = 0;
          }
          if (gin.randomnumbers.found) {
            printWarning(numWarnings, numErrors,
                    "Ignored GROMOS08 specific block RANDOMNUMBERS\n");
            gin.energymin.found = 0;
          }
        }
        // Ignore gromos96 specific blocks if gromos08 input is to be written
        // and vice versa (end)

        // Ignore md++ specific blocks if promd input is to be written
        // and vice versa (start)
        if (args::Arguments::inG96 == true) {
          if (!gromosXX) {
            if (gin.plist03.found) {
              printWarning(numWarnings, numErrors,
                      "Ignored old md++ specific block PLIST03\n");
              gin.plist03.found = 0;
            }
            if (gin.perturb03.found) {
              printWarning(numWarnings, numErrors,
                      "Ignored old md++ specific block PERTURB03\n");
              gin.perturb03.found = 0;
            }
          }
        } else {
          if (!gromosXX) { // Ignore md++ specific blocks
            if (gin.multibath.found) {
              if (gin.thermostat.found) {
                printWarning(numWarnings, numErrors,
                        "Ignored md++ specific block MULTIBATH\n");
                gin.multibath.found = 0;
              } else {
                printError(numWarnings, numErrors,
                        "You want to perform a PROMD run and you have specified the md++ specific\n"
                        "block MULTIBATH. Maybe you want to specify the THERMOSTAT block in stead?\n");
              }
            }
            if (gin.pressurescale.found) {
              if (gin.barostat.found) {
                printWarning(numWarnings, numErrors,
                        "Ignored md++ specific block PRESSURESCALE\n");
                gin.pressurescale.found = 0;
              } else {
                printError(numWarnings, numErrors,
                        "You want to perform a PROMD run and you have specified the md++ specific\n"
                        "block PRESSURESCALE. Maybe you want to specify the BAROSTAT block in stead?\n");
              }
            }
            if (gin.comtransrot.found) {
              if (gin.overalltransrot.found) {
                printWarning(numWarnings, numErrors,
                        "Ignored md++ specific block COMTRANSROT\n");
                gin.comtransrot.found = 0;
              } else {
                printError(numWarnings, numErrors,
                        "You want to perform a PROMD run and you have specified the md++ specific\n"
                        "block COMTRANSROT. Maybe you want to specify the OVERALLTRANSROT block in stead?\n");
              }
            }
            if (gin.ewarn.found) {
              printWarning(numWarnings, numErrors,
                      "Ignored md++ specific block EWARN\n");
              gin.ewarn.found = 0;
            }
            if (gin.constraint.found) {
              if (gin.geomconstraint.found) {
                printWarning(numWarnings, numErrors,
                        "Ignored md++ specific block CONSTRAINT\n");
                gin.constraint.found = 0;
              } else {
                printError(numWarnings, numErrors,
                        "You want to perform a PROMD run and you have specified the md++ specific\n"
                        "block CONSTRAINT. Maybe you want to specify the GEOMCONSTRAINT block in stead?\n");
              }
            }
            if (gin.pairlist.found) {
              if (gin.neighbourlist.found) {
                printWarning(numWarnings, numErrors,
                        "Ignored md++ specific block PAIRLIST\n");
                gin.pairlist.found = 0;
              } else {
                printError(numWarnings, numErrors,
                        "You want to perform a PROMD run and you have specified the md++ specific\n"
                        "block PAIRLIST. Maybe you want to specify the NEIGHBOURLIST block in stead?\n");
              }
            }
            if (gin.longrange.found) {
              if (gin.nonbonded.found) {
                printWarning(numWarnings, numErrors,
                        "Ignored md++ and GROMOS96 specific block LONGRANGE\n");
                gin.longrange.found = 0;
              } else {
                printError(numWarnings, numErrors,
                        "You want to perform a PROMD run and you have specified the md++ specific\n"
                        "block LONGRANGE. Maybe you want to specify the NONBONDED block in stead?\n");
              }
            }
            if (gin.cgrain.found) {
              printWarning(numWarnings, numErrors,
                      "Ignored md++ specific block CGRAIN\n");
              gin.cgrain.found = 0;
            }
            if (gin.rottrans.found) {
              printWarning(numWarnings, numErrors,
                      "Ignored md++ specific block ROTTRANS\n");
              gin.rottrans.found = 0;
            }
            if (gin.lambdas.found) {
              printWarning(numWarnings, numErrors,
                      "Ignored md++ specific block LAMBDAS\n");
              gin.lambdas.found = 0;
            }
            if (gin.perscale.found) {
              printWarning(numWarnings, numErrors,
                      "Ignored md++ specific block PERSCALE\n");
              gin.perscale.found = 0;
            }
            if (gin.replica.found) {
              printWarning(numWarnings, numErrors,
                      "Ignored md++ specific block REPLICA\n");
              gin.replica.found = 0;
            }
            if (gin.innerloop.found) {
              printWarning(numWarnings, numErrors,
                      "Ignored md++ specific block INNERLOOP\n");
              gin.innerloop.found = 0;
            }
            if (gin.integrate.found) {
              printWarning(numWarnings, numErrors,
                      "Ignored md++ specific block INTEGRATE\n");
              gin.integrate.found = 0;
            }
            if (gin.randomnumbers.found) {
              printWarning(numWarnings, numErrors,
                      "Ignored md++ specific block RANDOMNUMBERS\n");
              gin.randomnumbers.found = 0;
            }
          } else { // Ignore promd specific blocks
            if (gin.consistencycheck.found) {
              printWarning(numWarnings, numErrors,
                      "Ignored promd specific block CONSISTENCYCHECK\n");
              gin.consistencycheck.found = 0;
            }
            if (gin.thermostat.found) {
              if (gin.multibath.found) {
                printWarning(numWarnings, numErrors,
                        "Ignored promd specific block THERMOSTAT\n");
                gin.thermostat.found = 0;
              } else {
                printError(numWarnings, numErrors,
                        "You want to perform a md++ run and you have specified the PROMD specific\n"
                        "block THERMOSTAT. Maybe you want to specify the MULTIBATH block in stead?\n");
              }
            }
            if (gin.barostat.found) {
              if (gin.pressurescale.found) {
                printWarning(numWarnings, numErrors,
                        "Ignored promd specific block BAROSTAT\n");
                gin.barostat.found = 0;
              } else {
                printError(numWarnings, numErrors,
                        "You want to perform a md++ run and you have specified the PROMD specific\n"
                        "block BAROSTAT. Maybe you want to specify the PRESSURESCALE block in stead?\n");
              }
            }
            if (gin.virial.found) {
              printWarning(numWarnings, numErrors,
                      "Ignored promd specific block VIRIAL\n");
              gin.virial.found = 0;
            }
            if (gin.overalltransrot.found) {
              if (gin.comtransrot.found) {
                printWarning(numWarnings, numErrors,
                        "Ignored promd specific block OVERALLTRANSROT\n");
                gin.overalltransrot.found = 0;
              } else {
                printError(numWarnings, numErrors,
                        "You want to perform a md++ run and you have specified the PROMD specific\n"
                        "block OVERALLTRANSROT. Maybe you want to specify the COMTRANSROT block in stead?\n");
              }
            }
            if (gin.debug.found) {
              printWarning(numWarnings, numErrors,
                      "Ignored promd specific block DEBUG\n");
              gin.debug.found = 0;
            }
            if (gin.geomconstraint.found) {
              if (gin.constraint.found) {
                printWarning(numWarnings, numErrors,
                        "Ignored promd specific block GEOMCONSTRAINTS\n");
                gin.geomconstraint.found = 0;
              } else {
                printError(numWarnings, numErrors,
                        "You want to perform a md++ run and you have specified the PROMD specific\n"
                        "block GEOMCONSTRAINT. Maybe you want to specify the CONSTRAINT block in stead?\n");
              }
            }
            if (gin.neighbourlist.found) {
              if (gin.pairlist.found) {
                printWarning(numWarnings, numErrors,
                        "Ignored promd specific block NEIGHBOURLIST\n");
                gin.neighbourlist.found = 0;
              } else {
                printError(numWarnings, numErrors,
                        "You want to perform a md++ run and you have specified the PROMD specific\n"
                        "block NEIGHBOURLIST. Maybe you want to specify the PAIRLIST block in stead?\n");
              }
            }
            if (gin.nonbonded.found) {
              if (gin.longrange.found) {
                printWarning(numWarnings, numErrors,
                        "Ignored promd specific block NONBONDED\n");
                gin.nonbonded.found = 0;
              } else {
                printError(numWarnings, numErrors,
                        "You want to perform a md++ run and you have specified the PROMD specific\n"
                        "block NONBONDED. Maybe you want to specify the LONGRANGE block in stead?\n");
              }
            }
            if (gin.localelev.found) {
              printWarning(numWarnings, numErrors,
                      "Ignored promd specific block LOCALELEV\n");
              gin.localelev.found = 0;
            }
            if (gin.umbrella.found) {
              printWarning(numWarnings, numErrors,
                      "Ignored promd specific block UMBRELLA\n");
              gin.umbrella.found = 0;
            }
            if (gin.gromos96compat.found) {
              printWarning(numWarnings, numErrors,
                      "Ignored promd specific block GROMOS96COMPAT\n");
              gin.gromos96compat.found = 0;
            }
            if (gin.pathint.found) {
              printWarning(numWarnings, numErrors,
                      "Ignored promd specific block PATHINT\n");
              gin.pathint.found = 0;
            }
          }
        }
        // Ignore md++ specific blocks if promd input is to be written
        // and vice versa (end)

        // And check if all compulsory blocks have been specified
        if (!gin.system.found) {
          printError(numWarnings, numErrors, "Could not find SYSTEM block\n");
        }
        if (!gin.start.found && args::Arguments::inG96 == true) {
          printError(numWarnings, numErrors, "Could not find START block\n");
        }
        if (!gin.step.found) {
          printError(numWarnings, numErrors, "Could not find STEP block\n");
        }
        if (!gin.boundary.found && args::Arguments::inG96 == true) {
          printError(numWarnings, numErrors, "Could not find BOUNDARY block\n");
        }
        if (!gin.boundcond.found && args::Arguments::inG96 == false) {
          printError(numWarnings, numErrors, "Could not find BOUNDCOND block\n");
        }
        if (!gin.submolecules.found && args::Arguments::inG96 == true) {
          printError(numWarnings, numErrors, "Could not find SUBMOLECULES block\n");
        }
        if (!gin.print.found && args::Arguments::inG96 == true) {
          printError(numWarnings, numErrors, "Could not find PRINT block\n");
        }
        if (!gin.force.found) {
          printError(numWarnings, numErrors, "Could not find FORCE block\n");
        }
        if (args::Arguments::inG96 == true) {
          if (!gin.plist.found) {
            printError(numWarnings, numErrors, "Could not find PLIST block\n");
          }
        } else {
          if (!gromosXX && !gin.neighbourlist.found) {
            printError(numWarnings, numErrors, "Could not find NEIGHBOURLIST block\n");
          }
          if (gromosXX && !gin.pairlist.found) {
            printError(numWarnings, numErrors, "Could not find PAIRLIST block\n");
          }
        }
        if (args::Arguments::inG96 == false && !gromosXX) {
          if (!gin.nonbonded.found) {
            printError(numWarnings, numErrors, "Could not find NONBONDED block\n");
          }
        } else if (!gin.longrange.found) {
          printError(numWarnings, numErrors, "Could not find LONGRANGE block\n");
        }
        if (args::Arguments::inG96 == false) {
          if (!gromosXX && !gin.geomconstraint.found) {
            printError(numWarnings, numErrors, "Could not find GEOMCONSTRAINT block\n");
          }
          if (gromosXX && !gin.constraint.found) {
            printError(numWarnings, numErrors, "Could not find CONSTRAINT block\n");
          }
        }

        //SYSTEM block
        if (gin.system.found) {
          if (l_coord) {
            for (unsigned int i = 0; i < crd.blocks.size(); i++) {
              if (crd.blocks[i] == "POSITION" ||
                  crd.blocks[i] == "VELOCITY" ||
                  crd.blocks[i] == "REFPOSITION" ||
                  crd.blocks[i] == "REDPOSITION") {
                if (numTotalAtoms != crd.blockslength[i]) {
                  ostringstream os;
                  os << "From topology and SYSTEM block, I calculate "
                          << numTotalAtoms << " atoms in the system" << endl;
                  os << "But coordinate file has " << crd.blockslength[i]
                          << " lines in " << crd.blocks[i] << " block." << endl;
                  os << "Maybe NSM should be "
                          << (crd.blockslength[i] - numSoluteAtoms) / numSolventAtoms
                          << " ?" << endl;
                  printError(numWarnings, numErrors, os.str());
                }
              }
            }
          }
        }

        // START
        if (gin.start.found) {
          if (args::Arguments::inG96 == false) {
            printWarning(numWarnings, numErrors,
                    "Ignored GROMOS96 specific block START\n");
            gin.start.found = 0;
          } else if (l_coord) {
            int velblock = 0;
            for (unsigned int i = 0; i < crd.blocks.size(); i++)
              if (crd.blocks[i] == "VELOCITY") velblock = 1;
            if (velblock && gin.start.ntx == 1) {
              ostringstream os;
              os << "NTX = 1 in START block, which means that you don't want to "
                      << "read in the\n";
              os << "velocity block I found in the coordinate file.\n";
              printWarning(numWarnings, numErrors, os.str());
            }
            if (gin.start.ntx > 1 && !velblock) {
              ostringstream os;
              os << "NTX = " << gin.start.ntx << " in START block means that you "
                      << "want to read in the velocities from the coordinate file.\n";
              os << "But coordinate file does not have a VELOCITY block.\n";
              printError(numWarnings, numErrors, os.str());
            }
          }
        }

        // BOUNDARY
        double box[3];
        if (gin.boundary.found) {
          int boxblock = 0;
          if (gin.boundary.nrdbox && l_coord) {
            for (unsigned int i = 0; i < crd.blocks.size(); i++)
              if (crd.blocks[i] == "BOX" ||
                  crd.blocks[i] == "TRICLINICBOX") boxblock = 1;
            if (!boxblock) {
              string s = "NRDBOX = 1 in BOUNDARY block, ";
              s += "but no BOX block in coordinate file.\n";
              printError(numWarnings, numErrors, s);
            } else
              for (int i = 0; i < 3; i++)
                box[i] = crd.box[0];
          }
          if (gin.boundary.nrdbox == 0)
            for (int i = 0; i < 3; i++)
              box[i] = gin.boundary.box[i];
          if (gin.boundary.ntb < 0 && (gin.boundary.nrdbox == 0 || boxblock)) {
            if (box[0] != box[1] || box[0] != box[2] || box[1] != box[2]) {
              ostringstream os;
              os << "NTB = " << gin.boundary.ntb << " in BOUNDARY block "
                      << "means truncated octahedron.\n";
              os << "But boxdimensions (from ";
              if (boxblock) os << "coordinate file ";
              else os << "input file ";
              os << ") are not the same:\n";
              os << box[0] << "\t" << box[1] << "\t" << box[2] << endl;
              printError(numWarnings, numErrors, os.str());
            }
          }
        }

        // SUBMOLECULES
        if (gin.submolecules.found) {
          int na = 0;
          int error = 0;
          if (int(gin.submolecules.nsp.size()) != sys.numMolecules()) error = 1;
          else
            for (int i = 0; i < sys.numMolecules(); i++) {
              na += sys.mol(i).topology().numAtoms();
              if (na != gin.submolecules.nsp[i]) error = 1;
            }
          if (error) {
            ostringstream os;
            na = 0;
            os << "SUBMOLECULES block does not match topology; should be\n";
            os << "SUBMOLECULES\n";
            os << "#     NSPM  NSP(1.. NSPM)\n";
            os << setw(10) << sys.numMolecules() << endl;
            int countsbm = 0;

            for (int i = 0; i < sys.numMolecules(); i++) {
              na += sys.mol(i).topology().numAtoms();
              os << setw(6) << na;
              countsbm++;
              if (countsbm % 10 == 0) os << endl;

            }
            os << "\nEND\n";
            printWarning(numWarnings, numErrors, os.str());
          }
        }

        // TCOUPLE
        if (gin.tcouple.found) {
          if (gin.system.nsm == 0 && gin.tcouple.ntt[2] != 0) {
            string s = "There are no solvent ";
            s += "molecules specified in the SYSTEM block, but\n";
            s += "you do couple its temperature in the TCOUPLE ";
            s += "block (third line)\n";
            printError(numWarnings, numErrors, s);
          }
        }

        // PCOUPLE
        if (gin.pcouple.found) {
          if (gin.pcouple.ntp != 0 && abs(gin.boundary.ntb) != 2) {
            ostringstream os;
            int ntb;
            if (gin.boundary.ntb < 0) ntb = -2;
            else ntb = 2;
            os << "NTP = " << gin.pcouple.ntp << " in PCOUPLE block, but "
                    << "NTB = " << gin.boundary.ntb << " in BOUNDARY block, which "
                    << "means that we do not calculate the virial\n";
            os << "Set NTB = " << ntb << " to calculate the virial\n";
            printError(numWarnings, numErrors, os.str());
          }
          if (gin.pcouple.ntp == 2 && gin.boundary.ntb < 0) {
            ostringstream os;
            os << "NTP = " << gin.pcouple.ntp << " in PCOUPLE block specifies "
                    << "anisotropic pressure scaling.\n";
            os << "But NTB = " << gin.boundary.ntb << " in BOUNDARY block means "
                    << "truncated octahedral boundary conditions\n";
            printError(numWarnings, numErrors, os.str());
          }
          if (gin.tcouple.found) {
            double tautmax = 0;
            for (int i = 0; i < 3; i++)
              if (gin.tcouple.ntt[i] != 0 && gin.tcouple.taut[i] > tautmax)
                tautmax = gin.tcouple.taut[i];
            if (tautmax >= gin.pcouple.taup) {
              ostringstream os;
              os << "TAUP in PCOUPLE block (" << gin.pcouple.taup << ") is not "
                      << "larger than TAUT in TCOUPLE block (" << tautmax << ").\n";
              printWarning(numWarnings, numErrors, os.str());
            }
          }
        }

        //CENTROFMASS

        //WRITE

        //SHAKE
        if (((!gin.shake.found || gin.shake.ntc == 1) && gin.step.dt > 0.0005) ||
            ((gin.shake.found && gin.shake.ntc == 2) && gin.step.dt > 0.001) ||
            ((gin.shake.found && gin.shake.ntc == 3) && gin.step.dt > 0.002)) {
          ostringstream os;
          string comment;
          double suggest = 0.0005;
          if (!gin.shake.found) {
            comment = "no SHAKE block: no shake on solute";
            suggest = 0.0005;
            gin.shake.ntc = 1;
          } else if (gin.shake.ntc == 1) {
            comment = "no shake on solute";
            suggest = 0.0005;
          } else if (gin.shake.ntc == 2) {
            comment = "shake bonds with H";
            suggest = 0.001;
          } else if (gin.shake.ntc == 3) {
            comment = "shake all bonds";
            suggest = 0.002;
          }
          os << "DT in STEP block is set to " << gin.step.dt << ", which is "
                  << "considered to be too large\n";
          os << "if NTC = " << gin.shake.ntc << " in SHAKE block.\n";
          os << "For NTC = " << gin.shake.ntc << " (" << comment << ") rather "
                  << "use DT = " << suggest << ".\n";
          printWarning(numWarnings, numErrors, os.str());
        }

        // FORCE
        if (gin.force.found) {
          string comment;
          if (gin.shake.ntc == 2) comment = "shake bonds with H";
          if (gin.shake.ntc == 3) comment = "shake all bonds";

          if ((gin.shake.ntc > 1 && gin.force.ntf[0]) ||
              (gin.shake.ntc == 3 && gin.force.ntf[1])) {
            ostringstream os;
            os << "NTC = " << gin.shake.ntc << " in SHAKE block ("
                    << comment << ")\n";
            os << "so there is no need to calculate the forces due to these "
                    << "bonds, as indicated in the FORCE block.\n";
            printWarning(numWarnings, numErrors, os.str());
          }
          if ((gin.shake.ntc == 1 && (gin.force.ntf[0] == 0 || gin.force.ntf[1] == 0)) ||
              (gin.shake.ntc == 2 && gin.force.ntf[1] == 0)) {
            ostringstream os;
            os << "NTC = " << gin.shake.ntc << " in SHAKE block (";
            os << comment << ")\n";
            os << "But you do not want to calculate the forces on non-shaken "
                    << "bonds?\n";
            os << "ntf[0] = " << gin.force.ntf[0] << " ntf[1] = "
                    << gin.force.ntf[1] << endl;
            printWarning(numWarnings, numErrors, os.str());
          }
          if (gin.force.nre[gin.force.nre.size() - 1] != numTotalAtoms) {
            ostringstream os;
            os << "The last energy group in the FORCE block ("
                    << gin.force.nre[gin.force.nre.size() - 1] << ") should be equal "
                    << "to the total number of atoms in the system: "
                    << numTotalAtoms << ".\n";
            printError(numWarnings, numErrors, os.str());
          }
        }

        //PLIST
        if (gin.plist.found) {
          if (gin.plist.rcutp > gin.plist.rcutl) {
            ostringstream os;
            os << "In PLIST block RCUTP = " << gin.plist.rcutp << " and RCUTL = "
                    << gin.plist.rcutl << endl;
            os << "RCUTP should be less than or equal to RCUTL\n";
            printError(numWarnings, numErrors, os.str());
          }
          double minbox = 1e6;
          int trunc = 0;
          if (gin.boundary.ntb != 0) {

            if (gin.boundary.nrdbox == 0 || (gin.boundary.nrdbox == 1 && l_coord))
              for (int i = 0; i < 3; i++) if (box[i] < minbox) minbox = box[i];
            if (gin.boundary.ntb < 0) {
              trunc = 1;
              // check formula
              minbox *= 0.5 * 1.732051;
            }
            if (minbox < 2 * gin.plist.rcutl) {
              ostringstream os;
              os << "RCUTL in PLIST block is " << gin.plist.rcutl << endl;
              os << "for a ";
              if (trunc) os << "truncated octahedral ";
              else os << "rectangular ";
              os << "box with these dimensions\n";
              for (int i = 0; i < 3; i++) os << "\t" << box[i];
              if (gin.boundary.nrdbox == 1) os << "\t(from coordinate file)\n";
              else os << "\t(from input file)\n";
              os << "this is too long\n";
              printWarning(numWarnings, numErrors, os.str());
            }
          }
        }

        if (gin.plist03.found) { // PLIST03 (old gromosXX format)
          if (gin.plist03.rcutp > gin.plist03.rcutl) {
            ostringstream os;
            os << "In PLIST03 block RCUTP = " << gin.plist03.rcutp << " and RCUTL = "
                    << gin.plist03.rcutl << endl;
            os << "RCUTP should be less than or equal to RCUTL\n";
            printError(numWarnings, numErrors, os.str());
          }
          double minbox = 1e6;
          int trunc = 0;
          if (gin.boundary.ntb != 0) {

            if (gin.boundary.nrdbox == 0 || (gin.boundary.nrdbox == 1 && l_coord))
              for (int i = 0; i < 3; i++) if (box[i] < minbox) minbox = box[i];
            if (gin.boundary.ntb < 0) {
              trunc = 1;
              // check formula
              minbox *= 0.5 * 1.732051;
            }
            if (minbox < 2 * gin.plist03.rcutl) {
              ostringstream os;
              os << "RCUTL in PLIST03 block is " << gin.plist03.rcutl << endl;
              os << "for a ";
              if (trunc) os << "truncated octahedral ";
              else os << "rectangular ";
              os << "box with these dimensions\n";
              for (int i = 0; i < 3; i++) os << "\t" << box[i];
              if (gin.boundary.nrdbox == 1) os << "\t(from coordinate file)\n";
              else os << "\t(from input file)\n";
              os << "this is too long\n";
              printWarning(numWarnings, numErrors, os.str());
            }
          }
        }

        //LONGRANGE
        if (gin.longrange.found) {
          double rcutl = 0;
          if (gin.plist.found) rcutl = gin.plist.rcutl;
          else if (gin.plist03.found) rcutl = gin.plist03.rcutl;

          if ((gin.longrange.epsrf != 1.0) && (rcutl != fabs(gin.longrange.rcrf))) {
            ostringstream os;
            os << "We usually expect RCRF in the LONGRANGE block to be equal to "
                    << "RCUTL in the PLIST block\n";
            os << "You specified RCUTL = " << rcutl << " and RCRF = "
                    << gin.longrange.rcrf << ".\n";
            printWarning(numWarnings, numErrors, os.str());
          }
        }

        //POSREST
        if (gin.posrest.found && gin.posrest.ntr != 0) {
          int refblock = 0;
          int numref = 0;

          if (!gromosXX) {

            if (gin.posrest.nrdrx == 1 && l_coord) {
              for (unsigned int i = 0; i < crd.blocks.size(); i++) {
                if (crd.blocks[i] == "REFPOSITION") {
                  refblock = 1;
                  numref = crd.blockslength[i];
                }
              }
              if (!refblock) {
                string s = "NRDRX=1 in POSREST block, but no REFPOSITION block ";
                s += "in coordinate file\n";
                printError(numWarnings, numErrors, s);
              }
            } else {
              //first find out if the file exists
              if (!l_refpos && !gin.posrest.nrdrx) {
                if (!l_refpos) {
                  ifstream fin(filenames[refposfile].name(0).c_str());
                  if (fin) {
                    ostringstream os;
                    os << "No refpos-file specified, but I found "
                            << filenames[refposfile].name(0)
                            << " which I will use\n";
                    l_refpos = 1;
                    s_refpos = filenames[refposfile].name(0);
                    printWarning(numWarnings, numErrors, os.str());
                    fin.close();
                  } else {
                    ostringstream os;
                    os << "NRDRX = " << gin.posrest.nrdrx << " in POSREST block, but no "
                            << "refpos-file specified\n";
                    printError(numWarnings, numErrors, os.str());
                  }
                }
              } else if (l_refpos) {
                Ginstream irfp(s_refpos);
                fileInfo rfp;
                irfp >> rfp;
                irfp.close();
                for (unsigned int i = 0; i < rfp.blocks.size(); i++) {
                  if (rfp.blocks[i] == "REFPOSITION") {
                    refblock = 1;
                    numref = rfp.blockslength[i];
                  }
                }
                if (!refblock) {
                  ostringstream os;
                  os << "No REFPOSITION block in refpos file (" << s_refpos
                          << ")\n";
                  printError(numWarnings, numErrors, os.str());
                }
              }
            }

            if (refblock && numref != numTotalAtoms) {
              ostringstream os;
              os << "Number of atoms in REFPOSITION block in ";
              if (gin.posrest.nrdrx) os << s_coord;
              else os << s_refpos;
              os << " (" << numref << ")\n does not match total number of atoms ("
                      << numTotalAtoms << ")\n";
              printError(numWarnings, numErrors, os.str());
            }
          }

          if (!l_posresspec) {
            ifstream fin(filenames[posresspecfile].name(0).c_str());
            if (fin) {
              l_posresspec = 1;
              s_posresspec = filenames[posresspecfile].name(0);
              ostringstream os;
              os << "No posresspec-file specified, but I found"
                      << filenames[posresspecfile].name(0)
                      << " which I will use\n";
              printWarning(numWarnings, numErrors, os.str());
              fin.close();
            } else {
              ostringstream os;
              os << "NTR = " << gin.posrest.ntr << " in POSREST block, but no "
                      << "posresspec-file specified\n";
              printError(numWarnings, numErrors, os.str());
            }
          }
          if (l_posresspec) {
            Ginstream iprs(s_posresspec);
            fileInfo prs;
            iprs >> prs;
            iprs.close();
            int l_prs = 0;
            for (unsigned int i = 0; i < prs.blocks.size(); i++) {
              if ((!gromosXX) && prs.blocks[i] == "POSRESSPEC") l_prs = 1;
              if (gromosXX && prs.blocks[i] == "POSRES") l_prs = 1;
            }

            if (!l_prs) {
              string s = "No POSRESSPEC block in posresspec file (";
              if (gromosXX) s = "No POSRES block in posresspec file (";

              s += s_posresspec;
              s += ")\n";
              printError(numWarnings, numErrors, s);
            }
          }
        } else if (l_refpos || l_posresspec) {
          ostringstream os;
          if (l_refpos) os << "reference positions file ";
          if (l_refpos && l_posresspec) os << "and ";
          if (l_posresspec) os << "position restraints specification file ";
          os << "specified\n";
          os << "But no position res/constraining according to input file\n";
          printWarning(numWarnings, numErrors, os.str());
        }

        //PERTURB
        if ((gin.perturb.found && gin.perturb.ntg != 0) ||
            (gin.perturb03.found && gin.perturb03.ntg != 0)) {

          // check if perturbation topology is present
          if (!l_pttopo) {
            ifstream fin(filenames[pttopofile].name(0).c_str());
            if (fin) {
              l_pttopo = 1;
              s_pttopo = filenames[pttopofile].name(0);
              ostringstream os;
              os << "No perturbation topology specified, but I found "
                      << filenames[pttopofile].name(0)
                      << " which I will use\n";
              printWarning(numWarnings, numErrors, os.str());
              fin.close();
            } else {
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
          if (gin.perturb.found) {
            rlamfin = gin.perturb.rlam + gin.perturb.dlamt * gin.step.dt *
                    gin.step.nstlim;

            if (rlamfin > 1.0) {
              ostringstream os;
              os << "Using RLAM = " << gin.perturb.rlam << " and DLAMT = "
                      << gin.perturb.dlamt << " in the PERTURB block and NSTLIM = "
                      << gin.step.nstlim << " in the STEP block\n";
              os << "will lead to a final lambda value of " << rlamfin << endl;
              printWarning(numWarnings, numErrors, os.str());
            }
          } else {
            rlamfin = gin.perturb03.rlam + gin.perturb03.dlamt *
                    gin.step.dt * gin.step.nstlim;

            if (rlamfin > 1.0) {
              ostringstream os;
              os << "Using RLAM = " << gin.perturb.rlam << " and DLAMT = "
                      << gin.perturb.dlamt << " in the PERTURB block and NSTLIM = "
                      << gin.step.nstlim << " in the STEP block\n";
              os << "will lead to a final lambda value of " << rlamfin << endl;
              printWarning(numWarnings, numErrors, os.str());
            }
          }
        } else if (l_pttopo) {
          string s = "Perturbation topology specified, but no perturbation ";
          s += "according to the input file\n";
          printWarning(numWarnings, numErrors, s);
        }
      }

      // Now, write the stupid script
      if (numErrors != 0) {
        if (numErrors > 1) cout << "\n\nTHERE WERE " << numErrors << " ERRORS\n";
        else cout << "\n\nTHERE WAS 1 ERROR\n";

        if (args.count("force") == -1) {
          cout << "No script will be written\n";
          exit(1);
        } else {
          cout << "\nForcing script output (at your own risk)\n\n";
        }
      }
      if (numErrors == 0 && numWarnings == 0)
        cout << "OK" << endl << endl;


      //first check whether we should be in a different directory
      string subdir = simuldir;
      if (iter->second.dir != ".") {
        subdir += "/" + iter->second.dir;
      }
      mkdir(subdir.c_str(), 00755);
      chdir(subdir.c_str());
      cout << "Writing script: " << filenames[FILETYPE["script"]].name(0) << endl;

      ofstream fout(filenames[FILETYPE["script"]].name(0).c_str());
      fout.setf(ios::left, ios::adjustfield);
      fout << "#!/bin/sh" << endl;
      fout << "\n# first we set some variables\n";
      fout << "NAME=`whoami`\n";
      fout << "PROGRAM=" << args["bin"] << endl;
      fout << "SIMULDIR=" << simuldir << endl;
      fout << "\n# create temporary directory\n";
      fout << "WORKDIR=" << misc[0].name(0) << endl;
      fout << "mkdir -p ${WORKDIR}\n";
      fout << "cd       ${WORKDIR}\n";
      fout << "\n# set the input files\n";
      fout << "TOPO=${SIMULDIR}/" << s_topo << endl;
      fout << "IUNIT=${SIMULDIR}/";
      if (iter->second.dir != ".") fout << iter->second.dir << "/";
      fout << filenames[FILETYPE["input"]].name(0) << endl;

      if (iter != joblist.begin() || iter->second.dir != "." ||
          filenames[FILETYPE["input"]].name(0) != s_input) {
        // write the new input files
        cout << "     and input: " << filenames[FILETYPE["input"]].name(0);
        printInput(filenames[FILETYPE["input"]].name(0), gin);
      }

      cout << "\n--------------------------------------------------------------------------------" << endl;
      cout << endl;

      fout << "INPUTCRD=${SIMULDIR}/";
      if (iter == joblist.begin()) {
        fout << s_coord << endl;
      } else {
        if (iter->second.prev_id == -1) {
          //if(iter->second.dir!=".") fout << iter->second.dir << "/";
          fout << s_coord << endl;
        } else {
          jobinfo prevjob = joblist[iter->second.prev_id];
          if (prevjob.dir != ".") fout << prevjob.dir << "/";
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
          if (dual)
            filenames2[FILETYPE["coord"]].setInfo(systemname,
                  atof(iter->second.param["T"].c_str()),
                  atof(iter->second.param["DELTAT"].c_str()),
                  iter->first + dual_offset,
                  q);
        }
      }

      if (l_refpos) fout << "REFPOS=${SIMULDIR}/" << s_refpos << endl;
      if (l_posresspec) fout << "POSRESSPEC=${SIMULDIR}/"
              << s_posresspec << endl;
      if (l_disres) fout << "DISRES=${SIMULDIR}/" << s_disres << endl;
      if (l_dihres) fout << "DIHRES=${SIMULDIR}/" << s_dihres << endl;
      if (l_jvalue) fout << "JVALUE=${SIMULDIR}/" << s_jvalue << endl;
      if (l_ledih) fout << "LEDIH=${SIMULDIR}/" << s_ledih << endl;
      if (l_pttopo) fout << "PTTOPO=${SIMULDIR}/" << s_pttopo << endl;
      // any additional links?
      for (unsigned int k = 0; k < linkadditions.size(); k++)
        if (linkadditions[k] < 0)
          fout << linknames[k] << "=${SIMULDIR}/"
                << filenames[numFiletypes + k].name(0)
          << endl;

      fout << "\n#set the output files\n";
      fout << "OUNIT=" << filenames[FILETYPE["output"]].name(0) << endl;
      fout << "OUTPUTCRD=" << filenames[FILETYPE["coord"]].name(0) << endl;
      if (gin.write.ntwx || gin.writetraj.ntwx) fout << "OUTPUTTRX="
              << filenames[FILETYPE["outtrx"]].name(0)
        << endl;
      if (gin.write.ntwv || gin.writetraj.ntwv) fout << "OUTPUTTRV="
              << filenames[FILETYPE["outtrv"]].name(0)
        << endl;
      if (gin.write.ntwe || gin.writetraj.ntwe) fout << "OUTPUTTRE="
              << filenames[FILETYPE["outtre"]].name(0)
        << endl;
      if (gin.write.ntwg || gin.writetraj.ntwg) fout << "OUTPUTTRG="
              << filenames[FILETYPE["outtrg"]].name(0)
        << endl;
      if (gin.write.ntba > 0 || gin.writetraj.ntwb) fout << "OUTPUTBAE="
              << filenames[FILETYPE["outbae"]].name(0)
        << endl;

      if ((gin.write.ntba > 0 || gin.writetraj.ntwb) &&
          ((gin.perturb.found && gin.perturb.ntg > 0) ||
          (gin.perturb03.found && gin.perturb03.ntg > 0) ||
          (gin.perturbation.found && gin.perturbation.ntg > 0)))
        fout << "OUTPUTBAG="
              << filenames[FILETYPE["outbag"]].name(0)
        << endl;

      // any additional links?
      for (unsigned int k = 0; k < linkadditions.size(); k++)
        if (linkadditions[k] > 0)
          fout << linknames[k] << "="
                << filenames[numFiletypes + k].name(0) << endl;

      if (dual) {
        fout << "\n# output files of second job\n";
        fout << "OUNIT2=" << filenames2[FILETYPE["output"]].name(0) << endl;
        fout << "OUTPUTCRD2=" << filenames2[FILETYPE["coord"]].name(0) << endl;
        if (gin.write.ntwx || gin.writetraj.ntwx) fout << "OUTPUTTRX2="
                << filenames2[FILETYPE["outtrx"]].name(0)
          << endl;
        if (gin.write.ntwv || gin.writetraj.ntwv) fout << "OUTPUTTRV2="
                << filenames2[FILETYPE["outtrv"]].name(0)
          << endl;
        if (gin.write.ntwe || gin.writetraj.ntwe) fout << "OUTPUTTRE2="
                << filenames2[FILETYPE["outtre"]].name(0)
          << endl;
        if (gin.write.ntwg || gin.writetraj.ntwg) fout << "OUTPUTTRG2="
                << filenames2[FILETYPE["outtrg"]].name(0)
          << endl;
        if (gin.write.ntba > 0 || gin.writetraj.ntwb) fout << "OUTPUTBAE2="
                << filenames2[FILETYPE["outbae"]].name(0)
          << endl;

        if ((gin.write.ntba > 0 || gin.writetraj.ntwb) &&
            ((gin.perturb.found && gin.perturb.ntg > 0) ||
            (gin.perturb03.found && gin.perturb03.ntg > 0) ||
            (gin.perturbation.found && gin.perturbation.ntg > 0)))
          fout << "OUTPUTBAG2="
                << filenames2[FILETYPE["outbag"]].name(0)
          << endl;

        // any additional links?
        for (unsigned int k = 0; k < linkadditions.size(); k++)
          if (linkadditions2[k] > 0)
            fout << linknames2[k] << "="
                  << filenames2[numFiletypes + k].name(0) << endl;
      } // dual output


      if (misc[2].name(0) != "") {
        fout << "\n# first command\n";
        fout << misc[2].name(0) << "\n";
      }

      if (!gromosXX) {
        fout << "\n# link the files\n";
        fout << "rm -f fort.*\n";
        fout << setw(25) << "ln -s ${TOPO}" << " fort.20\n";
        fout << setw(25) << "ln -s ${INPUTCRD}" << " fort.21\n";
        if (l_refpos) fout << setw(25)
          << "ln -s ${REFPOS}" << " fort.22\n";
        if (l_posresspec) fout << setw(25)
          << "ln -s ${POSRESSPEC}" << " fort.23\n";
        if (l_disres) fout << setw(25)
          << "ln -s ${DISRES}" << " fort.24\n";
        if (l_dihres) fout << setw(25)
          << "ln -s ${DIHRES}" << " fort.25\n";
        if (l_jvalue) fout << setw(25)
          << "ln -s ${JVALUE}" << " fort.26\n";
        if (l_ledih) fout << setw(25)
          << "ln -s ${LEDIH}" << " fort.27\n";
        if (l_pttopo) fout << setw(25)
          << "ln -s ${PTTOPO}" << " fort.30\n";

        fout << setw(25) << "ln -s ${OUTPUTCRD}" << " fort.11\n";
        if (gin.write.ntwx || gin.writetraj.ntwx) fout << setw(25)
          << "ln -s ${OUTPUTTRX}" << " fort.12\n";
        if (gin.write.ntwv || gin.writetraj.ntwv) fout << setw(25)
          << "ln -s ${OUTPUTTRV}" << " fort.13\n";
        if (gin.write.ntwe || gin.writetraj.ntwe) fout << setw(25)
          << "ln -s ${OUTPUTTRE}" << " fort.15\n";
        if (gin.write.ntwg || gin.writetraj.ntwg) fout << setw(25)
          << "ln -s ${OUTPUTTRG}" << " fort.16\n";
        // any additional links
        for (unsigned int k = 0; k < linkadditions.size(); k++) {
          string s("ln -s ${" + linknames[k] + "}");
          fout << setw(25) << s << " fort." << abs(linkadditions[k])
                  << endl;
        }

        fout << "\n# run the program\n\n";
        fout << misc[3].name(0) << "${PROGRAM} < ${IUNIT} > ${OUNIT}\n";
      } else {

        fout << "\n\n";

        if (do_remd) {
          fout << "\n# run slave on single processor\n";
          fout << "OMP_NUM_THREADS=1\n\n";
        }
        fout << "MDOK=1\n\n";
        fout << misc[3].name(0) << "${PROGRAM}";

        fout << " \\\n\t" << setw(12) << "@topo" << " ${TOPO}";
        fout << " \\\n\t" << setw(12) << "@conf" << " ${INPUTCRD}";
        fout << " \\\n\t" << setw(12) << "@input" << " ${IUNIT}";
        if (l_pttopo) fout << " \\\n\t"
                << setw(12) << "@pttopo" << " ${PTTOPO}";
        if (l_posresspec) fout << " \\\n\t"
                << setw(12) << "@posres" << " ${POSRESSPEC}";
        if (l_disres) fout << " \\\n\t"
                << setw(12) << "@distrest" << " ${DISRES}";
        if (l_dihres) fout << " \\\n\t"
                << setw(12) << "@dihrest" << " ${DIHRES}";
        if (l_jvalue) fout << " \\\n\t"
                << setw(12) << "@jval" << " ${JVALUE}";

        fout << " \\\n\t" << setw(12) << "@fin" << " ${OUTPUTCRD}";
        if (gin.write.ntwx || gin.writetraj.ntwx) fout << " \\\n\t" << setw(12) << "@trj"
          << " ${OUTPUTTRX}";
        if (gin.write.ntwv || gin.writetraj.ntwv) fout << " \\\n\t" << setw(12) << "@trv"
          << " ${OUTPUTTRV}";
        if (gin.write.ntwe || gin.writetraj.ntwe) fout << " \\\n\t" << setw(12) << "@tre"
          << " ${OUTPUTTRE}";
        if (gin.write.ntwg || gin.writetraj.ntwg) fout << " \\\n\t" << setw(12) << "@trg"
          << " ${OUTPUTTRG}";
        if (gin.write.ntba > 0 || gin.writetraj.ntwb) fout << " \\\n\t" << setw(12) << "@bae"
          << " ${OUTPUTBAE}";

        if ((gin.write.ntba || gin.writetraj.ntwb) > 0 &&
            ((gin.perturb.found && gin.perturb.ntg > 0) ||
            (gin.perturb03.found && gin.perturb03.ntg > 0) ||
            (gin.perturbation.found && gin.perturbation.ntg > 0)))
          fout << " \\\n\t" << setw(12) << "@bag"
          << " ${OUTPUTBAG}";

        // any additional links
        for (unsigned int k = 0; k < linkadditions.size(); k++)
          fout << " \\\n\t@" << setw(11) << linknames[k]
                << " ${" << linknames[k] << "}";

        if (do_remd) {
          std::ostringstream os;
          os << "@slave " + hostname + " " << port;
          fout << "\\\n\t" << setw(25) << os.str();
        }

        fout << "\\\n\t" << setw(12) << ">" << " ${OUNIT}     || MDOK=0";

        if (dual) {
          fout << " &\n\n";

          if (do_remd) {
            fout << "\n# run slave on single processor\n";
            fout << "OMP_NUM_THREADS=1\n\n";
          }

          fout << misc2[3].name(0) << "${PROGRAM}";

          fout << " \\\n\t" << setw(12) << "@topo" << " ${TOPO}";
          fout << " \\\n\t" << setw(12) << "@conf" << " ${INPUTCRD}";
          fout << " \\\n\t" << setw(12) << "@input" << " ${IUNIT}";
          if (l_pttopo) fout << " \\\n\t"
                  << setw(12) << "@pttopo" << " ${PTTOPO}";
          if (l_posresspec) fout << " \\\n\t"
                  << setw(12) << "@posres" << " ${POSRESSPEC}";
          if (l_disres) fout << " \\\n\t"
                  << setw(12) << "@distrest" << " ${DISRES}";
          if (l_jvalue) fout << " \\\n\t"
                  << setw(12) << "@jval" << " ${JVALUE}";
          fout << " \\\n\t" << setw(12) << "@fin" << " ${OUTPUTCRD2}";
          if (gin.write.ntwx || gin.writetraj.ntwx) fout << " \\\n\t" << setw(12) << "@trj"
            << " ${OUTPUTTRX2}";
          if (gin.write.ntwv || gin.writetraj.ntwv) fout << " \\\n\t" << setw(12) << "@trv"
            << " ${OUTPUTTRV2}";
          if (gin.write.ntwe || gin.writetraj.ntwe) fout << " \\\n\t" << setw(12) << "@tre"
            << " ${OUTPUTTRE2}";
          if (gin.write.ntwg || gin.writetraj.ntwg) fout << " \\\n\t" << setw(12) << "@trg"
            << " ${OUTPUTTRG2}";
          if (gin.write.ntba > 0 || gin.writetraj.ntwb) fout << " \\\n\t" << setw(12) << "@bae"
            << " ${OUTPUTBAE2}";

          if ((gin.write.ntba || gin.writetraj.ntwb) > 0 &&
              ((gin.perturb.found && gin.perturb.ntg > 0) ||
              (gin.perturb03.found && gin.perturb03.ntg > 0) ||
              (gin.perturbation.found && gin.perturbation.ntg > 0)))
            fout << " \\\n\t" << setw(12) << "@bag"
            << " ${OUTPUTBAG2}";

          // any additional links
          for (unsigned int k = 0; k < linkadditions2.size(); k++)
            fout << " \\\n\t@" << setw(11) << linknames2[k]
                  << " ${" << linknames2[k] << "}";

          if (do_remd) {
            std::ostringstream os;
            os << "@slave " + hostname + " " << port;
            fout << "\\\n\t" << setw(25) << os.str();
          }

          fout << "\\\n\t" << setw(12) << ">" << " ${OUNIT2}";

          fout << "\n\nwait";
        }
        fout << "\n\n";
      }

      fout << "uname -a >> ${OUNIT}\n";
      if (dual) fout << "uname -a >> ${OUNIT2}\n";

      if (gin.write.ntwx || gin.write.ntwv || gin.write.ntwe || gin.write.ntwg ||
          gin.writetraj.ntwx || gin.writetraj.ntwv || gin.writetraj.ntwe || gin.writetraj.ntwg)
        fout << "\n# compress some files\n";
      if (gin.write.ntwx || gin.writetraj.ntwx) fout << "gzip ${OUTPUTTRX}\n";
      if (gin.write.ntwv || gin.writetraj.ntwv) fout << "gzip ${OUTPUTTRV}\n";
      if (gin.write.ntwe || gin.writetraj.ntwe) fout << "gzip ${OUTPUTTRE}\n";
      if (gin.write.ntwg || gin.writetraj.ntwg) fout << "gzip ${OUTPUTTRG}\n";
      if (gin.write.ntba > 0 || gin.writetraj.ntwb) fout << "gzip ${OUTPUTBAE}\n";
      if ((gin.write.ntba > 0 || gin.writetraj.ntwb) &&
          ((gin.perturb.found && gin.perturb.ntg > 0) ||
          (gin.perturb03.found && gin.perturb03.ntg > 0) ||
          (gin.perturbation.found && gin.perturbation.ntg > 0)))
        fout << "gzip ${OUTPUTBAG}\n";

      if (dual) {
        if (gin.write.ntwx || gin.writetraj.ntwx) fout << "gzip ${OUTPUTTRX2}\n";
        if (gin.write.ntwv || gin.writetraj.ntwv) fout << "gzip ${OUTPUTTRV2}\n";
        if (gin.write.ntwe || gin.writetraj.ntwe) fout << "gzip ${OUTPUTTRE2}\n";
        if (gin.write.ntwg || gin.writetraj.ntwg) fout << "gzip ${OUTPUTTRG2}\n";
        if (gin.write.ntba > 0 || gin.writetraj.ntwb) fout << "gzip ${OUTPUTBAE2}\n";
        if ((gin.write.ntba > 0 || gin.writetraj.ntwb) &&
            ((gin.perturb.found && gin.perturb.ntg > 0) ||
            (gin.perturb03.found && gin.perturb03.ntg > 0) ||
            (gin.perturbation.found && gin.perturbation.ntg > 0)))
          fout << "gzip ${OUTPUTBAG2}\n";
      }

      fout << "\n# copy the files back\n";
      fout << "OK=1\n";
      fout << setw(25) << "cp ${OUNIT}" << " ${SIMULDIR}";
      if (iter->second.dir != ".") fout << "/" << iter->second.dir;
      fout << " || OK=0\n";
      fout << setw(25) << "cp ${OUTPUTCRD}" << " ${SIMULDIR}";
      if (iter->second.dir != ".") fout << "/" << iter->second.dir;
      fout << " || OK=0\n";
      if (gin.write.ntwx || gin.writetraj.ntwx) {
        fout << setw(25) << "cp ${OUTPUTTRX}.gz" << " ${SIMULDIR}";
        if (iter->second.dir != ".") fout << "/" << iter->second.dir;
        fout << " || OK=0\n";
      }
      if (gin.write.ntwv || gin.writetraj.ntwv) {
        fout << setw(25) << "cp ${OUTPUTTRV}.gz" << " ${SIMULDIR}";
        if (iter->second.dir != ".") fout << "/" << iter->second.dir;
        fout << " || OK=0\n";
      }
      if (gin.write.ntwe || gin.writetraj.ntwe) {
        fout << setw(25) << "cp ${OUTPUTTRE}.gz" << " ${SIMULDIR}";
        if (iter->second.dir != ".") fout << "/" << iter->second.dir;
        fout << " || OK=0\n";
      }
      if (gin.write.ntwg || gin.writetraj.ntwg) {
        fout << setw(25) << "cp ${OUTPUTTRG}.gz" << " ${SIMULDIR}";
        if (iter->second.dir != ".") fout << "/" << iter->second.dir;
        fout << " || OK=0\n";
      }
      if (gin.write.ntba > 0 || gin.writetraj.ntwb) {
        fout << setw(25) << "cp ${OUTPUTBAE}.gz" << " ${SIMULDIR}";
        if (iter->second.dir != ".") fout << "/" << iter->second.dir;
        fout << " || OK=0\n";

        if ((gin.perturb.found && gin.perturb.ntg > 0) ||
            (gin.perturb03.found && gin.perturb03.ntg > 0) ||
            (gin.perturbation.found && gin.perturbation.ntg > 0)) {

          fout << setw(25) << "cp ${OUTPUTBAG}.gz" << " ${SIMULDIR}";
          if (iter->second.dir != ".") fout << "/" << iter->second.dir;
          fout << " || OK=0\n";
        }
      }

      // any additional links
      for (unsigned int k = 0; k < linkadditions.size(); k++) {
        if (linkadditions[k] > 0) {
          string s("cp ${" + linknames[k] + "}");
          string sgz("cp ${" + linknames[k] + "}.gz");
          fout << setw(25) << s << " ${SIMULDIR}";
          if (iter->second.dir != ".") fout << "/" << iter->second.dir;
          fout << " || OK=0\n";
        }
      }

      if (dual) {
        fout << setw(25) << "cp ${OUNIT2}" << " ${SIMULDIR}";
        if (iter->second.dir != ".") fout << "/" << iter->second.dir;
        fout << " || OK=0\n";
        fout << setw(25) << "cp ${OUTPUTCRD2}" << " ${SIMULDIR}";

        if (iter->second.dir != ".") fout << "/" << iter->second.dir;
        fout << " || OK=0\n";
        if (gin.write.ntwx || gin.writetraj.ntwx) {
          fout << setw(25) << "cp ${OUTPUTTRX2}.gz" << " ${SIMULDIR}";
          if (iter->second.dir != ".") fout << "/" << iter->second.dir;
          fout << " || OK=0\n";
        }
        if (gin.write.ntwv || gin.writetraj.ntwv) {
          fout << setw(25) << "cp ${OUTPUTTRV2}.gz" << " ${SIMULDIR}";
          if (iter->second.dir != ".") fout << "/" << iter->second.dir;
          fout << " || OK=0\n";
        }
        if (gin.write.ntwe || gin.writetraj.ntwe) {
          fout << setw(25) << "cp ${OUTPUTTRE2}.gz" << " ${SIMULDIR}";
          if (iter->second.dir != ".") fout << "/" << iter->second.dir;
          fout << " || OK=0\n";
        }
        if (gin.write.ntwg || gin.writetraj.ntwg) {
          fout << setw(25) << "cp ${OUTPUTTRG2}.gz" << " ${SIMULDIR}";
          if (iter->second.dir != ".") fout << "/" << iter->second.dir;
          fout << " || OK=0\n";
        }
        if (gin.write.ntba > 0 || gin.writetraj.ntwb) {
          fout << setw(25) << "cp ${OUTPUTBAE2}.gz" << " ${SIMULDIR}";
          if (iter->second.dir != ".") fout << "/" << iter->second.dir;
          fout << " || OK=0\n";

          if ((gin.perturb.found && gin.perturb.ntg > 0) ||
              (gin.perturb03.found && gin.perturb03.ntg > 0) ||
              (gin.perturbation.found && gin.perturbation.ntg > 0)) {

            fout << setw(25) << "cp ${OUTPUTBAG2}.gz" << " ${SIMULDIR}";
            if (iter->second.dir != ".") fout << "/" << iter->second.dir;
            fout << " || OK=0\n";
          }
        }

        // any additional links
        for (unsigned int k = 0; k < linkadditions2.size(); k++) {
          if (linkadditions2[k] > 0) {
            string s("cp ${" + linknames2[k] + "}");
            fout << setw(25) << s << " ${SIMULDIR}";

            if (iter->second.dir != ".") fout << "/" << iter->second.dir;
            fout << " || OK=0\n";
          }
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

      fout << "\n# stop if MD was not succesfull\n";
      fout << "if `test ${MDOK} -eq 0`; then\n";
      fout << "  exit\n";
      fout << "fi\n";

      fout << "\n# perform last command (usually submit next job)\n";
      // which job do we have to submit (also check in the earlier ones?)
      map<int, jobinfo>::const_iterator it = joblist.begin();
      while (it != to) {
        if (it->second.prev_id == iter->first) {
          setParam(gin, it->second);
          misc[1].setInfo(systemname,
                  atof(iter->second.param["ENDTIME"].c_str()),
                  gin.step.dt * gin.step.nstlim, it->first, q);
          if (it->first != iter->first) {

            fout << "cd ${SIMULDIR}";
            if (it->second.dir != ".") fout << "/" << it->second.dir;
            fout << "\n";
            fout << misc[1].name(0) << endl;
          }

        }
        ++it;
      }

      fout.close();
      chmod(filenames[FILETYPE["script"]].name(0).c_str(), 00755);
    }
  } catch (const gromos::Exception &e) {
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
    fout << "\tAutomatically generated input file\n\t";
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
    throw gromos::Exception("mk_script", "Reading of jobscript file failed. "
			    "No JOBSCRIPTS block");
  istringstream iss(buffer[1]);
  vector<string> head;
  string b;
  while((iss>>b) !=0) head.push_back(b);
  if(head[0]!="job_id" || head.back() !="run_after" 
     || head[head.size()-2]!="subdir")
    throw gromos::Exception("mk_script", "Reading of jobscript file failed.\n"
			    "First line syntax:\n"
			    "job_id PARAM PARAM ... subdir run_after");
  if(buffer.back().find("END")!=0)
	throw gromos::Exception("mk_script", "Jobscript file " 
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
	throw gromos::Exception("mk_script", "Template file " 
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
	  case outbaefile:     names[outbaefile].setTemplate(temp);     break;
	  case outbagfile:     names[outbagfile].setTemplate(temp);     break;
	  case scriptfile:     names[scriptfile].setTemplate(temp);     break;
	  case unknownfile:
	    printWarning(w,e, "Don't know how to handle template for "+sdum
			 +". Ingoring");
	}
      }
    }
    if(buffer.size() && first=="MISCELLANEOUS"){
      if(buffer[buffer.size()-1].find("END")!=0)
	throw gromos::Exception("mk_script", "Template file " + 
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
	if(sdum=="firstcommand") {
	  ostringstream os;
	  while(!iss.eof()){
	    iss >> sdum;
	    os << sdum << " ";
	  }
	  misc[2].setTemplate(os.str());
	}
	if(sdum=="mpicommand") {
	  ostringstream os;
	  while(!iss.eof()){
	    iss >> sdum;
	    os << sdum << " ";
	  }
	  misc[3].setTemplate(os.str());
	}
      }
      // re-set the standard lastcommand template, in case the script template
      // has changed
      if(!l_lastcommand)
	misc[1].setTemplate(submitcommand+names[FILETYPE["script"]].temp());
    
    }
    if(buffer.size() && first=="LINKADDITION"){
      if(buffer[buffer.size()-1].find("END")!=0)
	throw gromos::Exception("mk_script", "Template file " + 
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
    // SYSTEM
    if(iter->first=="NPM")
      gin.system.npm=atoi(iter->second.c_str());
    else if(iter->first=="NSM")
      gin.system.nsm=atoi(iter->second.c_str());
    // MINIMISE (g96) and ENERGYMIN (md++ and promd)
    else if(iter->first=="NTEM"){
      if(gin.minimise.found)
        gin.minimise.ntem=atoi(iter->second.c_str());
      else
        gin.energymin.ntem=atoi(iter->second.c_str());
    }
    else if(iter->first=="NCYC"){
      if(gin.minimise.found)
        gin.minimise.ncyc=atoi(iter->second.c_str());
      else
        gin.energymin.ncyc=atoi(iter->second.c_str());
    }
    else if(iter->first=="DELE"){
      if(gin.minimise.found)
        gin.minimise.dele=atof(iter->second.c_str());
      else
        gin.energymin.dele=atof(iter->second.c_str());
    }
    else if(iter->first=="DX0"){
      if(gin.minimise.found)
        gin.minimise.dx0=atof(iter->second.c_str());
      else
        gin.energymin.dx0=atof(iter->second.c_str());
    }
    else if(iter->first=="DXM"){
      if(gin.minimise.found)
        gin.minimise.dxm=atof(iter->second.c_str());
      else
        gin.energymin.dxm=atof(iter->second.c_str());
    }
    else if(iter->first=="NMIN")
      gin.energymin.nmin=atoi(iter->second.c_str());
    else if(iter->first=="FLIM")
      gin.energymin.flim=atof(iter->second.c_str());
    // STOCHASTIC (g96) and STOCHDYN (promd, md++)
    // START (g96) and INITIALISE (md++ and promd)
    else if(iter->first=="NTX")
      gin.start.ntx=atoi(iter->second.c_str());
    else if(iter->first=="INIT")
      gin.start.init=atoi(iter->second.c_str());
    else if(iter->first=="IG"){
      if(gin.start.found)
        gin.start.ig=atoi(iter->second.c_str());
      else
        gin.initialise.ig=atoi(iter->second.c_str());
    }
    else if(iter->first=="TEMPI"){
      if(gin.start.found)
        gin.start.tempi=atof(iter->second.c_str());
      else
        gin.initialise.tempi=atof(iter->second.c_str());
    }
    else if(iter->first=="HEAT")
      gin.start.heat=atof(iter->second.c_str());
    else if(iter->first=="NTXO")
      gin.start.ntx0=atoi(iter->second.c_str());
    else if(iter->first=="BOLTZ")
      gin.start.boltz=atof(iter->second.c_str());
    // READTRAJ (promd, md++)
    // CONSISTENCYCHECK (promd)
    // STEP (all)
    else if(iter->first=="NSTLIM")
      gin.step.nstlim=atoi(iter->second.c_str());
    else if(iter->first=="T")
      gin.step.t=atof(iter->second.c_str());
    else if(iter->first=="DT")
      gin.step.dt=atof(iter->second.c_str());
    // BOUNDARY (g96) and BOUNDCOND (promd, md++)
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
    // MULTICELL (promd, md++)
    // SUBMOLECULES (g96) do not implement here!
    // TCOUPLE (g96)
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
    // PCOUPLE (g96)
    else if(iter->first=="NTP")
      gin.pcouple.ntp=atoi(iter->second.c_str());
    else if(iter->first=="PRES0")
      gin.pcouple.pres0=atof(iter->second.c_str());
    else if(iter->first=="COMP")
      gin.pcouple.comp=atof(iter->second.c_str());
    else if(iter->first=="TAUP")
      gin.pcouple.taup=atof(iter->second.c_str());
    // THERMOSTAT (promd)
    // MULTIBATH (md++)
    // BAROSTAT (promd)
    // VIRIAL (promd)
    // PRESSURESCALE (md++)
    // CENTREOFMASS (g96)
    else if(iter->first=="NDFMIN")
      gin.centreofmass.ndfmin=atoi(iter->second.c_str());
    else if(iter->first=="NTCM")
      gin.centreofmass.ntcm=atoi(iter->second.c_str());
    else if(iter->first=="NSCM")
      gin.centreofmass.nscm=atoi(iter->second.c_str());
    // OVERALLTRANSROT (promd)
    // COMTRANSROT (md++)
    // PRINT (g96)
    else if(iter->first=="NTPR")
      gin.print.ntpr=atoi(iter->second.c_str());
    else if(iter->first=="NTPL")
      gin.print.ntpl=atoi(iter->second.c_str());
    else if(iter->first=="NTPP")
      gin.print.ntpp=atoi(iter->second.c_str());
    // PRINTOUT (promd, md++)
    // WRITE (g96)
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
    // WRITETRAJ (promd, md++)
    // EWARN (md++)
    // DEBUG (promd)
    // SHAKE (g96)
    else if(iter->first=="NTC")
      gin.shake.ntc=atoi(iter->second.c_str());
    else if(iter->first=="TOL")
      gin.shake.tol=atof(iter->second.c_str());
    // GEOMCONSTRAINTS (promd)
    // CONSTRAINT (md++)
    // FORCE (all)
    // COVALENTFORM (promd, md++)
    // PLIST (g96) and PLIST03 (old md++ format)
    else if(iter->first=="NTNB")
      gin.plist.ntnb=atoi(iter->second.c_str());
    else if(iter->first=="NSNB")
      gin.plist.nsnb=atoi(iter->second.c_str());
    else if(iter->first=="RCUTP")
      gin.plist.rcutp=atof(iter->second.c_str());
    else if(iter->first=="RCUTL")
      gin.plist.rcutl=atof(iter->second.c_str());
    // NEIGHBOURLIST (promd)
    // PAIRLIST (md++)
    // NONBONDED (promd)
    // LONGRANGE (md++, g96)
    else if(iter->first=="EPSRF")
      gin.longrange.epsrf=atof(iter->second.c_str());
    else if(iter->first=="APPAK")
      gin.longrange.appak=atof(iter->second.c_str());
    else if(iter->first=="RCRF")
      gin.longrange.rcrf=atoi(iter->second.c_str());
    // CGRAIN (md++)
    // POSREST (g96)
    else if(iter->first=="NTR")
      gin.posrest.ntr=atoi(iter->second.c_str());
    else if(iter->first=="CHO")
      gin.posrest.cho=atof(iter->second.c_str());
    else if(iter->first=="NRDRX")
      gin.posrest.nrdrx=atoi(iter->second.c_str());
    // POSITIONRES (promd, md++)
    // DISTREST (g96)
    else if(iter->first=="NTDR")
      gin.distrest.ntdr=atoi(iter->second.c_str());
    else if(iter->first=="CDIS")
      gin.distrest.cdis=atoi(iter->second.c_str());
    else if(iter->first=="DR0")
      gin.distrest.dr0=atoi(iter->second.c_str());
    else if(iter->first=="TAUDR")
      gin.distrest.taudr=atoi(iter->second.c_str());
    else if(iter->first=="NRDDR")
      gin.distrest.nrddr=atoi(iter->second.c_str());
    // DISTANCERES (promd, md++)
    // DIHEREST (g96) and DIHEDRALRES (promd, md++)
    // J-VALUE (g96)
    // JVALUERES (promd, md++)
    // LOCALELEVATION (g96)
    // LOCALELEV (promd, md++)
    // ROTTRANS (old md++)
    // PERTURB (g96), PERTURBATION (promd, md++) and PERTURB03 (old md++)
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
    else if(iter->first=="ALPHLJ"){
      if (gin.perturb.found)
	gin.perturb.alphlj=atof(iter->second.c_str());
      else
	gin.perturb03.alphlj=atof(iter->second.c_str());
    }
    else if(iter->first=="ALPHC"){
      if (gin.perturb.found)
	gin.perturb.alphc=atof(iter->second.c_str());
      else
	gin.perturb03.alphc=atof(iter->second.c_str());
    }
    else if(iter->first=="NLAM"){
      if (gin.perturb.found)
	gin.perturb.nlam=atoi(iter->second.c_str());
      else
	gin.perturb03.nlam=atoi(iter->second.c_str());
    }
    else if(iter->first=="MMU")
      gin.perturb.mmu=atoi(iter->second.c_str());
    // LAMBDAS (md++)
    // UMBRELLA (promd)
    // PERSCALE (md++)
    // REPLICA (md++)
    // INNERLOOP (md++)
    // FOURDIM (g96)
    // GROMOS96COMPAT (promd)
    // INTEGRATE (md++)
    
    else if(iter->first=="ENDTIME" || iter->first=="DELTAT"){}
    else
      throw gromos::Exception("mk_script", "Cannot automatically change "
			      +iter->first +" in input file");
  }
}

