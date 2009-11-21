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
 * A molecular dynamics simulation is usually performed by executing a
 * small script that combines all the necessary files and redirects the
 * output to the appropriate places. When simulations are
 * performed on a queue, such scripts become indispensable. Additionally, in
 * many simulation projects the user prepares similar input files and
 * scripts over and over again.
 * Program mk_script can either create a series of similar scripts that
 * run a sequential list of simulations (keyword \@script) or it can
 * create scripts for a more complex set of simulations that perform a
 * specific task (start-up, perturbation; keyword \@joblist).
 * Scripts for special cases such as REMD simulations (keyword \@remd)
 * can also be written.
 * 
 * GROMOS does not require specific filenames for specific types of files.
 * However, most users find it useful to retain some order in
 * their filenames. mk_script has a standard way of constructing
 * filenames that depends on the script number and the system name. The
 * user can specify another set of rules to create filenames through
 * the mk_script library file (keyword \@template). In this file,
 * machine-dependent modifications to the scripts that are to be written can also be
 * specified, such as the job submission command and which queue to use (keyword
 * \@queue). A standard location of the mk_script library file can be
 * specified through the environment variable MK_SCRIPT_TEMPLATE.
 *
 * Program mk_script can write input files for promd or md++ (keyword \@version).
 * The md input file (keyword \@files->input) should also be of the correct
 * format: mk_script cannot convert program-specific md input blocks into the
 * analogous blocks for the other version of GROMOS.
 *
 * In addition to writing out scripts and input files, mk_script performs a number
 * of tests on the given input files to prevent the user from submitting
 * a simulation that will fail within the first few seconds. In the messages produced by
 * these tests, a distinction is made between warnings and errors. A
 * warning is given for inconsistencies in the inputs that may lead to an erroneous simulation,
 * but could also be intentional. An error is produced by inconsistencies that will definitely
 * result in the program crashing. Note that passing the tests carried out by mk_script
 * does not guarantee that a simulation will work, as these checks are not exhaustive.
 * Some of the more common warnings and errors are listed below.
 * If there are no errors, all requested input files and scripts will
 * be written to disc. If there are errors, the scripts will not be
 * written, unless the user forces this (keyword \@force).
 *
 *
 * <b>Warnings:</b>
 * <ol>
 * <li>  the GROMOS binary specified in the mk_script input file cannot be found
 * <li>  the highest LAST atom number in the MULTIBATH block in the md input file is not
 *       equal to the total number of atoms calculated from the topology file and SYSTEM block
 * <li>  DT in the STEP block is larger than the suggested value according to
 *       the value of NTC in the CONSTRAINT block in the md input file (no SHAKE, NTC=1:
 *       DT=0.0005; SHAKE bonds with H, NTC=2: DT=0.001; SHAKE all bonds, NTC=3: DT=0.002)
 * <li>  the FORCE for bonds that are SHAKEn is calculated
 * <li>  the FORCE for bonds that are not SHAKEn is not calculated
 * <li>  the shortest edge-to-edge distance of the periodic box is less than
 *       twice the long-range cut-off RCUTL in the PAIRLIST block of the md input file
 * <li>  the reaction field cut-off distance RCRF in the NONBONDED block of the md input file
 *       is not equal to the long-range cut-off RCUTL in the PAIRLIST block
 * <li>  a perturbation topology was specified in the mk_script input file but no
 *       perturbation was requested in the md input file
 * <li>  the combination of RLAM and DLAMT in the PERTURBATION block and the number of
 *       steps from the STEP block in the md input file will lead to a lambda value larger than 1
 * </ol>
 *
 *
 * <b>Errors:</b>
 * <ol>
 * <li>  there is no SYSTEM block in the md input file
 * <li>  there is no VELOCITY block in the coordinate file, but NTIVEL in the
 *       INITIALISE block of the md input file specifies that the velocities should be read from file
 * <li>  non-zero NTISHI in the INITIALISE block of the md input file specifies that the lattice shifts
 *       should be initialised, but zero NTB in the BOUNDCOND block specifies a vacuum simulation
 * <li>  there is no LATTICESHIFT block in the coordinate file, but NTISHI in the
 *       INITIALISE block of the md input file specifies that the lattce shifts should be read from file
 * <li>  there is no GENBOX block in the coordinate file, but non-zero NTB in the BOUNDCOND block
 *       specifies a non-vacuum simulation
 * <li>  there is no STEP block in the md input file
 * <li>  there is no BOUNDCOND block in the md input file
 * <li>  the number of the last atom given in the FORCE block of the md input file is not equal to the total
 *       number of atoms calculated from the topology and SYSTEM block
 * <li>  the number of atoms calculated from the topology and SYSTEM block of the md input file is not equal
 *       to the number of atoms in the POSITION block of the coordinate file
 * <li>  there is no FORCE block in the md input file
 * <li>  there is no PAIRLIST block in the md input file
 * <li>  in the PAIRLIST block, the short-range cutoff RCUTP is larger than the
 *       long-range cutoff RCUTL
 * <li>  there is no NONBONDED block in the md input file
 * <li>  no position restraints specification file is specified in the mk_script input file,
 *       but position restraining is switched on in the md input file
 * <li>  the position restraints specification file does not contain a
 *       POSRESSPEC (promd) or POSRES (md++) block
 * <li>  no perturbation topology file is specified in the mk_script input file,
 *       but perturbation is switched on in the md input file
 * </ol>
 *
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@sys</td><td>&lt;system name&gt; </td></tr>
 * <tr><td> \@bin</td><td>&lt;GROMOS binary to use&gt; </td></tr>
 * <tr><td> \@version</td><td>&lt;md++ / promd&gt; </td></tr>
 * <tr><td> \@dir</td><td>&lt;where the files should be (all filenames to be given relative to this)&gt; </td></tr>
 * <tr><td> [\@script</td><td>&lt;first script&gt; &lt;number of scripts&gt;] </td></tr>
 * <tr><td> [\@joblist</td><td>&lt;joblist file&gt;] </td></tr>
 * <tr><td> \@files</td><td></td></tr>
 * <tr><td> [\@template</td><td>&lt;mk_script library file&gt;] </td></tr>
 * <tr><td> [\@queue</td><td>&lt;which queue to use in mk_script library file&gt;] </td></tr>
 * <tr><td> [\@remd</td><td>&lt;master / slave hostname port&gt; (replica exchange MD)] </td></tr>
 * <tr><td> [\@cmd</td><td>&lt;overwrite last command in mk_script library file&gt;] </td></tr>
 * <tr><td> [\@force</td><td>(write script regardless of errors)] </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  mk_script
    @sys         ex
    @bin         /usr/local/gromos/md++/bin/md
    @version     md++
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
#include "../src/gcore/Box.h"
#include "../src/gio/InTopology.h"
#include "../src/gio/Ginstream.h"

using namespace std;
using namespace gcore;
using namespace gio;
using namespace args;

void printIO(string b, string var, string val, string allow);

void printWarning(string s);
void printError(string s);

#include "mk_script.h"

void printInput(string ofile, input gin);
void readLibrary(string file, vector<filename> &names,
        vector<filename> &misc,
        vector<string> &linknames, vector<int> &linkadditions,
        string system, string queue, double t,
        double dt, int ns);
void readJobinfo(string file, map<int, jobinfo> &ji);
void setParam(input &gin, jobinfo const &job);

int main(int argc, char **argv) {

  Argument_List knowns;
  knowns << "sys" << "script" << "bin" << "dir" << "queue" << "remd"
          << "files" << "template" << "version" << "cmd" << "joblist"
          << "force";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@sys  <system name>\n";
  usage += "\t@bin           <gromos96 binary to use>\n";
  usage += "\t@dir           <where should the files be>\n";
  usage += "\t@version       <md++ / promd>\n";
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
  usage += "\t\t[xray        <xray restraints file>]\n";
  usage += "\t[@template     <template filename, absolute or relative to @dir>]\n";
  usage += "\t[@queue        <queue flags>]\n";
  usage += "\t[@remd         <master / slave hostname port> (replica exchange MD)]\n";
  usage += "\t[@cmd          <overwrite last command>]\n";
  usage += "\t[@force        (write script regardless of errors)]\n";

  try {

    Arguments args(argc, argv, knowns, usage);

    // This restriction is maybe a bit overdone, because we could in principle
    // of course still read the old topology etc. But the error checking 
    // becomes a burden. I would suggest to keep a separategromos96 version in
    // contrib if necessary
    if (args::Arguments::inG96 == true || args::Arguments::outG96 == true) {
      throw gromos::Exception("mk_script",
              "This program no longer supports the gromos96 formats");
    }

    // first get some input parameters
    int scriptNumber = 1, numScripts = 1;
    string simuldir;
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
                "Specified directory does not exist");
      } else
        simuldir = "`pwd`";

      iter = args.lower_bound("script");
      if (iter != args.upper_bound("script")) {
        scriptNumber = atoi(iter->second.c_str());
        ++iter;
      }
      if (iter != args.upper_bound("script"))
        numScripts = atoi(iter->second.c_str());

      if (numScripts < 0)
        throw Arguments::Exception("Can't deal with negative number of scripts in @script argument");
    }
    string systemname = args["sys"];

    // read in the library
    string libraryfile;
    if (args.lower_bound("template") == args.upper_bound("template")) {
      // try to get it from environment
      if (getenv("MK_SCRIPT_TEMPLATE")) {
        libraryfile = getenv("MK_SCRIPT_TEMPLATE");
      } else {
        throw gromos::Exception("mk_script", "Please give @template or set the "
                "MK_SCRIPT_TEMPLATE environment variable.");
      }
    } else { // supplied by @template
      libraryfile = args["template"];
    }

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

    // parse the files
    int l_coord = 0, l_topo = 0, l_input = 0, l_refpos = 0, l_posresspec = 0, l_xray = 0;
    int l_disres = 0, l_dihres = 0, l_jvalue = 0, l_ledih = 0, l_pttopo = 0;
    string s_coord, s_topo, s_input, s_refpos, s_posresspec, s_xray;
    string s_disres, s_dihres, s_jvalue, s_ledih, s_pttopo;
    for (Arguments::const_iterator iter = args.lower_bound("files"),
            to = args.upper_bound("files"); iter != to; ++iter) {
      switch (FILETYPE[iter->second]) {
        case anatrxfile: ++iter;
        break;
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
        case xrayfile:++iter;
          s_xray = iter->second;
          l_xray = 1;
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
          printWarning(iter->second + " not used");
          break;
        case outtrxfile: ++iter;
          printWarning(iter->second + " not used");
          break;
        case outtrvfile: ++iter;
          printWarning(iter->second + " not used");
          break;
        case outtrffile: ++iter;
          printWarning(iter->second + " not used");
        case outtrefile: ++iter;
          printWarning(iter->second + " not used");
          break;
        case outtrgfile: ++iter;
          printWarning(iter->second + " not used");
          break;
        case outbaefile: ++iter;
          printWarning(iter->second + " not used");
          break;
        case outbagfile: ++iter;
          printWarning(iter->second + " not used");
          break;
        case outtrsfile: ++iter;
          printWarning(iter->second + " not used");
          break;
        case scriptfile: ++iter;
          printWarning(iter->second + " not used");
          break;
        case unknownfile: printError("Don't know how to handle file "
                  + iter->second);
      }
    }

    // check which outformat we want (gromos96 or gromosXX)
    // in g08 only relevant for job scripts!
    bool gromosXX = false;
    if (args["version"] == "md++")
      gromosXX = true;
    else if (args["version"] == "promd")
      gromosXX = false;
    else
      throw gromos::Exception("mk_script", "Please specify version (md++ or promd).");


    // read topology
    if (!l_topo) {
      throw gromos::Exception("mk_script", "You have to specify a topology\n" + usage);
    }
    InTopology it(s_topo);
    System sys = it.system();
    fileInfo topo;
    {
      Ginstream itopo(s_topo);
      itopo >> topo;
      itopo.close();
    }


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

    // replace the %queue% variable?
    string queue = "";
    {
      ostringstream os;
      Arguments::const_iterator iter = args.lower_bound("queue"),
              to = args.upper_bound("queue");
      for (; iter != to; ++iter) {
        std::string s = iter->second;
        if (s.find("\\n") != string::npos)
          s.replace(s.find("\\n"), 2, "\n");
        else s += " ";

        // os << iter->second << " ";
        os << s;
      }
      queue = os.str();
    }

    // create names for automated file names
    vector<filename> filenames;
    vector<filename> misc;
    vector<int> linkadditions;
    vector<string> linknames;

    for (int i = 0; i < numFiletypes; i++) {
      filename newname(systemname, gin.step.t, gin.step.nstlim * gin.step.dt,
              scriptNumber, queue);
      filenames.push_back(newname);
    }
    // workdir lastcommand firstcommand mpicommand stopcommand
    for (int i = 0; i < 5; i++) {
      filename newname(systemname, gin.step.t, gin.step.nstlim * gin.step.dt,
              scriptNumber, queue);
      misc.push_back(newname);
    }

    // set the standard templates
    filenames[FILETYPE["script"]].setTemplate("%system%_%number%.run");
    filenames[FILETYPE["input"]].setTemplate("%system%_%number%.imd");
    filenames[FILETYPE["topo"]].setTemplate("%system%.top");
    filenames[FILETYPE["refpos"]].setTemplate("%system%_%number%.rpr");
    filenames[FILETYPE["posresspec"]].setTemplate("%system%_%number%.por");
    filenames[FILETYPE["anatrj"]].setTemplate("%system%_%number%.trc.gz");
    filenames[FILETYPE["xray"]].setTemplate("%system%_%number%.xrs");
    filenames[FILETYPE["disres"]].setTemplate("%system%_%number%.dsr");
    filenames[FILETYPE["pttopo"]].setTemplate("%system%.ptp");
    filenames[FILETYPE["dihres"]].setTemplate("%system%_%number%.dhr");
    filenames[FILETYPE["jvalue"]].setTemplate("%system%_%number%.jvr");
    filenames[FILETYPE["ledih"]].setTemplate("%system%_%number%.led");
    filenames[FILETYPE["coord"]].setTemplate("%system%_%number%.cnf");
    filenames[FILETYPE["output"]].setTemplate("%system%_%number%.omd");
    filenames[FILETYPE["outtrx"]].setTemplate("%system%_%number%.trc");
    filenames[FILETYPE["outtrv"]].setTemplate("%system%_%number%.trv");
    filenames[FILETYPE["outtrf"]].setTemplate("%system%_%number%.trf");
    filenames[FILETYPE["outtre"]].setTemplate("%system%_%number%.tre");
    filenames[FILETYPE["outtrg"]].setTemplate("%system%_%number%.trg");
    filenames[FILETYPE["outbae"]].setTemplate("%system%_%number%.bae");
    filenames[FILETYPE["outbag"]].setTemplate("%system%_%number%.bag");
    filenames[FILETYPE["outtrs"]].setTemplate("%system%_%number%.trs");

    // And here is a gromos-like function call!
    readLibrary(libraryfile, filenames, misc,
            linknames, linkadditions,
            systemname, queue, gin.step.t,
            gin.step.nstlim * gin.step.dt,
            scriptNumber);

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
        printWarning(os.str());
        s_coord = filenames[FILETYPE["coord"]].name(-1);
        l_coord = 1;
      } else {
        ostringstream os;
        os << "No coordinate file is specified, some checks are not performed\n";
        os << "Assuming it does not exist yet, I will use "
                << filenames[FILETYPE["coord"]].name(-1)
                << " in the script\n";
        s_coord = filenames[FILETYPE["coord"]].name(-1);

        printWarning(os.str());
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
        printWarning("Specified binary not found! "
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
                iter->first, queue);
      }
      for (unsigned int i = 0; i < misc.size(); i++) {
        misc[i].setInfo(systemname, gin.step.t, gin.step.dt * gin.step.nstlim,
                iter->first, queue);
      }

      // define a variable to check if we are in the first script or not
      bool first_script = iter == joblist.begin();

      // do some stuff for first script or some (more than the standard 3) para-
      // meters have changed (due to the jobfile)
      if (first_script || iter->second.param.size() != 3) {
        cout << "Performing checks for script " << iter->first << endl;
        cout << "--------------------------------------------" << endl;
        cout << endl;

        // Ignore md++ specific blocks if promd input is to be written
        // and vice versa (start)
        if (!gromosXX) { // Ignore md++ specific blocks
          if (gin.multibath.found) {
            if (gin.thermostat.found) {
              printWarning("Ignored md++ specific block MULTIBATH\n");
              gin.multibath.found = 0;
            } else {
              printError("You want to perform a PROMD run and you have specified the md++ specific\n"
                      "block MULTIBATH. Maybe you want to specify the THERMOSTAT block instead?\n");
            }
          }
          if (gin.pressurescale.found) {
            if (gin.barostat.found) {
              printWarning("Ignored md++ specific block PRESSURESCALE\n");
              gin.pressurescale.found = 0;
            } else {
              printError("You want to perform a PROMD run and you have specified the md++ specific\n"
                      "block PRESSURESCALE. Maybe you want to specify the BAROSTAT block instead?\n");
            }
          }
          if (gin.comtransrot.found) {
            if (gin.overalltransrot.found) {
              printWarning("Ignored md++ specific block COMTRANSROT\n");
              gin.comtransrot.found = 0;
            } else {
              printError("You want to perform a PROMD run and you have specified the md++ specific\n"
                      "block COMTRANSROT. Maybe you want to specify the OVERALLTRANSROT block instead?\n");
            }
          }
          if (gin.ewarn.found) {
            printWarning("Ignored md++ specific block EWARN\n");
            gin.ewarn.found = 0;
          }
          if (gin.constraint.found) {
            if (gin.geomconstraints.found) {
              printWarning("Ignored md++ specific block CONSTRAINT\n");
              gin.constraint.found = 0;
            } else {
              printError("You want to perform a PROMD run and you have specified the md++ specific\n"
                      "block CONSTRAINT. Maybe you want to specify the GEOMCONSTRAINTS block instead?\n");
            }
          }
          if (gin.pairlist.found) {
            if (gin.neighbourlist.found) {
              printWarning("Ignored md++ specific block PAIRLIST\n");
              gin.pairlist.found = 0;
            } else {
              printError("You want to perform a PROMD run and you have specified the md++ specific\n"
                      "block PAIRLIST. Maybe you want to specify the NEIGHBOURLIST block instead?\n");
            }
          }
          if (gin.cgrain.found) {
            printWarning("Ignored md++ specific block CGRAIN\n");
            gin.cgrain.found = 0;
          }
          if (gin.rottrans.found) {
            printWarning("Ignored md++ specific block ROTTRANS\n");
            gin.rottrans.found = 0;
          }
          if (gin.lambdas.found) {
            printWarning("Ignored md++ specific block LAMBDAS\n");
            gin.lambdas.found = 0;
          }
          if (gin.perscale.found) {
            printWarning("Ignored md++ specific block PERSCALE\n");
            gin.perscale.found = 0;
          }
          if (gin.polarise.found) {
            printWarning("Ignored md++ specific block POLARISE\n");
            gin.polarise.found = 0;
          }
          if (gin.replica.found) {
            printWarning("Ignored md++ specific block REPLICA\n");
            gin.replica.found = 0;
          }
          if (gin.innerloop.found) {
            printWarning("Ignored md++ specific block INNERLOOP\n");
            gin.innerloop.found = 0;
          }
          if (gin.integrate.found) {
            printWarning("Ignored md++ specific block INTEGRATE\n");
            gin.integrate.found = 0;
          }
          if (gin.randomnumbers.found) {
            printWarning("Ignored md++ specific block RANDOMNUMBERS\n");
            gin.randomnumbers.found = 0;
          }
          if (gin.multistep.found) {
            printWarning("Ignored md++ specific block MULTISTEP\n");
            gin.multistep.found = 0;
          }
        } else { // Ignore promd specific blocks
          if (gin.consistencycheck.found) {
            printWarning("Ignored promd specific block CONSISTENCYCHECK\n");
            gin.consistencycheck.found = 0;
          }
          if (gin.thermostat.found) {
            if (gin.multibath.found) {
              printWarning("Ignored promd specific block THERMOSTAT\n");
              gin.thermostat.found = 0;
            } else {
              printError("You want to perform a md++ run and you have specified the PROMD specific\n"
                      "block THERMOSTAT. Maybe you want to specify the MULTIBATH block instead?\n");
            }
          }
          if (gin.barostat.found) {
            if (gin.pressurescale.found) {
              printWarning("Ignored promd specific block BAROSTAT\n");
              gin.barostat.found = 0;
            } else {
              printError("You want to perform a md++ run and you have specified the PROMD specific\n"
                      "block BAROSTAT. Maybe you want to specify the PRESSURESCALE block instead?\n");
            }
          }
          if (gin.virial.found) {
            printWarning("Ignored promd specific block VIRIAL\n");
            gin.virial.found = 0;
          }
          if (gin.overalltransrot.found) {
            if (gin.comtransrot.found) {
              printWarning("Ignored promd specific block OVERALLTRANSROT\n");
              gin.overalltransrot.found = 0;
            } else {
              printError("You want to perform a md++ run and you have specified the PROMD specific\n"
                      "block OVERALLTRANSROT. Maybe you want to specify the COMTRANSROT block instead?\n");
            }
          }
          if (gin.debug.found) {
            printWarning("Ignored promd specific block DEBUG\n");
            gin.debug.found = 0;
          }
          if (gin.geomconstraints.found) {
            if (gin.constraint.found) {
              printWarning("Ignored promd specific block GEOMCONSTRAINTS\n");
              gin.geomconstraints.found = 0;
            } else {
              printError("You want to perform a md++ run and you have specified the PROMD specific\n"
                      "block GEOMCONSTRAINTS. Maybe you want to specify the CONSTRAINT block instead?\n");
            }
          }
          if (gin.neighbourlist.found) {
            if (gin.pairlist.found) {
              printWarning("Ignored promd specific block NEIGHBOURLIST\n");
              gin.neighbourlist.found = 0;
            } else {
              printError("You want to perform a md++ run and you have specified the PROMD specific\n"
                      "block NEIGHBOURLIST. Maybe you want to specify the PAIRLIST block instead?\n");
            }
          }
          // APE: LOCALELEV should not be ignored by md++
          // if (gin.localelev.found) {
          //  printWarning("Ignored promd specific block LOCALELEV\n");
          //  gin.localelev.found = 0;
          //}
          if (gin.umbrella.found) {
            printWarning("Ignored promd specific block UMBRELLA\n");
            gin.umbrella.found = 0;
          }
          if (gin.pathint.found) {
            printWarning("Ignored promd specific block PATHINT\n");
            gin.pathint.found = 0;
          }
        }
        // Ignore md++ specific blocks if promd input is to be written
        // and vice versa (end)

        // And check if all compulsory blocks have been specified
        if (!gin.system.found) {
          printError("Could not find SYSTEM block\n");
        }
        if (!gin.step.found) {
          printError("Could not find STEP block\n");
        }
        if (!gin.boundcond.found) {
          printError("Could not find BOUNDCOND block\n");
        }
        if (!gin.force.found) {
          printError("Could not find FORCE block\n");
        }
        if (gromosXX && !gin.constraint.found)
          printError("Could not find CONSTRAINT block\n");
        if (!gromosXX && !gin.geomconstraints.found) {
          printError("Could not find GEOMCONSTRAINTS block\n");
        }
        if (!gromosXX && !gin.neighbourlist.found)
          printError("Could not find NEIGHBOURLIST block\n");

        if (gromosXX && !gin.pairlist.found)
          printError("Could not find PAIRLIST block\n");
        if (!gin.nonbonded.found)
          printError("Could not find NONBONDED block\n");

        // The input restrictions were already checked when reading the blocks
        // Now, do the cross checks and logical errors
        //
        // Andreas:
        // True, BUT:
        // ==========
        // The input restrictions are not checked completely when using a joblist
        // file since the function setParam and readJobinfo just read but do NOT
        // check the input again
        // => unallowed values in the generated input files
        //
        // to avoid double error checkin all checks (restrictions and cross checks)
        // should be done in the following (hundreds of) lines:
        //
        // (add a block in alphabetical order, otherwise it gets a complete mess!)
        //
        // BAROSTAT
        if (gin.barostat.found) {
          bool npcpl_ok = false;
          int npcpl[6];
          for (int i = 0; i < 6; i++) npcpl[i] = gin.barostat.npcpl[i];
          if (npcpl[0] == 0 && npcpl[1] == 0 && npcpl[2] == 0 &&
                  npcpl[3] == 0 && npcpl[2] == 0 && npcpl[2] == 0) npcpl_ok = true;
          if (npcpl[0] == 1 && npcpl[1] == 1 && npcpl[2] == 1 &&
                  npcpl[3] == 1 && npcpl[2] == 1 && npcpl[2] == 1) npcpl_ok = true;
          if (npcpl[0] == 1 && npcpl[1] == 1 && npcpl[2] == 1 &&
                  npcpl[3] == 0 && npcpl[2] == 0 && npcpl[2] == 0) npcpl_ok = true;
          if (npcpl[0] >= 0 && npcpl[0] <= 3 &&
                  npcpl[1] >= 0 && npcpl[1] <= 3 &&
                  npcpl[2] >= 0 && npcpl[2] <= 3 &&
                  npcpl[0] != npcpl[1] && npcpl[0] != npcpl[2]) npcpl_ok = true;
          if (!npcpl_ok) {
            std::stringstream ss;
            ss << "Combination of NPCPL variables {" << npcpl[0] << ","
                    << npcpl[1] << "," << npcpl[2] << "," << npcpl[3] << ","
                    << npcpl[4] << "," << npcpl[5] << "} is not allowed";
            printError(ss.str());
          }
          if (gin.energymin.found && gin.energymin.ntem != 0 &&
                  gin.barostat.ntp != 0)
            printError("Cannot do pressure coupling during an energy minimisation");
          if (gin.readtraj.found && gin.readtraj.ntrd != 0 &&
                  gin.barostat.ntp != 0)
            printError("Cannot do pressure coupling when reading in a trajectory from file");
          if (gin.virial.found && gin.virial.ntv == 0 && gin.barostat.ntp != 0)
            printError("If virial is not calculated (NTV=0 in VIRIAL block), you cannot do pressure coupling (NTP!=0 in BAROSTAT block");
          if (gin.boundcond.found && gin.boundcond.ntb == 0 &&
                  gin.barostat.ntp != 0)
            printError("No pressure coupling possible under vacuum boundary conditions");
          if (gin.boundcond.found && gin.boundcond.ntb == -1 &&
                  gin.barostat.ntp != 0)
            if (npcpl[0] != 1 || npcpl[1] != 1 || npcpl[2] != 1 ||
                    npcpl[3] != 0 || npcpl[4] != 0 || npcpl[5] != 0)
              printError("Pressure coupling with NTB=-1 in BOUNDCOND block requires in BAROSTAT block\n NTB=0 or NPCPL={1,1,1,0,0,0}");
          if (gin.boundcond.found && gin.boundcond.ntb == 1 &&
                  gin.barostat.ntp != 0)
            if (npcpl[3] != 0 || npcpl[4] != 0 || npcpl[5] != 0)
              printError("Pressure coupling with NTB= 1 in BOUNDCOND block requires in BAROSTAT block\n NTB=0 or NPCPL={i,j,k,0,0,0} with i,j,k=0,1,2,3");
          if (gin.boundcond.found && gin.boundcond.ntb == 2 &&
                  gin.barostat.ntp != 0)
            if (npcpl[3] != 0 || npcpl[4] != 0 || npcpl[5] != 0)
              if (npcpl[0] != 1 || npcpl[1] != 1 || npcpl[2] != 1 ||
                      npcpl[3] != 1 || npcpl[4] != 1 || npcpl[5] != 1)
                printError("Pressure coupling with NTB=-1 in BOUNDCOND block requires in BAROSTAT block\n NTB=0 or NPCPL={i,j,k,0,0,0} with i,j,k=0,1,2,3 or NPCPL={1,1,1,1,1,1}");
          if (gin.initialise.found && gin.initialise.ntinhb != 0 &&
                  gin.barostat.ntp != 3)
            printError("NTINHB!=0 in INITIALISE block requires NTP=3 in BAROSTAT block");
          if (gin.positionres.found && gin.positionres.ntpors == 1 &&
                  gin.barostat.ntp == 0)
            printError("NTPORS==1 in POSITIONRES block requires NTP!=0 in BAROSTAT block");
          if (gin.nonbonded.found &&
                  (gin.nonbonded.nlrele <= -2 || gin.nonbonded.nlrele >= 2) &&
                  gin.nonbonded.nqeval == 0 &&
                  gin.barostat.ntp != 0)
            printWarning("Pressure coupling with a lattice sum method to calculate the nonbonded interactions and NQEVAL=0 in NONBONDED block, is maybe not very wise");
        }
        //
        // BOUNDCOND block
        if (gin.boundcond.found) {
	  if (gin.boundcond.ntb < -1 || gin.boundcond.ntb > 2) {
	    std::stringstream ss;
	    ss << gin.boundcond.ntb;
	    printIO("BOUNDCOND", "NTB", ss.str(), "-1,0,1,2");
	  }
	  if (gin.boundcond.ndfmin < 0) {
	    std::stringstream ss;
	    ss << gin.boundcond.ndfmin;
	    printIO("BOUNDCOND", "NDFMIN", ss.str(), ">=0");
	  }
          if ((!gin.overalltransrot.found ||
                  (gin.overalltransrot.found && gin.overalltransrot.ncmro == 0)) &&
                  (!gin.stochdyn.found ||
                  (gin.stochdyn.found && gin.stochdyn.ntsd == 0)) &&
                  (!gin.energymin.found ||
                  (gin.energymin.found && gin.energymin.ntem == 0)) &&
                  (!gin.readtraj.found ||
                  (gin.readtraj.found && gin.readtraj.ntrd == 0)) &&
                  (!gin.consistencycheck.found ||
                  (gin.consistencycheck.found && gin.consistencycheck.ntchk == 0)) &&
                  gin.boundcond.ntb == 0)
            printWarning("Running a vacuum simulation (NTB=0 in BOUNDCOND block) without NCMRO in OVERALLTRANSROT block may be dangerous");
          if (l_coord && gin.boundcond.ntb != 0) {

            if (gin.boundcond.ntb == -1)
              if (!(crd.box.K().abs() == crd.box.L().abs() && crd.box.K().abs() == crd.box.M().abs()))
                //if (crd.box[1] != crd.box[2] || crd.box[1] != crd.box[3] ||
                //    crd.box[2] != crd.box[3])
                printError("NTB=-1 in BOUNDCOND means truncated octahedron, but the box-lenghts are not identical");
            double minbox = 1e6;
            double maxcutoff = 0;
            if (crd.box.boxformat() == gcore::Box::genbox) {
              double a, b, c, alpha, beta, gamma, triclinicvolume;
              a = crd.box.K().abs();
              b = crd.box.L().abs();
              c = crd.box.M().abs();
              alpha = acos(crd.box.L().dot(crd.box.M()) / (crd.box.L().abs() * crd.box.M().abs()));
              beta = acos(crd.box.K().dot(crd.box.M()) / (crd.box.K().abs() * crd.box.M().abs()));
              gamma = acos(crd.box.K().dot(crd.box.L()) / (crd.box.K().abs() * crd.box.L().abs()));

              triclinicvolume = a * b * c *
                      (1.0 - cos(alpha) * cos(alpha) - cos(beta) * cos(beta) - cos(gamma) * cos(gamma)
                      + 2.0 * cos(alpha) * cos(beta) * cos(gamma));

              minbox = min(triclinicvolume / (a * b * sin(gamma)), min(
                      triclinicvolume / (a * c * sin(beta)),
                      triclinicvolume / (b * c * sin(alpha))));
            }
            if (gin.boundcond.ntb == -1) minbox *= 0.5 * 1.732051;
            if (gin.pairlist.found) maxcutoff = gin.pairlist.rcutl;
            if (gin.neighbourlist.found) maxcutoff = gin.neighbourlist.rltwpl;
            if (minbox < 2 * maxcutoff) {
              std::stringstream ss;
              ss << "The largest cutoff in the system (" << maxcutoff
                      << ") is larger than half the\nmaximum size according to "
                      << "the box (" << minbox << ")";
              printError(ss.str());
            }
          }
        }
        //
        // CGRAIN
        if (gin.cgrain.found) {
	  if (gin.cgrain.ntcgran < 0 || gin.cgrain.ntcgran > 2) {
	    std::stringstream ss;
	    ss << gin.cgrain.ntcgran;
	    printIO("CGRAIN", "NTCGRAN", ss.str(), "0,1,2");
	  }
	  if (gin.cgrain.eps < 0.0) {
	    std::stringstream ss;
	    ss << gin.cgrain.eps;
	    printIO("CGRAIN", "EPS", ss.str(), ">= 0.0");
	  }
	  if (gin.cgrain.eps < 0) {
	    printError("EPS in CGRAIN block must be >= 0");
	  }
	  if (gin.cgrain.ntcgran < 0 || gin.cgrain.ntcgran > 2) {
	    printError("NTCGRAIN in CGRAIN block must be 0, 1 or 2");
	  }
	}
        //
        // COMTRANSROT
        if (gin.comtransrot.found) {
          if (gin.energymin.found && gin.energymin.ntem != 0 &&
                  gin.comtransrot.nscm != 0)
            printError("You cannot remove centre-of-mass motion (NSCM != 0 in COMTRANSROT block) in an energy minimisation (NTEM != 0 in ENERGYMIN block).");
          if (gin.readtraj.found && gin.readtraj.ntrd != 0 &&
                  gin.comtransrot.nscm != 0)
            printError("You cannot remove centre-of-mass motion (NSCM != 0 in COMTRANSROT block) when reading in the trajectory from file (NTRD != 0 in READTRAJ block");
        }
        //
        // CONSISTENCYCHECK block
        if (gin.consistencycheck.found) {
	  for (unsigned int i = 0; i < gin.consistencycheck.nckf.size(); i++) {
	    if (gin.consistencycheck.nckf[i] <= 0) {
	      std::stringstream si, ss;
	      si << "NCKF[" << i + 1 << "]";
	      ss << gin.consistencycheck.nckf[i];
	      printIO("CONSISTENCYCHECK", si.str(), ss.str(), "> 0");
	    }
	  }
	  if (gin.consistencycheck.ntchk != 0 && gin.consistencycheck.ntchk != 1) {
	    std::stringstream ss;
	    ss << gin.consistencycheck.ntchk;
	    printIO("CONSISTENCYCHECK", "NTCHK", ss.str(), "0,1");
	  }
	  if (gin.consistencycheck.ntckf != 0 && gin.consistencycheck.ntckf != 1) {
	    std::stringstream ss;
	    ss << gin.consistencycheck.ntckf;
	    printIO("CONSISTENCYCHECK", "NTCKF", ss.str(), "0,1");
	  }
	  if (gin.consistencycheck.fdckf <= 0.0) {
	    std::stringstream ss;
	    ss << gin.consistencycheck.fdckf;
	    printIO("CONSISTENCYCHECK", "NDCKF", ss.str(), "> 0.0");
	  }
	  if (gin.consistencycheck.ntckv != 0 && gin.consistencycheck.ntckv != 1) {
	    std::stringstream ss;
	    ss << gin.consistencycheck.ntckv;
	    printIO("CONSISTENCYCHECK", "NTCKV", ss.str(), "0,1");
	  }
	  if (gin.consistencycheck.fdckv <= 0.0) {
	    std::stringstream ss;
	    ss << gin.consistencycheck.fdckv;
	    printIO("CONSISTENCYCHECK", "NDCKV", ss.str(), "> 0.0");
	  }
	  if (gin.consistencycheck.ntckt != 0 && gin.consistencycheck.ntckt != 1) {
	    std::stringstream ss;
	    ss << gin.consistencycheck.ntckt;
	    printIO("CONSISTENCYCHECK", "NTCKT", ss.str(), "0,1");
	  }
	  if (gin.consistencycheck.ntcke != 0 && gin.consistencycheck.ntcke != 1) {
	    std::stringstream ss;
	    ss << gin.consistencycheck.ntcke;
	    printIO("CONSISTENCYCHECK", "NTCKE", ss.str(), "0,1");
	  }
	  if (gin.consistencycheck.ntckr != 0 && gin.consistencycheck.ntckr != 1) {
	    std::stringstream ss;
	    ss << gin.consistencycheck.ntckr;
	    printIO("CONSISTENCYCHECK", "NTCKR", ss.str(), "0,1");
	  }
	  if (gin.consistencycheck.ntckl != 0 && gin.consistencycheck.ntckl != 1) {
	    std::stringstream ss;
	    ss << gin.consistencycheck.ntckl;
	    printIO("CONSISTENCYCHECK", "NTCKL", ss.str(), "0,1");
	  }
	  if (gin.consistencycheck.fdckl <= 0.0) {
	    std::stringstream ss;
	    ss << gin.consistencycheck.fdckl;
	    printIO("CONSISTENCYCHECK", "NDCKF", ss.str(), "> 0.0");
	  }
	  if (gin.boundcond.found && gin.boundcond.ntb != 1 &&
		  gin.consistencycheck.ntchk != 0 && gin.consistencycheck.ntckt != 0)
	    printError("NTCHK!=0 and NTCKT!=0 in CONSISTENCYCHECK block is only allowed for rectangular boundary conditions (NTB=1 in BOUNDCOND)");
	  if (gin.boundcond.found && gin.boundcond.ntb == 1 &&
		  gin.consistencycheck.ntckv != 0)
	    printWarning("if NTCKV!=0 in CONSISTENCYCHECK block, the boundary conditions will be re-set internally to triclinic");
	  if (gin.boundcond.found && gin.boundcond.ntb == 0 &&
		  gin.consistencycheck.ntckv == 1)
	    printWarning("NTCKV=1 in CONSISTENCYCHECK block and NTB=0 in BOUNDCOND block, this means thatFDCKV refers to a hypothetical cubic box with edge 1 nm");
	  int nckf_prev = 0;
	  for (unsigned int i = 0; i < gin.consistencycheck.nckf.size(); i++) {
	    if (gin.consistencycheck.nckf[i] >= numTotalAtoms) {
	      std::stringstream ss;
	      ss << "NCKF[" << i + 1 << "] is out of the range 1..NATTOT ("
		      << numTotalAtoms << ")";
	      printError(ss.str());
	    }
	    if (gin.consistencycheck.nckf[i] <= nckf_prev) {
	      printError("NCKF in CONSISTENCYCHECK block should come in ascending order");
	    }
	    nckf_prev = gin.consistencycheck.nckf[i];
	  }
	  if (gin.consistencycheck.ntckf == 1 &&
		  gin.consistencycheck.nckf.size() == 0)
	    printError("If NTCKF=1 in CONSISTENCYCHECK, you should also specify atoms");
	  if (gin.energymin.found && gin.energymin.ntem != 0 &&
		  gin.consistencycheck.ntchk != 0)
	    printError("You cannot use CONSISTENCYCHECK with a minimisation");
	  if (gin.stochdyn.found && gin.stochdyn.ntsd != 0 &&
		  gin.consistencycheck.ntchk != 0)
	    printError("You cannot use CONSISTENCYCHECK with Stochastic Dynamics");
	  if (gin.readtraj.found && gin.readtraj.ntrd != 0 &&
		  gin.consistencycheck.ntchk != 0)
	    printError("You cannot use CONSISTENCYCHECK when readin a trajectory");
	  if (gin.step.found && gin.step.nstlim != 1 &&
		  gin.consistencycheck.ntchk != 0)
	    printError("Doing a CONSISTENCYCHECK requires NSTLIM=1 in the STEP block");
	  if (gin.thermostat.found && gin.thermostat.ntt != 0 &&
		  gin.consistencycheck.ntchk != 0)
	    printError("CONSISTENCYCHECK is not implemented with NTT!=0 in THERMOSTAT block");
	  if (gin.virial.found && gin.virial.ntv == 0 &&
		  !(gin.consistencycheck.ntchk == 0 || gin.consistencycheck.ntckv == 0))
	    printError("if NTV=0 in VIRIAL block, you cannot do a virial CONSISTENCYCHECK. Set NTCHK=0 or NTCKV=0");
	  if (gin.overalltransrot.found && gin.overalltransrot.ncmtr != 0 &&
		  gin.consistencycheck.ntchk != 0)
	    printError("You cannot use CONSISTENCYCHECK when NCMTR!=0 in OVERALLTRANSROT block");
	  if (gin.overalltransrot.found && gin.overalltransrot.ncmro != 0 &&
		  gin.consistencycheck.ntchk != 0)
	    printError("You cannot use CONSISTENCYCHECK when NCMRO!=0 in OVERALLTRANSROT block");
	  if (gin.nonbonded.found && abs(gin.nonbonded.nlrele) == 2 &&
		  !(gin.consistencycheck.ntchk == 0 || gin.consistencycheck.ntcke == 0))
	    printError("if ABS(NLRELE)=2 in NONBONDED block, you cannot do CONSISTENCYCHECK. Set NTCHK=0 or NTCKE=0");
	  if (gin.nonbonded.found && gin.nonbonded.nlrele <= -2 &&
		  !(gin.consistencycheck.ntchk == 0 || gin.consistencycheck.ntckr == 0))
	    printError("if NLRELE <= -2 in NONBONDED block, you cannot do CONSISTENCYCHECK. Set NTCHK=0 or NTCKR=0");
	  if (gin.nonbonded.found && gin.nonbonded.na2clc == 4 &&
		  abs(gin.nonbonded.nlrele) < 3 &&
		  (gin.consistencycheck.ntchk != 1 || gin.consistencycheck.ntcke != 1))
	    printError("if NA2CLC=4 and ABS(NLRELE) <= 2 in NONBONDED block, then you have give NTCHK=1 and NTCKE=1 in the CONSISTENCYCHECK block");
	  if (gin.nonbonded.found &&
		  (gin.nonbonded.nrdgrd != 0 || gin.nonbonded.nwrgrd != 0) &&
		  gin.consistencycheck.ntchk != 0)
	    printError("You cannot use CONSISTENCYCHECK when NRDGRD!=0 or NWRGRD!=0 in NONBONDED block");
	  if ((!gin.perturbation.found ||
		  (gin.perturbation.found && gin.perturbation.ntg == 0)) &&
		  gin.consistencycheck.ntchk != 0 && gin.consistencycheck.ntckl != 0)
	    printError("You cannot give NTCHK!=0 and NTCKL!=0 in CONSISTENCYCHECK block with NTG=0 in PERTURBATION block");
	}
	//
	// CONSTRAINT block
	if (gin.constraint.found) {

	  if (gin.constraint.ntc < 0 || gin.constraint.ntc > 4) {
	    stringstream ss;
	    ss << gin.constraint.ntc;
	    printIO("CONSTRAINT", "NTC", ss.str(), "0..4");
	  }
	  if (gin.constraint.ntcp < 1 || gin.constraint.ntcp > 3) {
	    stringstream ss;
	    ss << gin.constraint.ntcp;
	    printIO("CONSTRAINT", "NTCP", ss.str(), "1..3");
	  }
	  // number of option parameters for NTCP
	  int p = 1;
	  if(gin.constraint.ntcp == 3){
	    p = 3;
	  }
	  for (int i = 0; i < p; ++i) {
	    if (gin.constraint.ntcp0[p] < 0) {
	      stringstream ss, is;
	      is << "NTCP0[" << i + 1 << "]";
	      ss << gin.constraint.ntcp0[i];
	      printIO("CONSTRAINT", is.str(), ss.str(), ">= 0");
	    }
	  }
	  if (gin.constraint.ntcs < 1 || gin.constraint.ntcs > 4) {
	    stringstream ss;
	    ss << gin.constraint.ntcs;
	    printIO("CONSTRAINT", "NTCS", ss.str(), "1..4");
	  }
	  // number of option parameters for NTCS
	  p = 1;
	  if (gin.constraint.ntcs == 3) {
	    p = 3;
	  } else if (gin.constraint.ntcs == 4) {
	    p = 0;
	  }
	  for (int i = 0; i < p; ++i) {
	    if (gin.constraint.ntcs0[i] < 0) {
	      stringstream ss, is;
	      is << "NTCS0[" << i + 1 << "]";
	      ss << gin.constraint.ntcs0[i];
	      printIO("CONSTRAINT", is.str(), ss.str(), ">= 0");
	    }
	  }
          if (gin.system.found && gin.system.npm == 0 && gin.constraint.ntc != 1)
            printError("No solute molecules (NPM=0 in SYSTEM block), what do you want to constrain (NTC!=0 in CONSTRAINT block)");

          if ((gin.constraint.ntc == 1 && gin.step.dt > 0.0005) ||
                  (gin.constraint.ntc == 2 && gin.step.dt > 0.001) ||
                  (gin.constraint.ntc == 3 && gin.step.dt > 0.002) ||
                  (gin.constraint.ntc == 4 && gin.step.dt > 0.0005)) {
            ostringstream os;
            string comment;
            double suggest = 0.0005;
            if (gin.constraint.ntc == 1) {
              comment = "no constraints on solute";
              suggest = 0.0005;
            } else if (gin.constraint.ntc == 2) {
              comment = "constraints on bonds with H";
              suggest = 0.001;
            } else if (gin.constraint.ntc == 3) {
              comment = "constraints on all bonds";
              suggest = 0.002;
            } else if (gin.constraint.ntc == 4) {
              comment = "constraints on some bonds, and no constraints on other";
              suggest = 0.0005;
            }
            os << "DT in STEP block is set to " << gin.step.dt << ", which is "
                    << "considered to be too large\n";
            os << "if NTC = " << gin.constraint.ntc
                    << " in CONSTRAINT block.\n";
            os << "For NTC = " << gin.constraint.ntc << " (" << comment
                    << ") rather "
                    << "use DT = " << suggest << ".\n";
            printWarning(os.str());
          }
        }
        //
        // COVALENTFORM block
        if (gin.covalentform.found) {
	  if (gin.covalentform.ntbbh != 0 && gin.covalentform.ntbbh != 1) {
	    std::stringstream ss;
	    ss << gin.covalentform.ntbbh;
	    printIO("COVALENTFORM", "NTBBH", ss.str(), "0,1");
	  }
	  if (gin.covalentform.ntbah != 0 && gin.covalentform.ntbah != 1) {
	    std::stringstream ss;
	    ss << gin.covalentform.ntbah;
	    printIO("COVALENTFORM", "NTBAH", ss.str(), "0,1");
	  }
	  if (gin.covalentform.ntbdn != 0 && gin.covalentform.ntbdn != 1) {
	    std::stringstream ss;
	    ss << gin.covalentform.ntbdn;
	    printIO("COVALENTFORM", "NTBDN", ss.str(), "0,1");
	  }
          bool bondtype = false;
          bool bondstretchtype = false;
          bool bondangletype = false;
          bool bondanglebendtype = false;
          bool dihedraltype = false;
          bool torsdihedraltype = false;
          for (unsigned int i = 0; i < topo.blocks.size(); i++) {
            if (topo.blocks[i] == "BONDTYPE" && topo.blockslength[i] != 0)
              bondtype = true;
            if (topo.blocks[i] == "BONDSTRETCHTYPE" && topo.blockslength[i] != 0)
              bondstretchtype = true;
            if (topo.blocks[i] == "BONDANGLETYPE" && topo.blockslength[i] != 0)
              bondangletype = true;
            if (topo.blocks[i] == "BONDANGLEBENDTYPE" && topo.blockslength[i] != 0)
              bondanglebendtype = true;
            if (topo.blocks[i] == "DIHEDRALTYPE" && topo.blockslength[i] != 0)
              dihedraltype = true;
            if (topo.blocks[i] == "TORSDIHEDRALTYPE" && topo.blockslength[i] != 0)
              torsdihedraltype = true;
          }
          if (bondtype && !bondstretchtype && gin.covalentform.ntbbh != 0) {
            printError("Topology has only a BONDTYPE block, this means that NTBBH in COVALENTFORM block should be 0");
          }
          if (bondangletype && !bondanglebendtype && gin.covalentform.ntbah != 0) {
            printError("Topology has only a BONDANGLETYPE block, this means that NTBAH in COVALENTFORM block should be 0");
          }
          if (dihedraltype && !torsdihedraltype && gin.covalentform.ntbdn != 1)
            printError("Topology has only DIHEDRALTYPE block, this means that NTBDN in COVALENTFORM block should be 1");
          if (torsdihedraltype && gin.covalentform.ntbdn == 1)
            printWarning("NTBDN=1 in COVALENTFORM block and you give a TORSDIHEDRALTYPE block. Make sure that the phaseshifts are 0 or 180 degree");
          if (bondtype && bondstretchtype &&
                  (gin.force.ntf[0] != 0 || gin.force.ntf[1] != 0))
            printWarning("BONDTYPE block in topology will be ignored in favour of BONDSTRETCHTYPE block");
          if (bondangletype && bondanglebendtype &&
                  (gin.force.ntf[2] != 0 || gin.force.ntf[3] != 0))
            printWarning("BONDANGLETYPE block in topology will be ignored in favour of BONDANGLEBENDTYPE block");
          if (dihedraltype && torsdihedraltype &&
                  (gin.force.ntf[5] != 0 || gin.force.ntf[6] != 0))
            printWarning("DIHEDRALTYPE block in topology will be ignored in favour of TORSDIHEDRALTYPE block");
	}
	//
	// DEBUG
	if (gin.debug.found) {
	  for (unsigned int i = 0; i < gin.debug.routines.size(); ++i) {
	    if (gin.debug.routines[i].iiideo < 0) {
	      std::stringstream ss, si;
	      ss << gin.debug.routines[i].iiideo;
	      si << "IIIDEO[" << i + 1 << "]";
	      printIO("DEBUG", si.str(), ss.str(), "> 0");
	    }
	  }
	}
	//
	// DIHEDRALRES block
	if (gin.dihedralres.found) {
	  if (gin.dihedralres.ntdlr < 0 || gin.dihedralres.ntdlr > 3) {
	    std::stringstream ss;
	    ss << gin.dihedralres.ntdlr;
	    printIO("DIHEDRALRES", "NTDLR", ss.str(), "0,1,2,3");
	  }
	  if (gin.dihedralres.cdlr < 0.0) {
	    std::stringstream ss;
	    ss << gin.dihedralres.cdlr;
	    printIO("DIHEDRALRES", "CDLR", ss.str(), ">= 0.0");
	  }
          if (gin.dihedralres.ntdlr != 0) {
            if (!l_dihres)
              printError("Dihedral restraining specified in DIHEDRALRES block, but no DIHEDRALRESSPEC file specified");
            else {
              fileInfo dihres;
              Ginstream idihres(s_dihres);
              idihres >> dihres;
              idihres.close();
              bool dihedralresspec = false;
              for (unsigned int i = 0; i < dihres.blocks.size(); i++)
                if (dihres.blocks[i] == "DIHEDRALRESSPEC")
                  dihedralresspec = true;
              if (dihedralresspec == false)
                printError("Dihedral restraining specified in DIHEDRALRES block, but specification file does not contain a DIHEDRALRESSPEC block");
            }
          }
        }
        //
	// DISTANCERES block
	if (gin.distanceres.found) {
	  if (gin.distanceres.ntdir < -2 || gin.distanceres.ntdir > 3) {
	    std::stringstream ss;
	    ss << gin.distanceres.ntdir;
	    printIO("DISTANCERES", "NTDIR", ss.str(), "-2,-1,0,1,2,3");
	  }
	  if (gin.distanceres.ntdira != 0 && gin.distanceres.ntdira != 1) {
	    std::stringstream ss;
	    ss << gin.distanceres.ntdir;
	    printIO("DISTANCERES", "NTDIRA", ss.str(), "0,1");
	  }
	  if (gin.distanceres.cdir < 0.0) {
	    std::stringstream ss;
	    ss << gin.distanceres.cdir;
	    printIO("DISTANCERES", "CDIR", ss.str(), ">= 0.0");
	  }
	  if (gin.distanceres.dir0 < 0.0) {
	    std::stringstream ss;
	    ss << gin.distanceres.dir0;
	    printIO("DISTANCERES", "DIR0", ss.str(), ">= 0.0");
	  }
	  if (gin.distanceres.taudir < 0.0) {
	    std::stringstream ss;
	    ss << gin.distanceres.taudir;
	    printIO("DISTANCERES", "TAUDIR", ss.str(), ">= 0.0");
	  }
          if (gin.distanceres.ntdira == 1 && gin.distanceres.ntdir >= 0)
            printError("NTDIRA=1 in DISTANCERES block is only allowed if NTDIR < 0");
          if (gin.distanceres.ntdir != 0) {
            if (!l_disres)
              printError("Distance restraining specified in DISTANCERES block, but no DISTANCERESSPEC file specified");
            else {
              fileInfo disres;
              Ginstream idisres(s_disres);
              idisres >> disres;
              idisres.close();
              bool distanceresspec = false;
              for (unsigned int i = 0; i < disres.blocks.size(); i++)
		if (disres.blocks[i] == "DISTANCERESSPEC")
		  distanceresspec = true;
	      if (distanceresspec == false)
		printError("Distance restraining specified in DISTANCERES block, but specification file does not contain a DISTANCERESSPEC block");
	    }
	  }
	  if (gin.distanceres.ntdira == 1) {
	    bool disresexpave = false;
	    for (unsigned int i = 0; i < crd.blocks.size(); i++)
	      if (crd.blocks[i] == "DISRESEXPAVE") disresexpave = true;
	    if (disresexpave == false)
	      printError("Distance averages to be read from startup file (NTDIRA=1 in DISTANCERES block), but no DISRESEXPAVE block in coordinate file");
	  }
	  if (gin.distanceres.ntdir < 0 && gin.distanceres.ntdira == 0)
	    printWarning("Time-averaged distance restraining, but no reading of initial averages from file. This is unwise if you are doing a continuation run");
	}
	//
	// ENERGYMIN block
	if (gin.energymin.found) {
	  if (gin.energymin.ntem < 0 || gin.energymin.ntem > 2) {
	    std::stringstream ss;
	    ss << gin.energymin.ntem;
	    printIO("ENERGYMIN", "NTEM", ss.str(), "0,1,2");
	  }
	  if (gin.energymin.ncyc <= 0) {
	    std::stringstream ss;
	    ss << gin.energymin.ncyc;
	    printIO("ENERGYMIN", "NCYC", ss.str(), "> 0");
	  }
	  if (gin.energymin.dele <= 0.0) {
	    std::stringstream ss;
	    ss << gin.energymin.dele;
	    printIO("ENERGYMIN", "DELE", ss.str(), "> 0.0");
	  }
	  if (gin.energymin.dx0 <= 0.0) {
	    std::stringstream ss;
	    ss << gin.energymin.dx0;
	    printIO("ENERGYMIN", "DX0", ss.str(), "> 0.0");
	  }
	  if (gin.energymin.dxm <= 0.0) {
	    std::stringstream ss;
	    ss << gin.energymin.dxm;
	    printIO("ENERGYMIN", "DXM", ss.str(), "> 0.0");
	  }
	  if (gin.energymin.nmin <= 0) {
	    std::stringstream ss;
	    ss << gin.energymin.nmin;
	    printIO("ENERGYMIN", "NMIN", ss.str(), "> 0");
	  }
	  if (gin.energymin.flim < 0.0) {
	    std::stringstream ss;
	    ss << gin.energymin.flim;
	    printIO("ENERGYMIN", "FLIM", ss.str(), ">= 0.0");
	  }
	}
        if (gin.energymin.found && gin.energymin.ntem != 0) {
          if (gin.energymin.dx0 > gin.energymin.dxm)
            printError("DX0 should be <= DXM in ENERGYMIN block");
          if (gin.stochdyn.found && gin.stochdyn.ntsd != 0)
            printError("NTEM!=0 in ENERGYMIN block does not work with NTSD!=0 in STOCHDYN block");
          if (gin.readtraj.found && gin.readtraj.ntrd != 0)
            printError("NTEM!=0 in ENERGYMIN block does not work with NTRD!=0 in READTRAJ block");
          if (gin.virial.found && gin.virial.ntv != 0)
            printError("NTEM!=0 in ENERGYMIN block does not work with NTV!=0 in VIRIAL block");
          if (gin.thermostat.found && gin.thermostat.ntt != 0)
            printError("NTEM!=0 in ENERGYMIN block does not work with NTT!=0 in THERMOSTAT block");
          if (gin.multibath.found && gin.multibath.temp0.size())
            printWarning("NTEM!=0 in ENERGYMIN block will ignore the temperature coupling");
          if (gin.barostat.found && gin.barostat.ntp != 0)
            printError("NTEM!=0 in ENERGYMIN block does not work with NTP!=0 in BAROSTAT block");
          if (gin.pressurescale.found && gin.pressurescale.couple != 0)
            printError("NTEM!=0 in ENERGYMIN block does not work with COUPLE!=off(0) in PRESSURESCALE block");
          if (gin.overalltransrot.found && gin.overalltransrot.ncmtr != 0)
            printError("NTEM!=0 in ENERGYMIN block does not work with NCMTR!=0 in OVERALLTRANSROT block");
          if (gin.overalltransrot.found && gin.overalltransrot.ncmro != 0)
            printError("NTEM!=0 in ENERGYMIN block does not work with NCMRO!=0 in OVERALLTRANSROT block");
          if (gin.localelev.found && gin.localelev.ntles != 0)
            printError("NTEM!=0 in ENERGYMIN block does not work with NTLES!=0 in LOCALELEV block");
          if (gin.initialise.found && gin.initialise.ntivel != 0)
            printError("Boltzmann weighted velocities cannot be generated for an energy minimisation, set NTIVEL=0 in INITIALISE block");
          if (gin.initialise.found && gin.initialise.ntishk > 1)
            printError("You cannot shake velocities (NTISHK > 1 in INITIALISE block) for a minimisation");
        }
        //
        // FORCE block
        if (gin.force.found) {
	  for (int i = 0; i < 10; i++) {
	    if (gin.force.ntf[i] != 0 && gin.force.ntf[i] != 1) {
	      std::stringstream ss, si;
	      ss << gin.force.ntf[i];
	      si << "NTF[" << i + 1 << "]";
	      printIO("FORCE", si.str(), ss.str(), "0,1");
	    }
	  }
	  for (unsigned int i = 0; i < gin.force.nre.size(); i++) {
	    if (gin.force.nre[i] <= 0) {
	      std::stringstream ss, si;
	      si << "NRE[" << i + 1 << "]";
	      ss << gin.force.nre[i];
	      printIO("FORCE", si.str(), ss.str(), "> 0");
	    }
	  }
          int cnt = 0;
          for (unsigned int i = 0; i < gin.force.nre.size(); i++) {
            if (gin.force.nre[i] <= cnt) {
              std::stringstream ss;
              ss << "Values of NRE in FORCE block not in ascending order:\n"
                      << gin.force.nre[i] << " <= " << cnt;
              printError(ss.str());
            }
          }
          if (gin.force.nre.size() &&
                  gin.force.nre[gin.force.nre.size() - 1] != numTotalAtoms) {
            std::stringstream ss;
            ss << "Last value of NRE in FORCE block should be equal to total"
                    << " number of atoms:"
                    << gin.force.nre[gin.force.nre.size() - 1] << " != "
                    << numTotalAtoms;
            printError(ss.str());
          }
          if ((!gin.nonbonded.found ||
                  (gin.nonbonded.found && gin.nonbonded.nlrele == 0)) &&
                  gin.force.ntf[8] == 1)
            printError("No nonbonded calculation (NLRELE=0 in NONBONDED block) does not agree with NTF[9]=1 in FORCE block");
          if (gin.nonbonded.found && gin.nonbonded.nlrele != 0 &&
                  gin.force.ntf[8] == 0)
            printError("Nonbonded calculation specified (NLRELE!=0 in NONBONDED block), but NTF[9]=0 in FORCE block");
          if (gin.nonbonded.found && gin.nonbonded.nlrlj != 0 &&
                  gin.force.ntf[9] == 0)
            printError("NLRLJ in NONBONDED requests a LJ correction, but LJ not calculated according to NTF[10] in FORCE block");
          bool nbonh = false, nbon = false, ntheh = false, nthe = false;
          bool nqhih = false, nqhi = false, nphih = false, nphi = false;
          for (unsigned int i = 0; i < topo.blocks.size(); i++) {
            if (topo.blocks[i] == "BONDH" && topo.blockslength[i] >= 1)
              nbonh = true;
            if (topo.blocks[i] == "BOND" && topo.blockslength[i] >= 1)
              nbon = true;
            if (topo.blocks[i] == "BONDANGLEH" && topo.blockslength[i] >= 1)
              ntheh = true;
            if (topo.blocks[i] == "BONDANGLE" && topo.blockslength[i] >= 1)
              nthe = true;
            if (topo.blocks[i] == "IMPDIHEDRALH" && topo.blockslength[i] >= 1)
              nqhih = true;
            if (topo.blocks[i] == "IMPDIHEDRAL" && topo.blockslength[i] >= 1)
              nqhi = true;
            if (topo.blocks[i] == "DIHEDRALH" && topo.blockslength[i] >= 1)
              nphih = true;
            if (topo.blocks[i] == "DIHEDRAL" && topo.blockslength[i] >= 1)
              nphi = true;
          }
          if (gin.force.ntf[0] == 1 && nbonh == false)
            printError("NTF[1]=1 in FORCE block, but no bonds in BONDH block");
          if (gin.force.ntf[1] == 1 && nbon == false)
            printError("NTF[2]=1 in FORCE block, but no bonds in BOND block");
          if (gin.force.ntf[2] == 1 && ntheh == false)
            printError("NTF[3]=1 in FORCE block, but no angles in BONDANGLEH block");
          if (gin.force.ntf[3] == 1 && nthe == false)
            printError("NTF[4]=1 in FORCE block, but no angles in BONDANGLE block");
          if (gin.force.ntf[4] == 1 && nqhih == false)
            printError("NTF[5]=1 in FORCE block, but no impropers in IMPDIHEDRALH block");
          if (gin.force.ntf[5] == 1 && nqhi == false)
            printError("NTF[6]=1 in FORCE block, but no impropers in IMPDIHEDRAL block");
          if (gin.force.ntf[6] == 1 && nphih == false)
            printError("NTF[7]=1 in FORCE block, but no dihedrals in DIHEDRALH block");
          if (gin.force.ntf[7] == 1 && nphi == false)
            printError("NTF[8]=1 in FORCE block, but no dihedrals in DIHEDRAL block");
          if (gin.geomconstraints.found && gin.geomconstraints.ntcph == 1 &&
                  !nbonh)
            printError("NTCPH=1 in GEOMCONSTRAINTS block, but no bonds in BONDH block");
          if (gin.geomconstraints.found && gin.geomconstraints.ntcpn == 1 &&
                  !nbon)
            printError("NTCPN=1 in GEOMCONSTRAINTS block, but no bonds in BOND block");
          if (gin.constraint.found && gin.constraint.ntc == 2 &&
                  !nbonh)
            printError("NTC==2 in CONSTRAINT block, but no bonds in BONDH block");
          if (gin.constraint.found && gin.constraint.ntc > 2 &&
                  !(nbonh || nbon))
            printError("NTC>2 in CONSTRAINT block, but no bonds in BOND or BONDH block");
          if (gin.geomconstraints.found && gin.geomconstraints.ntcph == 1 &&
                  nbonh && gin.force.ntf[0] == 1)
            printWarning("NTF[1]=1 in FORCE block, but bond lengths are constraint");
          if (gin.geomconstraints.found && gin.geomconstraints.ntcpn == 1 &&
                  nbon && gin.force.ntf[1] == 1)
            printWarning("NTF[2]=1 in FORCE block, but bond lengths are constraint");
          if (gin.constraint.found && gin.constraint.ntc == 2 &&
                  nbon && gin.force.ntf[0])
            printWarning("NTF[1]=1 in FORCE block, but bond lengths are constraint");
          if (gin.constraint.found && gin.constraint.ntc > 2 &&
                  nbon && gin.force.ntf[1])
            printWarning("NTF[2]=1 in FORCE block, but bond lengths are constraint");
          if (gin.geomconstraints.found && gin.geomconstraints.ntcph == 0 &&
                  nbonh && gin.force.ntf[0] == 0)
            printWarning("NTF[1]=0 in FORCE block, and bond lengths are not constraint");
          if (gin.geomconstraints.found && gin.geomconstraints.ntcpn == 0 &&
                  nbon && gin.force.ntf[1] == 0)
            printWarning("NTF[2]=0 in FORCE block, and bond lengths are not constraint");
          if (gin.constraint.found && gin.constraint.ntc < 2 &&
                  nbon && gin.force.ntf[0] == 0)
            printWarning("NTF[1]=0 in FORCE block, and bond lengths are not constraint");
          if (gin.constraint.found && gin.constraint.ntc > 3 &&
                  nbon && gin.force.ntf[2] == 0)
            printWarning("NTF[2]=0 in FORCE block, and bond lengths are not constraint");
          if (ntheh && gin.force.ntf[2] == 0)
            printWarning("NTF[3]=0 in FORCE block, but BONDANGLEH block not empty");
          if (nthe && gin.force.ntf[3] == 0)
            printWarning("NTF[4]=0 in FORCE block, but BONDANGLE block not empty");
          if (nqhih && gin.force.ntf[4] == 0)
            printWarning("NTF[3]=0 in FORCE block, but IMPDIHEDRALH block not empty");
          if (nqhi && gin.force.ntf[5] == 0)
            printWarning("NTF[4]=0 in FORCE block, but IMPDIHEDRAL block not empty");
          if (nphih && gin.force.ntf[6] == 0)
            printWarning("NTF[3]=0 in FORCE block, but DIHEDRALH block not empty");
          if (nphi && gin.force.ntf[7] == 0)
            printWarning("NTF[4]=0 in FORCE block, but DIHEDRAL block not empty");
          if (gin.neighbourlist.found &&
                  (gin.neighbourlist.nuirin != 0 || gin.neighbourlist.nusrin != 0) &&
                  gin.force.ntf[8] == 0 && gin.force.ntf[9] == 0)
            printWarning("According to NEIGHBOURLIST block you want to make a pairlist, but NTF[9] = NTF[10] = 0 in FORCE block");
        }
        //
        // GEOMCONSTRAINTS
	if (gin.geomconstraints.found) {
	  if (gin.geomconstraints.ntcph != 0 && gin.geomconstraints.ntcph != 1) {
	    std::stringstream ss;
	    ss << gin.geomconstraints.ntcph;
	    printIO("GEOMCONSTRAINTS", "NTCPH", ss.str(), "0,1");
	  }
	  if (gin.geomconstraints.ntcpn != 0 && gin.geomconstraints.ntcpn != 1) {
	    std::stringstream ss;
	    ss << gin.geomconstraints.ntcpn;
	    printIO("GEOMCONSTRAINTS", "NTCPN", ss.str(), "0,1");
	  }
	  if (gin.geomconstraints.ntcs != 0 && gin.geomconstraints.ntcs != 1) {
	    std::stringstream ss;
	    ss << gin.geomconstraints.ntcs;
	    printIO("GEOMCONSTRAINTS", "NTCS", ss.str(), "0,1");
	  }
	  if (gin.geomconstraints.shktol <= 0.0) {
	    std::stringstream ss;
	    ss << gin.geomconstraints.shktol;
	    printIO("GEOMCONSTRAINTS", "SHKTOL", ss.str(), "> 0.0");
	  }
          if (gin.system.found && gin.system.npm != 0 &&
                  (gin.geomconstraints.ntcph == 1 || gin.geomconstraints.ntcpn == 1))
            printError("NPM=0 in SYSTEM block, so what do you want to shake in GEOMCONSTRAINTS block?");
	}
	//
	// GROMOS96COMPAT
	if (gin.gromos96compat.found) {
	  if (gin.gromos96compat.ntnb96 != 0 && gin.gromos96compat.ntnb96 != 1) {
	    std::stringstream ss;
	    ss << gin.gromos96compat.ntnb96;
	    printIO("GROMOS96COMPAT", "NTNB96", ss.str(), "0,1");
	  }
	  if (gin.gromos96compat.ntr96 != 0 && gin.gromos96compat.ntr96 != 1) {
	    std::stringstream ss;
	    ss << gin.gromos96compat.ntr96;
	    printIO("GROMOS96COMPAT", "NTR96", ss.str(), "0,1");
	  }
	  if (gin.gromos96compat.ntp96 != 0 && gin.gromos96compat.ntp96 != 1) {
	    std::stringstream ss;
	    ss << gin.gromos96compat.ntp96;
	    printIO("GROMOS96COMPAT", "NTP96", ss.str(), "0,1");
	  }
	  if (gin.gromos96compat.ntg96 != 0 && gin.gromos96compat.ntg96 != 1) {
	    std::stringstream ss;
	    ss << gin.gromos96compat.ntg96;
	    printIO("GROMOS96COMPAT", "NTG96", ss.str(), "0,1");
	  }
	}
	//
  // INITIALISE
        if (gin.initialise.found) {
	  if (gin.initialise.ntivel != 0 && gin.initialise.ntivel != 1) {
	    std::stringstream ss;
	    ss << gin.initialise.ntivel;
	    printIO("INITIALISE", "NTIVEL", ss.str(), "0,1");
	  }
	  if (gin.initialise.ntishk < 0 || gin.initialise.ntivel > 3) {
	    std::stringstream ss;
	    ss << gin.initialise.ntishk;
	    printIO("INITIALISE", "NTISHK", ss.str(), "0,1,2,3");
	  }
	  if (gin.initialise.ntinht != 0 && gin.initialise.ntinht != 1) {
	    std::stringstream ss;
	    ss << gin.initialise.ntinht;
	    printIO("INITIALISE", "NTINHT", ss.str(), "0,1");
	  }
	  if (gin.initialise.ntinhb != 0 && gin.initialise.ntinhb != 1) {
	    std::stringstream ss;
	    ss << gin.initialise.ntinhb;
	    printIO("INITIALISE", "NTINHB", ss.str(), "0,1");
	  }
	  if (gin.initialise.ntishi != 0 && gin.initialise.ntishi != 1) {
	    std::stringstream ss;
	    ss << gin.initialise.ntishi;
	    printIO("INITIALISE", "NTISHI", ss.str(), "0,1");
	  }
	  if (gin.initialise.ntirtc != 0 && gin.initialise.ntirtc != 1) {
	    std::stringstream ss;
	    ss << gin.initialise.ntirtc;
	    printIO("INITIALISE", "NTIRTC", ss.str(), "0,1");
	  }
	  if (gin.initialise.nticom < 0 || gin.initialise.nticom > 3) {
	    std::stringstream ss;
	    ss << gin.initialise.nticom;
	    printIO("INITIALISE", "NTICOM", ss.str(), "0,1,2,3");
	  }
	  if (gin.initialise.ntisti != 0 && gin.initialise.ntisti != 1) {
	    std::stringstream ss;
	    ss << gin.initialise.ntisti;
	    printIO("INITIALISE", "NTISTI", ss.str(), "0,1");
	  }
	  if (gin.initialise.ig <= 0) {
	    std::stringstream ss;
	    ss << gin.initialise.ig;
	    printIO("INITIALISE", "IG", ss.str(), "> 0");
	  }
	  if (gin.initialise.tempi < 0.0) {
	    std::stringstream ss;
	    ss << gin.initialise.tempi;
	    printIO("INITIALISE", "TEMPI", ss.str(), ">= 0");
	  }
          if (gin.boundcond.found && gin.boundcond.ntb != 0 &&
                  (gin.initialise.nticom != 0 || gin.initialise.nticom > 1))
            printError("NTICOM=0,1 in INITISALISE block is only allowed for vacuum boundary conditions (NTB=0 in BOUNDCOND block");
          if (gin.boundcond.found && gin.boundcond.ntb == 0 &&
                  gin.initialise.ntishi != 0)
            printError("NTISHI!=0 in INITIALISE block is not allowed for vacuum boundary conditions (NTB!=0 in BOUNDCOND block)");
          if (gin.geomconstraints.found && gin.geomconstraints.ntcph == 0 &&
                  gin.geomconstraints.ntcpn == 0 && gin.system.found &&
                  gin.system.nsm == 0 && gin.initialise.ntishk != 0)
            printError("You have turned SHAKE off and have no solvent NTISHK in INITIALISE block should also be 0");
          if (gin.thermostat.found && gin.thermostat.ntt != 3 &&
                  gin.initialise.ntinht != 0)
            printError("If NTT!=3 in THERMOSTAT block, NTINHT in INITIALISE block should be 0");
          if (gin.multibath.found && gin.multibath.algorithm <= 1 &&
                  gin.initialise.ntinht != 0)
            printError("You want to initialise the Nose-Hoover variables (NTINHT=1 in INITIALISE block), but you are not doing Nose-Hoover-chains temperature coupling");
          if (gin.barostat.found && gin.barostat.ntp != 3 &&
                  gin.initialise.ntinhb != 0)
            printError("NTP!=3 in BAROSTAT block required NTINHB!=0 in INITIALISE block");

          if ((!gin.stochdyn.found ||
                  (gin.stochdyn.found && gin.stochdyn.ntsd == 0)) &&
                  gin.initialise.ntisti != 0)
            printError("You are not doing SD, but want to generate Stochastic integrals (NTISTI in INITIALISE block)");
          bool velblock = false;
          bool nhtvariable = false;
          bool nhpvariable = false;
          bool latticeshift = false;
          bool rototransref = false;
          bool stochint = false;
          bool lehistory = false;
          bool pertdata = false;
          bool box = false;
          // only do these checks for the first script!
          if (!first_script) {
            velblock = true;
            nhtvariable = true;
            nhpvariable = true;
            latticeshift = true;
            rototransref = true;
            stochint = true;
            lehistory = true;
            pertdata = true;
            box = true;
          }
          for (unsigned int i = 0; i < crd.blocks.size(); i++) {
            if ((crd.blocks[i] == "VELOCITY" || crd.blocks[i] == "VELOCITYRED") &&
                    crd.blockslength[i] != 0) velblock = true;
            if (crd.blocks[i] == "NHTVARIABLES" && crd.blockslength[i] != 0)
              nhtvariable = true;
            if (crd.blocks[i] == "NHPVARIABLES" && crd.blockslength[i] != 0)
              nhpvariable = true;
            if (crd.blocks[i] == "LATTICESHIFTS" && crd.blockslength[i] != 0)
              latticeshift = true;
            if (crd.blocks[i] == "ROTOTRANSREF" && crd.blockslength[i] != 0)
              rototransref = true;
            if (crd.blocks[i] == "STOCHINT" && crd.blockslength[i] != 0)
              stochint = true;
            if (crd.blocks[i] == "LEHISTORY" && crd.blockslength[i] != 0)
              lehistory = true;
            if (crd.blocks[i] == "PERTDATA" && crd.blockslength[i] != 0)
              pertdata = true;
            if (crd.blocks[i] == "GENBOX" && crd.blockslength[i] != 0)
              box = true;
          }
          if ((!gin.energymin.found ||
                  (gin.energymin.found && gin.energymin.ntem == 0)) &&
                  gin.initialise.ntivel == 0 && !velblock)
            //printError("You want to read velocities from file, but no VELOCITY(RED) block available");
          if ((gin.thermostat.found && gin.thermostat.ntt == 3) &&
                  gin.initialise.ntinht == 0 && !nhtvariable)
            printError("You want to read the Nose-Hoover thermostat variables from file, but no NHTVARIABLE block in coordinate file");
          if ((gin.multibath.found && gin.multibath.algorithm > 0) &&
                  gin.initialise.ntinht == 0 && !nhtvariable)
            printError("You want to read the Nose-Hoover thermostat variables from file, but not NHTVARIABLE block in coordinate file");
          if (gin.barostat.found && gin.barostat.ntp == 3 &&
                  gin.initialise.ntinhb == 0 && !nhpvariable)
            printError("You want to read the Nose-Hoover barostat variables from file, but no NHPVARIABLE block in coordinate file");
          if (gin.boundcond.found && gin.boundcond.ntb != 0 &&
                  gin.initialise.ntishi == 0 && !latticeshift)
            printError("You want to read the lattice shift vectors from file, but no LATTICESHIFT block in coordinate file");
          if (gin.overalltransrot.found && gin.overalltransrot.ncmro != 0 &&
                  gin.initialise.ntirtc == 0 && !rototransref)
            printError("You want to read the initial positions and orientations from file but no ROTOTRANSREF block in coordinate file");
          if (gin.stochdyn.found && gin.stochdyn.ntsd != 0 &&
                  gin.initialise.ntisti == 0 && !stochint)
            printError("You want to read the stocahstic integrals from file, but no STOCHINT block in coordinate file");
          if (gin.localelev.found && gin.localelev.ntles != 0 &&
                  gin.localelev.ntlesa == 1 && !lehistory)
            printError("You want to read the average local-elevation data from file, but no STOCHINT block in coordinate file");
          if (gin.perturbation.found && gin.perturbation.ntg != 0 &&
                  gin.perturbation.nrdgl == 1 && !pertdata)
            printError("You want to read the lambda value from file, but not PERTDATA block in coordinate file");
          if (gin.boundcond.found && gin.boundcond.ntb != 0 && !box)
            printError("Non-vacuum boundary conditions, but no GENBOX block in coordinate file");
        }
        //
	// INNERLOOP
	if (gin.innerloop.found) {
	  if (gin.innerloop.ntilm < 0 || gin.innerloop.ntilm > 4) {
	    std::stringstream ss;
	    ss << gin.innerloop.ntilm;
	    printIO("INNERLOOP", "NTILM", ss.str(), "0,1,2,3,4");
	  }
	  if (gin.innerloop.ntils != 0 && gin.innerloop.ntils != 1) {
	    std::stringstream ss;
	    ss << gin.innerloop.ntils;
	    printIO("INNERLOOP", "NTILS", ss.str(), "0,1");
	  }
	}
	//
	// INTEGRATE
	if (gin.integrate.found) {
	  if (gin.integrate.nint != 0 && gin.integrate.nint != 1) {
	    std::stringstream ss;
	    ss << gin.integrate.nint;
	    printIO("INTEGRATE", "NINT", ss.str(), "0,1");
	  }
	}
	//
	// JVALUES
	if (gin.jvalueres.found) {
	  if (gin.jvalueres.ntjvr < -3 || gin.jvalueres.ntjvr > 2) {
	    std::stringstream ss;
	    ss << gin.jvalueres.ntjvr;
	    printIO("JVALUERES", "NTJVR", ss.str(), "-3,-2,-1,0,1,2");
	  }
	  if (gin.jvalueres.ntjvra != 0 && gin.jvalueres.ntjvra != 1) {
	    std::stringstream ss;
	    ss << gin.jvalueres.ntjvra;
	    printIO("JVALUERES", "NTJVRA", ss.str(), "0,1");
	  }
	  if (gin.jvalueres.cjvr < 0.0) {
	    std::stringstream ss;
	    ss << gin.jvalueres.cjvr;
	    printIO("JVALUERES", "CJVR", ss.str(), ">= 0.0");
	  }
	  if (gin.jvalueres.taujvr < 0.0) {
	    std::stringstream ss;
	    ss << gin.jvalueres.taujvr;
	    printIO("JVALUERES", "TAUJVR", ss.str(), ">= 0.0");
	  }
	  if (gin.jvalueres.le != 0 && gin.jvalueres.le != 1) {
	    std::stringstream ss;
	    ss << gin.jvalueres.le;
	    printIO("JVALUERES", "LE", ss.str(), "0,1");
	  }
	  if (gin.jvalueres.ngrid <= 0) {
	    std::stringstream ss;
	    ss << gin.jvalueres.ngrid;
	    printIO("JVALUERES", "NGRID", ss.str(), "> 0");
	  }
	  if (gin.jvalueres.delta < 0.0) {
	    std::stringstream ss;
	    ss << gin.jvalueres.delta;
	    printIO("JVALUERES", "DELTA", ss.str(), ">= 0.0");
	  }
	  if (gin.jvalueres.write < 0) {
	    std::stringstream ss;
	    ss << gin.jvalueres.write;
	    printIO("JVALUERES", "NTWJV", ss.str(), ">= 0");
	  }
	}
        //
        // LAMBDAS block
        if (gin.lambdas.found && gin.lambdas.ntil == 1) {
          if (!gin.perturbation.found ||
                  (gin.perturbation.found && gin.perturbation.ntg == 0))
            printWarning("Specification of LAMBDAS block with NTIL!=0 is only effective if NTG!=0 in PERTURBATION block");
          int maxnlg = 0;
          if (gin.force.found) maxnlg = gin.force.nre.size();
          for (unsigned int i = 0; i < gin.lambdas.lambints.size(); i++) {
            if (gin.lambdas.lambints[i].nilg1 > gin.lambdas.lambints[i].nilg2) {
              std::stringstream ss;
              ss << "NILG1 = " << gin.lambdas.lambints[i].nilg1 << "\nNILG2 = "
                      << gin.lambdas.lambints[i].nilg2
                      << "\nNILG1 > NILG2 in LAMBDAS block is not allowed";
              printError(ss.str());
            }

            if (gin.lambdas.lambints[i].nilg1 > maxnlg) {
              std::stringstream ss;
              ss << "NILG1 = " << gin.lambdas.lambints[i].nilg1
                      << " in LAMBDAS block is larger than NEGR ("
                      << maxnlg << ") in FORCE block";
              printError(ss.str());
            }
            if (gin.lambdas.lambints[i].nilg2 > maxnlg) {
              std::stringstream ss;
              ss << "NILG2 = " << gin.lambdas.lambints[i].nilg2
                      << " in LAMBDAS block is larger than NEGR ("
                      << maxnlg << ") in FORCE block";
              printError(ss.str());
            }
          }
        }
        //
        // LOCALELEV block
	if (gin.localelev.found) {
	  if (gin.localelev.ntles < 0 || gin.localelev.ntles > 1) {
	    std::stringstream ss;
	    ss << gin.localelev.ntles;
	    printIO("LOCALELEV", "NTLES", ss.str(), "0, 1");
	  }
	  if (gin.localelev.ntlesa < 0 || gin.localelev.ntlesa > 2) {
	    std::stringstream ss;
	    ss << gin.localelev.ntlesa;
	    printIO("LOCALELEV", "NTLESA", ss.str(), "0..2");
	  }
	  for (map<int, int>::iterator it = gin.localelev.nlepid_ntlerf.begin();
		  it != gin.localelev.nlepid_ntlerf.end(); ++it) {
	    if (it->second < 0 || it->second > 1) {
	      std::stringstream ss, is;
	      is << "NTLEF[" << it->first << "]";
	      ss << it->second;
	      printIO("LOCALELEV", is.str(), ss.str(), "0, 1");
	    }
	  }
	}
        if (gin.localelev.found && gin.localelev.ntles != 0) {
          if (gin.energymin.ntem != 0) {
            stringstream msg;
            msg << "NTEM != 0 in ENERGYMIN block requires NTLES = 0 in LOCALELEV block"
                    << ", read: NTLES = " << gin.localelev.ntles;
            printError("NTEM != 0 in ENERGYMIN block requires NTLES = 0 in LOCALELEV block");
          }
        }

        /*
         * Andreas:
         * these are the "old" test, speak with Halvor and implemnet the new one!

        if (gin.localelev.found && gin.localelev.ntles != 0) {
          if (!l_ledih)
            printError("Local elevation simulation, but no localelevspec file given");
          else {
            fileInfo le;
            Ginstream ile(s_ledih);
            ile >> le;
            bool localelevspec = false;
            for (unsigned int i = 0; i < le.blocks.size(); i++) {
              if (le.blocks[i] == "LOCALELEVSPEC" && le.blockslength[i] != 0)
                localelevspec = true;
            }
            if (!localelevspec) {
              std::stringstream ss;
              ss << "Local elevation simulation, but file " << s_ledih
                      << "does not contains LOCALELEVSPEC block";
              printError(ss.str());
            }
          }
          if (gin.localelev.ntlesa == 2)
            printError("mk_script does not know how to handle the 'IOLEUS' file, please report this error");
          if (gin.localelev.ntles == 1 && gin.localelev.ntlesa == 0)
            printWarning("Zero initial averages in LOCALELEV block is not wise if a continuation run is being done");
          if (gin.localelev.ntlefr == 1 && gin.localelev.ntlesa == 0)
            printError("if NTLEFR=1 in LOCALELEV block, NTLESA cannot be 0");
          if (gin.localelev.ntles == 2 && gin.localelev.ntlesa != 2)
            printError("if NTLES=2 in LOCALELEV block, NTLESA should be 2");
          if (gin.localelev.ntles == 2 && gin.localelev.ntlefr != 1)
            printError("if NTLES=2 in LOCALELEV block, NTLEFR should be 1");
          if (gin.localelev.ntlesa == 2 && gin.localelev.ntles != 2)
            printError("if NTLESA=2 in LOCALELEV block, NTLES should be 2");
          if (gin.localelev.ntlefu == 0 && gin.localelev.wles != 1.0)
            printError("if NTLEFU=0 in LOCALELEV block, WLES should be 1.0");
          if (gin.localelev.ntlefu == 0 && gin.localelev.rles != 1.0)
            printError("if NTLEFU=0 in LOCALELEV block, RLES should be 1.0");
        }*/
        //
        // MULTIBATH block
	if (gin.multibath.found) {
	  if (gin.multibath.algorithm < 0 || gin.multibath.algorithm > 2) {
	    stringstream ss;
	    ss << gin.multibath.algorithm;
	    printIO("MULTIBATH", "algorithm", ss.str(), "1..3");
	  }
	  for (int i = 0; i < gin.multibath.nbaths; ++i) {
	    if (gin.multibath.temp0[i] < 0.0) {
	      std::stringstream ss, si;
	      si << "TEMP0[" << i + 1 << "]";
	      ss << gin.multibath.temp0[i];
	      printIO("MULTIBATH", si.str(), ss.str(), ">= 0.0");
	    }
	    if (gin.multibath.tau[i] < 0.0) {
	      std::stringstream ss, si;
	      si << "TAU[" << i + 1 << "]";
	      ss << gin.multibath.tau[i];
	      printIO("MULTIBATH", si.str(), ss.str(), ">= 0.0, -1");
	    }
	  }
	  for (int i = 0; i < gin.multibath.dofset; ++i) {
	    if (gin.multibath.last[i] < 0) {
	      std::stringstream ss, si;
	      si << "LAST[" << i + 1 << "]";
	      ss << gin.multibath.last[i];
	      printIO("MULTIBATH", si.str(), ss.str(), ">= 0");
	    }
	    if (gin.multibath.combath[i] < 1) {
	      std::stringstream ss, si;
	      si << "COMBATH[" << i + 1 << "]";
	      ss << gin.multibath.combath[i];
	      printIO("MULTIBATH", si.str(), ss.str(), ">= 1");
	    }
	    if (gin.multibath.irbath[i] < 1) {
	      std::stringstream ss, si;
	      si << "IRBATH[" << i + 1 << "]";
	      ss << gin.multibath.irbath[i];
	      printIO("MULTIBATH", si.str(), ss.str(), ">= 1");
	    }
	  }
          int mxlast = 0;
          for (unsigned int i = 0; i < gin.multibath.last.size(); i++) {
            if (gin.multibath.last[i] > mxlast) mxlast = gin.multibath.last[i];
          }
          if (mxlast != numTotalAtoms) {
            std::stringstream ss;
            ss << "Highest occuring LAST atom in MULTIBATH ("
                    << mxlast << ") should be equal to the total number of atoms ("
                    << numTotalAtoms << ")";
            printWarning(ss.str());
          }
          if (gin.pressurescale.found && gin.pressurescale.couple == 2) {
            for (unsigned int i = 0; i < gin.multibath.tau.size(); i++) {
              if (gin.pressurescale.taup <= gin.multibath.tau[i]) {
                std::stringstream ss;
                ss << "TAUP=" << gin.pressurescale.taup
                        << " in PRESSURESCALE block is larger than TAU[" << i + 1
                        << "]=" << gin.multibath.tau[i] << "in MULTIBATH block";
                printWarning(ss.str());
              }
            }
          }
          // THERE ARE PROBABLY MORE TESTS TO BE DONE
        }
        //
        // MULTICELL block
	if (gin.multicell.found) {
	  if (gin.multicell.ntm < 0 || gin.multicell.ntm > 1) {
	    std::stringstream ss;
	    ss << gin.multicell.ntm;
	    printIO("MULTICELL", "NTM", ss.str(), "0,1");
	  }
	  if (gin.multicell.ncella < 1) {
	    std::stringstream ss;
	    ss << gin.multicell.ncella;
	    printIO("MULTICELL", "NCELLA", ss.str(), ">= 1");
	  }
	  if (gin.multicell.ncellb < 1) {
	    std::stringstream ss;
	    ss << gin.multicell.ncellb;
	    printIO("MULTICELL", "NCELLB", ss.str(), ">= 1");
	  }
	  if (gin.multicell.ncellc < 1) {
	    std::stringstream ss;
	    ss << gin.multicell.ncellc;
	    printIO("MULTICELL", "NCELLC", ss.str(), ">= 1");
	  }
          if (!gromosXX) {
	    if (gin.multicell.tolpx <= 0.0) {
	      std::stringstream ss;
	      ss << gin.multicell.tolpx;
	      printIO("MULTICELL", "TOLPX", ss.str(), "> 0.0");
	    }
	    if (gin.multicell.tolpv <= 0.0) {
	      std::stringstream ss;
	      ss << gin.multicell.tolpv;
	      printIO("MULTICELL", "TOLPV", ss.str(), "> 0.0");
	    }
	    if (gin.multicell.tolpf <= 0.0) {
	      std::stringstream ss;
	      ss << gin.multicell.tolpf;
	      printIO("MULTICELL", "TOLPF", ss.str(), "> 0.0");
	    }
	    if (gin.multicell.tolpfw <= 0.0) {
	      std::stringstream ss;
	      ss << gin.multicell.tolpfw;
	      printIO("MULTICELL", "TOLPFW", ss.str(), "> 0.0");
  	    }
            if (gin.boundcond.found &&
                    (gin.boundcond.ntb < 1 || gin.boundcond.ntb > 2) &&
                    gin.multicell.ntm != 0) {
              printError("NTM!=0 in MULTICELL block is only allowed for NTB=1,2 in BOUNDCOND block");
            }
          }
        }
        //
        //MULTISTEP
        if(gin.multistep.found) {
	  if (gin.multistep.steps < 0) {
	    std::stringstream ss;
	    ss << gin.multistep.steps;
	    printIO("MULTISTEP", "STEPS", ss.str(), ">= 0");
	  }
	  if (gin.multistep.boost != 0 && gin.multistep.boost != 1) {
	    std::stringstream ss;
	    ss << gin.multistep.boost;
	    printIO("MULTISTEP", "BOOST", ss.str(), "0, 1");
	  }
          if(gin.multistep.steps < 0) {
            printError("STEPS in MULTISTEP block must not be < 0");
          }
          if(gin.multistep.boost != 0 && gin.multistep.boost != 1) {
            printError("BOOST in MULTISTEP block must be 0 or 1");
          }
        }
        //
        // NEIGHBOURLIST block
        if (gin.neighbourlist.found) {
	  if (gin.neighbourlist.nmprpl < 0 || gin.neighbourlist.nmprpl > 3) {
	    std::stringstream ss;
	    ss << gin.neighbourlist.nmprpl;
	    printIO("NEIGBOURLIST", "NMPRPL", ss.str(), "0,1,2,3");
	  }
	  if (gin.neighbourlist.nuprpl < 0) {
	    std::stringstream ss;
	    ss << gin.neighbourlist.nuprpl;
	    printIO("NEIGBOURLIST", "NMPRPL", ss.str(), ">= 0");
	  }
	  if (gin.neighbourlist.rcprpl < 0) {
	    std::stringstream ss;
	    ss << gin.neighbourlist.rcprpl;
	    printIO("NEIGBOURLIST", "RCPRPL", ss.str(), ">= 0.0");
	  }
	  if (gin.neighbourlist.grprpl < 0 || gin.neighbourlist.grprpl > 1.0) {
	    std::stringstream ss;
	    ss << gin.neighbourlist.grprpl;
	    printIO("NEIGBOURLIST", "GRPRPL", ss.str(), ">= 0.0 and <= 1.0");
	  }
	  if (gin.neighbourlist.nmtwpl < 0 || gin.neighbourlist.nmtwpl > 3) {
	    std::stringstream ss;
	    ss << gin.neighbourlist.nmtwpl;
	    printIO("NEIGBOURLIST", "NMTWPL", ss.str(), "0,1,2,3");
	  }
	  if (gin.neighbourlist.nutwpl < 0) {
	    std::stringstream ss;
	    ss << gin.neighbourlist.nutwpl;
	    printIO("NEIGBOURLIST", "NUTWPL", ss.str(), ">= 0");
	  }
	  if (gin.neighbourlist.rstwpl < 0) {
	    std::stringstream ss;
	    ss << gin.neighbourlist.rstwpl;
	    printIO("NEIGBOURLIST", "RSTWPL", ss.str(), ">= 0");
	  }
	  if (gin.neighbourlist.rltwpl < gin.neighbourlist.rstwpl) {
	    std::stringstream ss, si;
	    ss << gin.neighbourlist.nmprpl;
	    si << ">= RSTWPL (" << gin.neighbourlist.rstwpl << ")";
	    printIO("NEIGBOURLIST", "RLTWPL", ss.str(), si.str());
	  }
	  if (gin.neighbourlist.nuirin < 0) {
	    std::stringstream ss;
	    ss << gin.neighbourlist.nuirin;
	    printIO("NEIGBOURLIST", "NUIRIN", ss.str(), ">= 0");
	  }
	  if (gin.neighbourlist.nusrin < 0) {
	    std::stringstream ss;
	    ss << gin.neighbourlist.nusrin;
	    printIO("NEIGBOURLIST", "NMPRPL", ss.str(), ">= 0");
	  }
	  if (gin.neighbourlist.nmtwin < 0 || gin.neighbourlist.nmtwin > 2) {
	    std::stringstream ss;
	    ss << gin.neighbourlist.nmtwin;
	    printIO("NEIGBOURLIST", "NMTWIN", ss.str(), "0,1,2");
	  }
	  if (gin.neighbourlist.rctwin < 0.0) {
	    std::stringstream ss;
	    ss << gin.neighbourlist.rctwin;
	    printIO("NEIGBOURLIST", "RCTWIN", ss.str(), ">= 0.0");
	  }
	  if (gin.neighbourlist.ncgcen < -2) {
	    std::stringstream ss;
	    ss << gin.neighbourlist.nmprpl;
	    printIO("NEIGBOURLIST", "NCGCEN", ss.str(), ">= -2");
	  }
          if (gin.boundcond.found && gin.boundcond.ntb == 0 &&
                  gin.neighbourlist.nmprpl != 0)
            printError("NMPRPL!=0 in NEIGHBOURLIST block cannot be done under vacuum boundary conditions (NTB=0 in BOUNDCOND block");
          if (gin.neighbourlist.nmprpl == 0 && gin.neighbourlist.nuprpl != 0)
            printError("NMPRPL=0 in NEIGHBOURLIST requires NUPRPL=0");
          if (gin.neighbourlist.nmprpl == 0 && gin.neighbourlist.rcprpl != 0.0)
            printError("NMPRPL=0 in NEIGHBOURLIST requires RCPRPL=0.0");
          if (gin.neighbourlist.nmprpl == 0 && gin.neighbourlist.grprpl != 0.0)
            printError("NMPRPL=0 in NEIGHBOURLIST requires GRPRPL=0.0");
          if (gin.neighbourlist.nmprpl != 0 && gin.neighbourlist.nuprpl == 0)
            printError("NMPRPL!=0 in NEIGHBOURLIST requires NUPRPL!=0");
          if (gin.neighbourlist.nmprpl != 0 && gin.neighbourlist.rcprpl == 0.0)
            printError("NMPRPL!=0 in NEIGHBOURLIST requires RCPRPL!=0.0");
          if (gin.neighbourlist.nmprpl != 0 && gin.neighbourlist.grprpl == 0.0)
            printError("NMPRPL!=0 in NEIGHBOURLIST requires GRPRPL!=0.0");
          if (gin.neighbourlist.nmprpl != 0 &&
                  gin.neighbourlist.rcprpl < gin.neighbourlist.rltwpl)
            printError("NMPRPL!=0 in NEIGHBOURLIST requires RCPRPL>=RLTWPL");
          if (gin.neighbourlist.nmprpl == 1 && gin.neighbourlist.nmtwpl == 3)
            printError("NMPRPL=1 in NEIGHBOURLIST requires NMTWPL!=3");
          if ((gin.neighbourlist.nmprpl == 2 || gin.neighbourlist.nmprpl == 3) &&
                  !(gin.neighbourlist.nmtwpl == 0 || gin.neighbourlist.nmtwpl == 3))
            printError("NMPRPL=2 or 3 in NEIGHBOURLIST block required NMTWPL=0 or 3");
          if (gin.neighbourlist.nmtwpl == 0 && gin.neighbourlist.nutwpl != 0)
            printError("NMTWPL=0 in NEIGHBOURLIST block requires NUTWPL=0");
          if (gin.neighbourlist.nmtwpl == 0 && gin.neighbourlist.rstwpl != 0.0)
            printError("NMTWPL=0 in NEIGHBOURLIST block requires RSTWPL=0.0");
          if (gin.neighbourlist.nmtwpl == 0 && gin.neighbourlist.rltwpl != 0.0)
            printError("NMTWPL=0 in NEIGHBOURLIST block requires RLTWPL=0.0");
          if (gin.neighbourlist.nmtwpl != 0 && gin.neighbourlist.nutwpl == 0)
            printError("NMTWPL!=0 in NEIGHBOURLIST block requires NUTWPL!=0");
          if (gin.neighbourlist.nmtwpl != 0 && gin.neighbourlist.rstwpl == 0.0)
            printError("NMTWPL!=0 in NEIGHBOURLIST block requires RSTWPL!=0.0");
          if (gin.neighbourlist.nmtwpl != 0 && gin.neighbourlist.rltwpl == 0.0)
            printError("NMTWPL!=0 in NEIGHBOURLIST block requires RLTWPL!=0.0");
          if (gin.neighbourlist.rltwpl < gin.neighbourlist.rstwpl)
            printError("RLTWPL should be larger than RSTWPL in NEIGHBOURLIST block");
          if (gin.neighbourlist.nuirin != 0 && gin.neighbourlist.nmtwpl == 0)
            printError("NUIRIN!=0 in NEIGHBOURLIST block requires NMTWPL!=0");
          if (gin.neighbourlist.nuirin != 0 &&
                  gin.neighbourlist.rltwpl <= gin.neighbourlist.rstwpl)
            printError("NUIRIN!=0 in NEIGHBOURLIST block requires RLTWPL> RSTWPL");
          if (gin.neighbourlist.nuirin != 0 && gin.neighbourlist.nusrin == 0)
            printError("NUIRIN!=0 in NEIGHBOURLIST block requires NUSRIN!=0");
          if (gin.neighbourlist.nusrin != 0 && gin.neighbourlist.nmtwpl == 0)
            printError("NUSRIN!=0 in NEIGHBOURLIST block requires NMTWPL!=0");
          if (gin.neighbourlist.nmtwin == 0 && gin.neighbourlist.rctwin != 0.0)
            printError("NMTWIN=0 in NEIGHBOURLIST block requires RCTWIN=0.0");
          if (gin.neighbourlist.nmtwin != 0 && gin.neighbourlist.rctwin == 0.0)
            printError("NMTWIN!=0 in NEIGHBOURLIST block requires RCTWIN!=0.0");
          if (gin.neighbourlist.nmtwin != 0 &&
                  gin.neighbourlist.rctwin > gin.neighbourlist.rltwpl)
            printError("NMTWIN!=0 in NEIGHBOURLIST block requires RCTWIN <= RLTWPL");
          if (gin.neighbourlist.nmtwin == 1 && gin.neighbourlist.nmtwpl != 1)
            printError("NMTWIN=1 in NEIGHBOURLIST block requires NMTWPL=1");
          if (gin.neighbourlist.nmtwin == 2 &&
                  !(gin.neighbourlist.nmtwpl == 2 || gin.neighbourlist.nmtwpl == 3))
            printError("NMTWIN=2 in NEIGHBOURLIST block requires NMTWPL=2 or 3");
          if (gin.gromos96compat.found && gin.gromos96compat.ntnb96 != 0) {
            if (!gromosXX) {
              if (gin.neighbourlist.nmprpl != 0)
                printError("NTNB96!=0 in GROMOS96COMPAT block requires NMPRPL=0 in NEIGHBOURLIST block");
              if (gin.neighbourlist.nmtwpl != 0 && gin.neighbourlist.nmtwpl != 1)
                printError("NTNB96!=0 in GROMOS96COMPAT block requires NMTWPL=0 or 1 in NEIGHBOURLIST block");
              if (gin.neighbourlist.nmtwin != 0)
                printError("NTNB96!=0 in GROMOS96COMPAT block requires NMTWIN=0 in NEIGHBOURLIST block");
            }
          }
          if (gin.gromos96compat.found && gin.gromos96compat.ntnb96 == 0) {
            if (gin.neighbourlist.nmtwpl != 0 && gin.neighbourlist.nmtwpl != 2)
              printError("NTNB96=0 in GROMOS96COMPAT block requires NMTWPL=0 or 2 in NEIGHBOURLIST block");
          }
          if (gin.neighbourlist.nmprpl > 1)
            printError("NMPRPL > 1 in NEIGHBOURLIST block not yet implemented");
          if (gin.neighbourlist.nuprpl < 0)
            printError("NUPRPL < 0 in NEIGHBOURLIST block not yet implemented");
          if (gin.neighbourlist.nmtwpl == 3)
            printError("NMTWPL=3 in NEIGHBOURLIST block not yet implemented");
          if (gin.neighbourlist.nmtwin == 1)
            printError("NMTWIN=1 in NEIGHBOURLIST block not yet implemented");
          if (gin.force.found &&
                  (gin.force.ntf[8] != 0 || gin.force.ntf[9] != 0) &&
                  gin.neighbourlist.rltwpl > gin.neighbourlist.rstwpl &&
                  gin.neighbourlist.nuirin == 0)
            printError("NTF[9]!=0 or NTF[10]!=0 with RLTWPL > RSTWPL requires NUIRIN!=0 in NEIGHBOURLIST block");
          if (gin.force.found &&
                  (gin.force.ntf[8] != 0 || gin.force.ntf[9] != 0) &&
                  gin.neighbourlist.nusrin == 0)
            printError("NTF[9]!=0 or NTF[10]!=0 requires NUSRIN!=0 in NEIGHBOURLIST block");
          if (gin.neighbourlist.nmprpl != 0 && gin.neighbourlist.nmtwpl == 0)
            printWarning("NMPRPL!=0 and NMTWPL=0 in NEIGHBOURLIST block, this means that you make a PR-pairlist, but don't use it");
          if (gin.neighbourlist.nmtwpl != 0 &&
                  gin.neighbourlist.rltwpl > gin.neighbourlist.rstwpl &&
                  gin.neighbourlist.nuirin == 0)
            printWarning("NMTWPL!=0, RLTWPL > RSTWPL and NUIRIN=0 in NEIGHBOURLIST block, this means that you make an IR-pairlist, but don't use it");
          if (gin.neighbourlist.nmtwpl != 0 &&
                  gin.neighbourlist.nusrin == 0)
            printWarning("NMTWPL!=0 and NUSRIN=0 in NEIGHBOURLIST block, this means that you make an SR-pairlist, but don't use it");
          if (gin.force.found &&
                  gin.force.ntf[8] == 0 && gin.force.ntf[9] == 0 &&
                  gin.neighbourlist.nuirin == 0)
            printWarning("NUIRIN=0 in NEIGHBOURLIST blocl, but NTF[9]=NTF[10] in FORCE block, this means that you calculate the IR interactions, but they are 0");
          if (gin.force.found &&
                  gin.force.ntf[8] == 0 && gin.force.ntf[9] == 0 &&
                  gin.neighbourlist.nusrin == 0)
            printWarning("NUSRIN=0 in NEIGHBOURLIST blocl, but NTF[9]=NTF[10] in FORCE block, this means that you calculate the SR interactions, but they are 0");
          if (gin.neighbourlist.nmprpl != 0 &&
                  gin.neighbourlist.nmtwpl != 0 &&
                  gin.neighbourlist.nuprpl != gin.neighbourlist.nutwpl)
            printWarning("NUPRPL!=NUTWPL in NEIGHBOURLIST block, this leads to asynchronous updates of the PR-pairlist and the TW pairlist.");
          if (gin.neighbourlist.nmtwpl != 0 &&
                  gin.neighbourlist.nuirin != 0 &&
                  gin.neighbourlist.nutwpl != gin.neighbourlist.nuirin)
            printWarning("NUTWPL!=NUIRIN in NEIGHBOURLIST block, this leads to asynchronous updates of the TW-pairlist and the IR interactions.");
          if (gin.neighbourlist.nusrin != 0 && gin.neighbourlist.nusrin != 1)
            printWarning("NUSRIN!=1 means that you do not update the short-range interactions every step. Are you sure?");
          if (gin.readtraj.found && gin.readtraj.ntrd != 0 &&
                  gin.neighbourlist.nutwpl != 0 && gin.neighbourlist.nutwpl != 1)
            printWarning("When reading a trajectory, NUTWPL!=1 in NEIGHBOURLIST is not wise");
          if (gin.neighbourlist.ncgcen != 0)
            printWarning("Using NCGCEN!=0 is not according to the GROMOS convention");
        }
        //
        // NONBONDED block	
        if (gin.nonbonded.found) {
	  if (gin.nonbonded.nlrele < -4 || gin.nonbonded.nlrele > 4) {
	    std::stringstream ss;
	    ss << gin.nonbonded.nlrele;
	    printIO("NONBONDED", "NLRELE", ss.str(), "-4,-3,-2,-1,0,1,2,3,4");
	  }
	  if (gin.nonbonded.appak < 0.0) {
	    std::stringstream ss;
	    ss << gin.nonbonded.appak;
	    printIO("NONBONDED", "APPAK", ss.str(), ">= 0.0");
	  }
	  if (gin.nonbonded.rcrf < 0.0) {
	    std::stringstream ss;
	    ss << gin.nonbonded.rcrf;
	    printIO("NONBONDED", "RCRF", ss.str(), ">= 0.0");
	  }
	  if (gin.nonbonded.epsrf < 1.0 && gin.nonbonded.epsrf != 0.0) {
	    std::stringstream ss;
	    ss << gin.nonbonded.epsrf;
	    printIO("NONBONDED", "EPSRF", ss.str(), "0.0 (off) or >= 1.0");
	  }
	  if (gin.nonbonded.nshape < -1 || gin.nonbonded.nshape > 10) {
	    std::stringstream ss;
	    ss << gin.nonbonded.nshape;
	    printIO("NONBONDED", "NSHAPE", ss.str(), "-1 .. 10");
	  }
	  if (gin.nonbonded.ashape <= 0.0) {
	    std::stringstream ss;
	    ss << gin.nonbonded.ashape;
	    printIO("NONBONDED", "ASHAPE", ss.str(), "> 0.0");
	  }
	  if (gin.nonbonded.na2clc < 0 || gin.nonbonded.na2clc > 4) {
	    std::stringstream ss;
	    ss << gin.nonbonded.na2clc;
	    printIO("NONBONDED", "NA2CLC", ss.str(), "0,1,2,3,4");
	  }
	  if (gin.nonbonded.tola2 <= 0.0) {
	    std::stringstream ss;
	    ss << gin.nonbonded.tola2;
	    printIO("NONBONDED", "TOLA2", ss.str(), "> 0.0");
	  }
	  if (gin.nonbonded.epsls < 1.0 && gin.nonbonded.epsls != 0.0) {
	    std::stringstream ss;
	    ss << gin.nonbonded.epsls;
	    printIO("NONBONDED", "EPSLS", ss.str(), "0.0 (off) or >= 1.0");
	  }
	  if (gin.nonbonded.nkx <= 0) {
	    std::stringstream ss;
	    ss << gin.nonbonded.nkx;
	    printIO("NONBONDED", "NKX", ss.str(), "> 0");
	  }
	  if (gin.nonbonded.nky <= 0) {
	    std::stringstream ss;
	    ss << gin.nonbonded.nky;
	    printIO("NONBONDED", "NKY", ss.str(), "> 0");
	  }
	  if (gin.nonbonded.nkz <= 0) {
	    std::stringstream ss;
	    ss << gin.nonbonded.nkz;
	    printIO("NONBONDED", "NKZ", ss.str(), "> 0");
	  }
	  if (gin.nonbonded.kcut <= 0.0) {
	    std::stringstream ss;
	    ss << gin.nonbonded.kcut;
	    printIO("NONBONDED", "KCUT", ss.str(), "> 0.0");
	  }
	  if (gin.nonbonded.ngx <= 0) {
	    std::stringstream ss;
	    ss << gin.nonbonded.ngx;
	    printIO("NONBONDED", "NGX", ss.str(), "> 0");
	  }
	  if (gin.nonbonded.ngy <= 0) {
	    std::stringstream ss;
	    ss << gin.nonbonded.ngy;
	    printIO("NONBONDED", "NGY", ss.str(), "> 0");
	  }
	  if (gin.nonbonded.ngz <= 0) {
	    std::stringstream ss;
	    ss << gin.nonbonded.ngz;
	    printIO("NONBONDED", "NGZ", ss.str(), "> 0");
	  }
	  if (gin.nonbonded.nasord < 1 || gin.nonbonded.nasord > 5) {
	    std::stringstream ss;
	    ss << gin.nonbonded.nasord;
	    printIO("NONBONDED", "NASORD", ss.str(), "1,2,3,4,5");
	  }
	  if (gin.nonbonded.nfdord < 0 || gin.nonbonded.nfdord > 5) {
	    std::stringstream ss;
	    ss << gin.nonbonded.nfdord;
	    printIO("NONBONDED", "NFDORD", ss.str(), "0,1,2,3,4,5");
	  }
	  if (gin.nonbonded.nalias <= 0) {
	    std::stringstream ss;
	    ss << gin.nonbonded.nalias;
	    printIO("NONBONDED", "NALIAS", ss.str(), "> 0");
	  }
	  if (gin.nonbonded.nqeval < 0) {
	    std::stringstream ss;
	    ss << gin.nonbonded.nqeval;
	    printIO("NONBONDED", "NQEVAL", ss.str(), ">= 0");
	  }
	  if (gin.nonbonded.faccur <= 0.0) {
	    std::stringstream ss;
	    ss << gin.nonbonded.faccur;
	    printIO("NONBONDED", "FACCUR", ss.str(), "> 0");
	  }
	  if (gin.nonbonded.nrdgrd < 0 || gin.nonbonded.nrdgrd > 1) {
	    std::stringstream ss;
	    ss << gin.nonbonded.nrdgrd;
	    printIO("NONBONDED", "NRDGRD", ss.str(), "0,1");
	  }
	  if (gin.nonbonded.nwrgrd < 0 || gin.nonbonded.nwrgrd > 1) {
	    std::stringstream ss;
	    ss << gin.nonbonded.nwrgrd;
	    printIO("NONBONDED", "NWRGRD", ss.str(), "0,1");
	  }
	  if (gin.nonbonded.nlrlj < 0 || gin.nonbonded.nlrlj > 1) {
	    std::stringstream ss;
	    ss << gin.nonbonded.nlrlj;
	    printIO("NONBONDED", "NLRLJ", ss.str(), "0,1");
	  }
	  if (gin.nonbonded.slvdns <= 0.0) {
	    std::stringstream ss;
	    ss << gin.nonbonded.slvdns;
	    printIO("NONBONDED", "SLVDNS", ss.str(), "> 0.0");
	  }
          if (gin.boundcond.found && gin.boundcond.ntb == 0 &&
                  abs(gin.nonbonded.nlrele) > 1)
            printError("ABS(NLRELE) > 1 in NONBONDED block cannot be combined with vacuum boundary conditions (NTB=0 in BOUNDCOND block");
          if (gin.nonbonded.na2clc == 1 && abs(gin.nonbonded.nlrele) != 2)
            printError("NA2CLC=1 in NONBONDED block, requires NLRELE=-2 or 2");
          if (gin.nonbonded.na2clc == 4 &&
                  !(abs(gin.nonbonded.nlrele) == 3 || abs(gin.nonbonded.nlrele) == 4 ||
                  (gin.consistencycheck.found && gin.consistencycheck.ntchk == 1 &&
                  gin.consistencycheck.ntcke == 1)))
            printError("NA2CLC=4 in NONBONDED block, requires abs(NLRELE)=3 or 4 or NTCHK=1 and NTCKE=1");
          if (gin.nonbonded.ngx % 2 != 0)
            printError("NGX in NONBONDED block should be even");
          if (gin.nonbonded.ngy % 2 != 0)
            printError("NGY in NONBONDED blockshould be even");
          if (gin.nonbonded.ngz % 2 != 0)
            printError("NGZ in NONBONDED blockshould be even");
          if (gin.nonbonded.nlrele == 4 || gin.nonbonded.nlrele == -4)
            printError("NLRELE=-4 or 4 in NONBONDED block is not yet implemented");
          if (gin.nonbonded.nlrele != 0 && gin.nonbonded.nlrele != 1 &&
                  gin.pairlist.found &&
                  gin.nonbonded.ashape > gin.pairlist.rcutp)
            printWarning("NLRELE!=0,1 and ASHAPE > RCUTP");
          if (gin.nonbonded.nlrele != 0 && gin.nonbonded.nlrele != 1 &&
                  gin.neighbourlist.found &&
                  gin.nonbonded.ashape > gin.neighbourlist.rstwpl)
            printWarning("NLRELE!=0,1 and ASHAPE > RSTWPL");
          if (gin.nonbonded.na2clc == 0 &&
                  abs(gin.nonbonded.nlrele) != 1 && gin.nonbonded.nlrele != 0)
            printWarning("NA2CLC=0 in NONBONDED block is not very wise.");
        }
        //
        // OVERALLTRANSROT block
        if (gin.overalltransrot.found) {
	  if (gin.overalltransrot.ncmtr < 0 || gin.overalltransrot.ncmtr > 1) {
	    std::stringstream ss;
	    ss << gin.overalltransrot.ncmtr;
	    printIO("OVERALLTRANSROT", "NCMTR", ss.str(), "0,1");
	  }
	  if (gin.overalltransrot.ncmro < 0 || gin.overalltransrot.ncmro > 1) {
	    std::stringstream ss;
	    ss << gin.overalltransrot.ncmro;
	    printIO("OVERALLTRANSROT", "NCMRO", ss.str(), "0,1");
	  }
          if (gin.boundcond.found && gin.boundcond.ntb != 0 &&
                  gin.overalltransrot.ncmro != 0)
            printError("NCMRO!=0 in OVERALLTRANSROT block is only allowed for vacuum boundary conditions (NTB=0 in BOUNDCOND block)");
          int ndfmin = 0;
          if (gin.overalltransrot.ncmtr != 0) ndfmin += 3;
          if (gin.overalltransrot.ncmro != 0) ndfmin += 3;
          if (gin.boundcond.found && gin.boundcond.ndfmin != ndfmin) {
            std::stringstream ss;
            ss << "NCMTR=" << gin.overalltransrot.ncmtr << " and NCMRO="
                    << gin.overalltransrot.ncmro << " in OVERALLTRANSROT block "
                    << "should lead to NDFMIN=" << ndfmin << " in BOUNDCOND block\n"
                    << "You gave " << gin.boundcond.ndfmin;
            printWarning(ss.str());
          }
        }
        if (!gromosXX && ((!gin.overalltransrot.found ||
                (gin.overalltransrot.found && gin.overalltransrot.ncmtr == 0)) &&
                (!gin.stochdyn.found ||
                (gin.stochdyn.found && gin.stochdyn.ntsd == 0)) &&
                (!gin.energymin.found ||
                (gin.energymin.found && gin.energymin.ntem == 0)) &&
                (!gin.readtraj.found ||
                (gin.readtraj.found && gin.readtraj.ntrd == 0)) &&
                (!gin.consistencycheck.found ||
                (gin.consistencycheck.found && gin.consistencycheck.ntchk == 0))))
          printWarning("Running a simulation with NCMTR=0 in OVERALLTRANSROT block may be unwise");
        //
        // PAIRLIST block
        if (gin.pairlist.found) {
	  if (gin.pairlist.algorithm < 0 || gin.pairlist.algorithm > 1) {
	    stringstream ss;
	    ss << gin.pairlist.algorithm;
	    printIO("PAIRLIST", "algorithm", ss.str(), "0..1");
	  }
	  if (gin.pairlist.size <= 0.0) {
	    stringstream ss;
	    ss << gin.pairlist.size;
	    printIO("PAIRLIST", "SIZE", ss.str(), "> 0.0");
	  }
	  if (gin.pairlist.nsnb <= 0) {
	    std::stringstream ss;
	    ss << gin.pairlist.nsnb;
	    printIO("PAIRLIST", "NSNB", ss.str(), "> 0");
	  }
	  if (gin.pairlist.rcutl <= 0.0) {
	    std::stringstream ss;
	    ss << gin.pairlist.rcutl;
	    printIO("PAIRLIST", "RCUTL", ss.str(), "> 0.0");
	  }
	  if (gin.pairlist.rcutp <= 0.0) {
	    std::stringstream ss;
	    ss << gin.pairlist.rcutp;
	    printIO("PAIRLIST", "RCUTP", ss.str(), "> 0.0");
	  }
	  if (gin.pairlist.type < 0 || gin.pairlist.type > 1) {
	    stringstream ss;
	    ss << gin.pairlist.type;
	    printIO("PAIRLIST", "TYPE", ss.str(), "0..1");
	  }
          if (gin.pairlist.rcutp > gin.pairlist.rcutl)
            printError("In PAIRLIST block, RCUTP > RCUTL is not allowed");
          if (gin.nonbonded.found && gin.nonbonded.rcrf != 0.0 &&
                  gin.nonbonded.rcrf != gin.pairlist.rcutl)
            printWarning("Usually we set RCRF in NONBONDED block equal to RCUTL in PAIRLIST block");
        }
        //
	// PATHINT block
	if (gin.pathint.found) {
	  if (gin.pathint.ntpi < 0 || gin.pathint.ntpi > 1) {
	    std::stringstream ss;
	    ss << gin.pathint.ntpi;
	    printIO("PATHINT", "NTPI", ss.str(), "0,1");
	  }
	}
	if (gin.pathint.found && gin.pathint.ntpi == 1) {
	  if (gin.perturbation.found && gin.perturbation.ntg != 0)
	    printError("NTPI=1 in PATHINT block cannot be done with NTG!=0 in PERTURBATION block");
	}
	//
	// PERSCALE block
	if (gin.perscale.found) {
	  if (gin.perscale.restype != 1) {
	    stringstream ss;
	    ss << gin.perscale.restype;
	    printIO("PERSCALE", "RESTYPE", ss.str(), "0..1");
	  }
	  if (gin.perscale.kdih < 0.0) {
	    std::stringstream ss;
	    ss << gin.perscale.kdih;
	    printIO("PERSCALE", "KDIH", ss.str(), ">= 0.0");
	  }
	  if (gin.perscale.kj < 0.0) {
	    std::stringstream ss;
	    ss << gin.perscale.kj;
	    printIO("PERSCALE", "KJ", ss.str(), ">= 0.0");
	  }
	  if (gin.perscale.t <= 0) {
	    std::stringstream ss;
	    ss << gin.perscale.t;
	    printIO("PERSCALE", "T", ss.str(), "> 0");
	  }
	  if (gin.perscale.diff < 0.0) {
	    std::stringstream ss;
	    ss << gin.perscale.diff;
	    printIO("PERSCALE", "DIFF", ss.str(), ">= 0.0");
	  }
	  if (gin.perscale.ratio <= 0.0) {
	    std::stringstream ss;
	    ss << gin.perscale.ratio;
	    printIO("PERSCALE", "RATIO", ss.str(), "> 0.0");
	  }
	  if (gin.perscale.read < 0 || gin.perscale.read > 1) {
	    std::stringstream ss;
	    ss << gin.perscale.read;
	    printIO("PERSCALE", "READ", ss.str(), "0,1");
	  }
          if (!gin.jvalueres.found ||
                  (gin.jvalueres.found && gin.jvalueres.ntjvr == 0))
            printError("You can only use the PERSCALE block when J-value restraining (NTJVR!=0 in JVALUERES block)");
        }
        //
        // PERTURBATION block
	if (gin.perturbation.found) {
	  if (gin.perturbation.ntg < 0 || gin.perturbation.ntg > 1) {
	    std::stringstream ss;
	    ss << gin.perturbation.ntg;
	    printIO("PERTURBATION", "NTG", ss.str(), "0,1");
	  }
	  if (gin.perturbation.nrdgl < 0 || gin.perturbation.nrdgl > 1) {
	    std::stringstream ss;
	    ss << gin.perturbation.nrdgl;
	    printIO("PERTURBATION", "NRDGL", ss.str(), "0,1");
	  }
	  if (gin.perturbation.rlam < 0.0 || gin.perturbation.rlam > 1.0) {
	    std::stringstream ss;
	    ss << gin.perturbation.rlam;
	    printIO("PERTURBATION", "RLAM", ss.str(), ">= 0.0 and <= 1.0");
	  }
	  if (gin.perturbation.dlamt < 0.0) {
	    std::stringstream ss;
	    ss << gin.perturbation.dlamt;
	    printIO("PERTURBATION", "DLAMT", ss.str(), ">= 0.0");
	  }
	  if (gin.perturbation.alphlj < 0.0) {
	    std::stringstream ss;
	    ss << gin.perturbation.alphlj;
	    printIO("PERTURBATION", "ALPHLJ", ss.str(), ">= 0.0");
	  }
	  if (gin.perturbation.alphc < 0.0) {
	    std::stringstream ss;
	    ss << gin.perturbation.alphc;
	    printIO("PERTURBATION", "ALPHC", ss.str(), ">= 0.0");
	  }
	  if (gin.perturbation.nlam <= 0) {
	    std::stringstream ss;
	    ss << gin.perturbation.nlam;
	    printIO("PERTURBATION", "NLAM", ss.str(), "> 0");
	  }
	  if (gin.perturbation.nscale < 0 || gin.perturbation.nscale > 2) {
	    std::stringstream ss;
	    ss << gin.perturbation.nscale;
	    printIO("PERTURBATION", "NSCALE", ss.str(), "0,1,2");
	  }
	}
        if (gin.perturbation.found && gin.perturbation.ntg != 0) {
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
              printWarning(os.str());
              fin.close();
            } else {
              ostringstream os;
              os << "NTG = " << gin.perturbation.ntg
                      << " in PERTURBATION block, but "
                      << "no perturbation topology specified\n";

              printError(os.str());
            }
          }
          double rlamfin;
          rlamfin = gin.perturbation.rlam +
                  gin.perturbation.dlamt * gin.step.dt * gin.step.nstlim;
          if (rlamfin > 1.0) {
            ostringstream os;
            os << "Using RLAM = " << gin.perturbation.rlam << " and DLAMT = "
                    << gin.perturbation.dlamt
                    << " in the PERTURB block and NSTLIM = "
                    << gin.step.nstlim << " in the STEP block\n";
            os << "will lead to a final lambda value of " << rlamfin << endl;
            printWarning(os.str());
          }
        }
        if ((!gin.perturbation.found ||
                (gin.perturbation.found && gin.perturbation.ntg == 0)) &&
                l_pttopo) {
          string s = "Perturbation topology specified, but no perturbation ";
          s += "according to the input file\n";
          printWarning(s);
        }
        //
        // POLARISE block
	if (gin.polarise.found) {
	  if (gin.polarise.cos != 0 && gin.polarise.cos != 1) {
	    std::stringstream ss;
	    ss << gin.polarise.cos;
	    printIO("POLARISE", "COS", ss.str(), "0,1");
	  }
	  if (gin.polarise.efield != 0 && gin.polarise.efield != 1) {
	    std::stringstream ss;
	    ss << gin.polarise.efield;
	    printIO("POLARISE", "EFIELD", ss.str(), "0,1");
	  }
	  if (gin.polarise.minfield <= 0.0) {
	    std::stringstream ss;
	    ss << gin.polarise.minfield;
	    printIO("POLARISE", "MINFIELD", ss.str(), "> 0.0");
	  }
	  if (gin.polarise.damp != 0 && gin.polarise.damp != 1) {
	    std::stringstream ss;
	    ss << gin.polarise.damp;
	    printIO("POLARISE", "DAMP", ss.str(), "0,1");
	  }
	  if (gin.polarise.write < 0) {
	    std::stringstream ss;
	    ss << gin.polarise.write;
	    printIO("POLARISE", "WRITE", ss.str(), "> 0");
	  }
	 	}
        if (gin.polarise.found && gin.polarise.cos == 1) {
          // add checks here
        }
        //
        // POSITIONRES block
	if (gin.positionres.found) {
	  if (gin.positionres.ntpor < 0 || gin.positionres.ntpor > 3) {
	    std::stringstream ss;
	    ss << gin.positionres.ntpor;
	    printIO("POSITIONRES", "NTPOR", ss.str(), "0,1,2,3");
	  }
	  if (gin.positionres.ntporb < 0 || gin.positionres.ntporb > 1) {
	    std::stringstream ss;
	    ss << gin.positionres.ntporb;
	    printIO("POSITIONRES", "NTPORB", ss.str(), "0,1");
	  }
	  if (gin.positionres.ntpors < 0 || gin.positionres.ntpors > 1) {
	    std::stringstream ss;
	    ss << gin.positionres.ntpors;
	    printIO("POSITIONRES", "NTPORS", ss.str(), "0,1");
	  }
	  if (gin.positionres.cpor < 0) {
	    std::stringstream ss;
	    ss << gin.positionres.cpor;
	    printIO("POSITIONRES", "CPOR", ss.str(), ">= 0.0");
	  }
	 	}
        if (gin.positionres.found && gin.positionres.ntpor != 0) {
          if (gin.positionres.ntpor == 2 && gin.positionres.ntporb != 1)
            printError("NTPOR=2 in POSITIONRES block requires that NTPORB=1");
          if (!l_posresspec)
            printError("Position restraining, but no file specified");
          else {
            Ginstream ipos(s_posresspec);
            fileInfo posres;
            ipos >> posres;
            ipos.close();
            bool posresspec = false;
            for (unsigned int i = 0; i < posres.blocks.size(); i++) {
              if (posres.blocks[i] == "POSRESSPEC" && posres.blockslength[i] != 0)
                posresspec = true;
            }
            if (!posresspec) {
              std::stringstream ss;
              ss << "No POSRESSPEC block in file " << s_posresspec;

              printError(ss.str());
            }
          }
        }
        //
        // PRESSURESCALE block
        if (gin.pressurescale.found) {

	  if (gin.pressurescale.couple < 0 || gin.pressurescale.couple > 2) {
	    stringstream ss;
	    ss << gin.pressurescale.couple;
	    printIO("PRESSURESCALE", "COUPLE", ss.str(), "0..2");
	  }
	  if (gin.pressurescale.scale < 0 || gin.pressurescale.scale > 4) {
	    stringstream ss;
	    ss << gin.pressurescale.scale;
	    printIO("PRESSURESCALE", "SCALE", ss.str(), "off(0),iso(1),aniso(2),full(3),semianiso(4)");
	  }
	  if (gin.pressurescale.comp <= 0.0) {
	    std::stringstream ss;
	    ss << gin.pressurescale.comp;
	    printIO("PRESSURESCALE", "COMP", ss.str(), "> 0.0");
	  }
	  if (gin.pressurescale.taup <= 0.0) {
	    std::stringstream ss;
	    ss << gin.pressurescale.taup;
	    printIO("PRESSURESCALE", "COMP", ss.str(), "> 0.0");
	  }
	  if (gin.pressurescale.virial < 0 || gin.pressurescale.virial > 2) {
	    stringstream ss;
	    ss << gin.pressurescale.virial;
	    printIO("PRESSURESCALE", "VIRIAL", ss.str(), "0..2");
	  }
          if (gin.pressurescale.x_semi < 0 || gin.pressurescale.x_semi > 2) {
	    stringstream ss;
	    ss << gin.pressurescale.x_semi;
	    printIO("PRESSURESCALE", "X_SEMI", ss.str(), "0..2");
	  }
          if (gin.pressurescale.y_semi < 0 || gin.pressurescale.y_semi > 2) {
	    stringstream ss;
	    ss << gin.pressurescale.y_semi;
	    printIO("PRESSURESCALE", "Y_SEMI", ss.str(), "0..2");
	  }
          if (gin.pressurescale.z_semi < 0 || gin.pressurescale.z_semi > 2) {
	    stringstream ss;
	    ss << gin.pressurescale.z_semi;
	    printIO("PRESSURESCALE", "Z_SEMI", ss.str(), "0..2");
	  }
          if (gin.boundcond.found && gin.boundcond.ntb == 0 &&
                  gin.pressurescale.couple == 2)
            printError("You cannot do pressure scaling if NTB=0 (vacuum) in BOUNDCOND block");
          if (gin.readtraj.found && gin.readtraj.ntrd != 0 &&
                  gin.pressurescale.couple == 2)
            printError("You cannot do pressure scaling when reading in a trajectory (NTRD!=0 in READTRAJ block)");
        }
        //
        // PRINTOUT block
        if (gin.printout.found) {
	  if (gin.printout.ntpr < 0) {
	    std::stringstream ss;
	    ss << gin.printout.ntpr;
	    printIO("PRINTOUT", "NTPR", ss.str(), ">= 0");
	  }
	  if (gin.printout.ntpp < 0 || gin.printout.ntpp > 1) {
	    std::stringstream ss;
	    ss << gin.printout.ntpp;
	    printIO("PRINTOUT", "NTPP", ss.str(), "0,1");
	  }
          if (gin.step.found && gin.step.nstlim < gin.printout.ntpr &&
                  gin.printout.ntpr != 0)
            printWarning("Printing of output less often than the number of steps");
	}
	//
	// RANDOMNUMBERS block
	if (gin.randomnumbers.found) {
	  if (gin.randomnumbers.ntrng < 0 || gin.randomnumbers.ntrng > 1) {
	    std::stringstream ss;
	    ss << gin.randomnumbers.ntrng;
	    printIO("RANDOMNUMBERS", "NTRNG", ss.str(), "0,1");
	  }
	  if (gin.randomnumbers.ntgsl < -1) {
	    std::stringstream ss;
	    ss << gin.randomnumbers.ntgsl;
	    printIO("RANDOMNUMBERS", "NTGSL", ss.str(), "> -1");
	  }
	}
        //
        // READTRAJ block
        if (gin.readtraj.found) {
	  if (gin.readtraj.ntrd < 0 || gin.readtraj.ntrd > 1) {
	    std::stringstream ss;
	    ss << gin.readtraj.ntrd;
	    printIO("READTRAJ", "NTRD", ss.str(), "0,1");
	  }
	  if (gin.readtraj.ntrn < 1 || gin.readtraj.ntrn > 18) {
	    std::stringstream ss;
	    ss << gin.readtraj.ntrn;
	    printIO("READTRAJ", "NTRN", ss.str(), "1..18");
	  }
	  if (gin.readtraj.ntrb < 0 || gin.readtraj.ntrb > 1) {
	    std::stringstream ss;
	    ss << gin.readtraj.ntrb;
	    printIO("READTRAJ", "NTRB", ss.str(), "0,1");
	  }
	  if (gin.readtraj.ntshk < 0 || gin.readtraj.ntshk > 1) {
	    std::stringstream ss;
	    ss << gin.readtraj.ntshk;
	    printIO("READTRAJ", "NTSHK", ss.str(), "0,1");
	  }
          if (gin.boundcond.found && gin.boundcond.ntb == 0 &&
                  gin.readtraj.ntrd != 0 && gin.readtraj.ntrb != 0)
            printError("Reading of boxdimensions from trajectory (NTRB in READTRAJ block) not possible with vacuum boundary conditions (NTB=0 in BOUNDCOND block)");
          if (gin.initialise.found && gin.initialise.ntivel != 0 &&
                  gin.readtraj.ntrd != 0)
            printError("If NTRD!=0 in READTRAJ block, you should not generate new velocities (NTIVEL in INITIALISE block)");
          if (gin.initialise.found && gin.initialise.ntishk > 1 &&
                  gin.readtraj.ntrd != 0)
            printError("If NTRD!=0 in READTRAJ block, you should not shake velocities (NTISHK in INITIALISE block)");
          if (gin.overalltransrot.found && gin.overalltransrot.ncmtr != 0)
            printError("If NTRD!=0 in READTRAJ block, you cannot use NCMTR!=0 in OVERALLTRANSROT block");
          if (gin.overalltransrot.found && gin.overalltransrot.ncmro != 0)
            printError("If NTRD!=0 in READTRAJ block, you cannot use NCMRO!=0 in OVERALLTRANSROT block");
          if (gin.virial.found && gin.virial.ntv != 0)
            printError("Reading of trajectory (NTRD!=0 in READTRAJ block), does not allow for a calculation of the virial (NTV!=0 in VIRIAL block)");
          if (gin.thermostat.found && gin.thermostat.ntt != 0)
            printError("Reading of trajectory (NTRD!=0 in READTRAJ block), does not allow for temperature scaling (NTT!=0 in THERMOSTAT block)");
          if (gin.writetraj.found && gin.writetraj.ntwv != 0)
            printError("When reading a trajectory file (NTRD!=0 in READTRAJ block), you cannot write a velocity trajectory (NTWV!=0 in WRITETRAJ block)");
        }
        //
        // REPLICA block
        if (gin.replica.found) {
	  for (unsigned int i = 0; i < gin.replica.ret.size(); ++i) {
	    if (gin.replica.ret[i] < 0) {
	      std::stringstream ss, si;
	      ss << gin.replica.ret[i];
	      si << "RET[" << i + 1 << "]";
	      printIO("REPLICA", si.str(), ss.str(), ">= 0.0");
	    }
	  }
	  if (gin.replica.lrescale < 0 || gin.replica.lrescale > 1) {
	    std::stringstream ss;
	    ss << gin.replica.lrescale;
	    printIO("REPLICA", "LRESCALE", ss.str(), "0,1");
	  }
	  for (unsigned int i = 0; i < gin.replica.relam.size(); ++i) {
	    if (gin.replica.relam[i] < 0) {
	      std::stringstream ss, si;
	      ss << gin.replica.relam[i];
	      si << "RELAM[" << i + 1 << "]";
	      printIO("REPLICA", si.str(), ss.str(), ">= 0.0");
	    }
	    if (gin.replica.rets[i] < 0) {
	      std::stringstream ss, si;
	      ss << gin.replica.rets[i];
	      si << "RETS[" << i + 1 << "]";
	      printIO("REPLICA", si.str(), ss.str(), ">= 0.0");
	    }
	  }
	  if (gin.replica.nretrial < 0) {
	    std::stringstream ss;
	    ss << gin.replica.nretrial;
	    printIO("REPLICA", "NRETRIAL", ss.str(), ">= 0");
	  }
	  if (gin.replica.nrequil < 0) {
	    std::stringstream ss;
	    ss << gin.replica.nrequil;
	    printIO("REPLICA", "NREQUIL", ss.str(), ">= 0");
	  }
	  if (gin.replica.nrejob <= 0) {
	    std::stringstream ss;
	    ss << gin.replica.nrejob;
	    printIO("REPLICA", "NREJOB", ss.str(), "> 0");
	  }
	  if (gin.replica.nrewrt < 0) {
	    std::stringstream ss;
	    ss << gin.replica.nrewrt;
	    printIO("REPLICA", "NREWRT", ss.str(), ">= 0");
	  }
          if (!gin.perturbation.found ||
                  ((gin.perturbation.found && gin.perturbation.ntg == 0) &&
                  gin.replica.relam.size()))
            printError("You cannot do a Hamiltonian REMD simulation without specifying a perturbation");
	}
	//
	// ROTTRANS block
	if (gin.rottrans.found) {
	  if (gin.rottrans.rtc < 0 || gin.rottrans.rtc > 1) {
	    std::stringstream ss;
	    ss << gin.rottrans.rtc;
	    printIO("ROTTRANS", "RTC", ss.str(), "0,1");
	  }
	  if (gin.rottrans.rtclast <= 0) {
	    std::stringstream ss;
	    ss << gin.rottrans.rtclast;
	    printIO("ROTTRANS", "RTCLAST", ss.str(), "> 0");
	  }
	}
	//
	// STEP block
	if (gin.step.found) {
	  if (gin.step.nstlim < 0) {
	    std::stringstream ss;
	    ss << gin.step.nstlim;
	    printIO("STEP", "NSTLIM", ss.str(), ">= 0");
	  }
	  if (gin.step.t < 0.0) {
	    std::stringstream ss;
	    ss << gin.step.t;
	    printIO("STEP", "T", ss.str(), ">= 0.0");
	  }
	  if (gin.step.dt <= 0) {
	    std::stringstream ss;
	    ss << gin.step.dt;
	    printIO("STEP", "DT", ss.str(), "> 0.0");
	  }
	}
	//
        // STOCHDYN block
	if (gin.stochdyn.found) {
	  if (gin.stochdyn.ntsd < 0 || gin.stochdyn.ntsd > 1) {
	    std::stringstream ss;
	    ss << gin.stochdyn.ntsd;
	    printIO("STOCHDYN", "NTSD", ss.str(), "0,1");
	  }
	  if (gin.stochdyn.ntfr < 0 || gin.stochdyn.ntfr > 3) {
	    std::stringstream ss;
	    ss << gin.stochdyn.ntfr;
	    printIO("STOCHDYN", "NTFR", ss.str(), "0,1,2,3");
	  }
	  if (gin.stochdyn.nsfr <= 0) {
	    std::stringstream ss;
	    ss << gin.stochdyn.nsfr;
	    printIO("STOCHDYN", "NSFR", ss.str(), "> 0");
	  }
	  if (gin.stochdyn.nbref <= 0) {
	    std::stringstream ss;
	    ss << gin.stochdyn.nbref;
	    printIO("STOCHDYN", "NBREF", ss.str(), "> 0");
	  }
	  if (gin.stochdyn.rcutf < 0) {
	    std::stringstream ss;
	    ss << gin.stochdyn.rcutf;
	    printIO("STOCHDYN", "RCUTF", ss.str(), ">= 0.0");
	  }
	  if (gin.stochdyn.cfric < 0) {
	    std::stringstream ss;
	    ss << gin.stochdyn.cfric;
	    printIO("STOCHDYN", "CFRIC", ss.str(), ">= 0.0");
	  }
	  if (gin.stochdyn.tempsd < 0) {
	    std::stringstream ss;
	    ss << gin.stochdyn.tempsd;
	    printIO("STOCHDYN", "TEMPSD", ss.str(), ">= 0.0");
	  }
	}
        if (gin.stochdyn.found && gin.stochdyn.ntsd != 0) {
          if (gin.overalltransrot.found && gin.overalltransrot.ncmtr != 0)
            printWarning("NCMTR!=0 in OVERALLTRANSROT block with SD is not advisable");
          if (gin.overalltransrot.found && gin.overalltransrot.ncmro != 0)
            printWarning("NCMRO!=0 in OVERALLTRANSROT block with SD is not advisable");
          if (gin.readtraj.found && gin.readtraj.ntrd != 0)
            printError("Reading of a trajectory file (NTRD!=0 in READTRAJ) block is not possible with SD (NTSD!=0 in STOCHDYN block)");
        }
        //
        //SYSTEM block
        if (gin.system.found) {
	  if (gin.system.npm < 0) {
	    std::stringstream ss;
	    ss << gin.system.npm;
	    printIO("SYSTEM", "NPM", ss.str(), ">= 0");
	  }
	  if (gin.system.nsm < 0) {
	    std::stringstream ss;
	    ss << gin.system.nsm;
	    printIO("SYSTEM", "NSM", ss.str(), ">= 0");
	  }
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
                  printError(os.str());
                }
              }
            }
          }
          if (gin.system.npm == 0 && gin.system.nsm == 0)
            printError("NPM and NSM in SYSTEM block cannot both be zero");
          if (gromosXX && gin.system.npm > 1)
            printWarning("for md++, NPM in SYSTEM block cannot be > 1");

          if (gin.multicell.found && gin.multicell.ntm != 0) {
            int num = gin.multicell.ncella *
                    gin.multicell.ncellb *
                    gin.multicell.ncellc;
            if (!gromosXX && num % gin.system.npm != 0)
              printError("if NTM != 0 in MULTICELL block, then NCELLX * NCELLY * NCELLZ has to be a multiple of NPM");
            if (!gromosXX && num % gin.system.nsm != 0)
              printError("if NTM != 0 in MULTICELL block, then NCELLX * NCELLY * NCELLZ has to be a multiple of NSM");
          }
          if (gin.localelev.found && gin.localelev.ntles != 0 && gin.system.npm != 1)
            printError("if NTLES != 0 in LOCALELEV block, NPM in SYSTEM block has to be 1");
          if (gin.positionres.found && gin.positionres.ntpor == 3 &&
                  gin.system.npm != 1)
            printError("if NTPOR = 3 in POSITIONRES block, NPM in SYSTEM block has to be 1");
          if (gin.distanceres.found && gin.distanceres.ntdir != 0 &&
                  gin.system.npm == 0)
            printError("if NTDIR != 0 in DISTANCERES block, NPM in SYSTEM block cannot be 0");
          if (gin.dihedralres.found && gin.dihedralres.ntdlr != 0 &&
                  gin.system.npm == 0)
            printError("if NTDLR != 0 in DIHEDRALRES block, NPM in SYSTEM block cannot be 0");
        }
        //
        // THERMOSTAT block
	if (gin.thermostat.found) {
	  if (gin.thermostat.ntt < 0 || gin.thermostat.ntt > 1) {
	    std::stringstream ss;
	    ss << gin.thermostat.ntt;
	    printIO("THERMOSTAT", "NTT", ss.str(), "0,1");
	  }
	  for (unsigned int i = 0; i < gin.thermostat.baths.size(); ++i) {
	    if (gin.thermostat.baths[i].ntbtyp < 0 || gin.thermostat.baths[i].ntbtyp > 3) {
	      std::stringstream ss, si;
	      ss << gin.thermostat.baths[i].ntbtyp;
	      si << "NTBTYP[" << i + 1 << "]";
	      printIO("THERMOSTAT", si.str(), ss.str(), "0,1,2,3");
	    }
	    if (gin.thermostat.baths[i].tembth < 0.0) {
	      std::stringstream ss, si;
	      ss << gin.thermostat.baths[i].tembth;
	      si << "TEMBTH[" << i + 1 << "]";
	      printIO("THERMOSTAT", si.str(), ss.str(), ">= 0.0");
	    }
	    for (unsigned int j = 0; j < gin.thermostat.baths[i].taubth.size(); ++j) {
	      if (gin.thermostat.baths[i].taubth[j] < 0.0) {
		std::stringstream ss, si;
		ss << gin.thermostat.baths[i].taubth[j];
		si << "TAUBTH[" << i + 1 << "," << j + 1 << "]";
		printIO("THERMOSTAT", si.str(), ss.str(), ">= 0.0");
	      }
	    }
	  }
	  for (unsigned int i = 0; i < gin.thermostat.dofgroups.size(); ++i) {
	    // class ithermostat::dofgroup d;
	    if (gin.thermostat.dofgroups[i].ntscpl < 0 || gin.thermostat.dofgroups[i].ntscpl > (int) gin.thermostat.baths.size()) {
	      std::stringstream ss, si, sj;
	      ss << gin.thermostat.dofgroups[i].ntscpl;
	      si << "NTSCPL[" << i + 1 << "]";
	      sj << "0..NTBTH (" << gin.thermostat.baths.size() << ")";
	      printIO("THERMOSTAT", si.str(), ss.str(), sj.str());
	    }
	    if (gin.thermostat.dofgroups[i].ntstyp < 0 || gin.thermostat.dofgroups[i].ntstyp > 2) {
	      std::stringstream ss, si;
	      ss << gin.thermostat.dofgroups[i].ntstyp;
	      si << "NTSTYP[" << i + 1 << "]";
	      printIO("THERMOSTAT", si.str(), ss.str(), "0,1,2");
	    }
	    if (gin.thermostat.dofgroups[i].ntscns < 0 || gin.thermostat.dofgroups[i].ntscns > 1) {
	      std::stringstream ss, si;
	      ss << gin.thermostat.dofgroups[i].ntscns;
	      si << "NTSCNS[" << i + 1 << "]";
	      printIO("THERMOSTAT", si.str(), ss.str(), "0,1");
	    }
	    if (gin.thermostat.dofgroups[i].ntsgt < -2) {
	      std::stringstream ss, si;
	      ss << gin.thermostat.dofgroups[i].ntsgt;
	      si << "NTSGT[" << i + 1 << "]";
	      printIO("THERMOSTAT", si.str(), ss.str(), ">= -2");
	    }
	    int gt = 0;
	    if (gin.thermostat.dofgroups[i].ntsgt <= 0) gt = 0;
	    else gt = gin.thermostat.dofgroups[i].ntsgt;
	    for (int k = 0; k < gt; k++) {
	      if (gin.thermostat.dofgroups[i].ntsgtg[k] <= 0) {
		std::stringstream ss, si;
		ss << gin.thermostat.dofgroups[i].ntsgtg[k];
		si << "NTSGTG[" << i + 1 << "," << k + 1 << "]";
		printIO("THERMOSTAT", si.str(), ss.str(), ">= 1");
	      }
	    }
	  }
	}
        if (gin.thermostat.found && gin.thermostat.ntt != 0) {
          if (gin.thermostat.dofgroups.size() == 0)
            printError("You cannot to temperature coupling with NTSET=0 in THERMOSTAT block");
          // Someone who actually understands this block may want to add
          // more tests.
        }
        //
	// UMBRELLA block
	if (gin.umbrella.found) {
	  if (gin.umbrella.ntus < 0 || gin.umbrella.ntus > 1) {
	    std::stringstream ss;
	    ss << gin.umbrella.ntus;
	    printIO("UMBRELLA", "NTUS", ss.str(), "0,1");
	  }
	  if (gin.umbrella.uscst1 < 0.0) {
	    std::stringstream ss;
	    ss << gin.umbrella.uscst1;
	    printIO("UMBRELLA", "USCST1", ss.str(), ">= 0.0");
	  }
	  if (gin.umbrella.uscst2 < 0.0) {
	    std::stringstream ss;
	    ss << gin.umbrella.uscst2;
	    printIO("UMBRELLA", "USCST2", ss.str(), ">= 0.0");
	  }
	  if (gin.umbrella.usref1 < 0.0) {
	    std::stringstream ss;
	    ss << gin.umbrella.usref1;
	    printIO("UMBRELLA", "USREF1", ss.str(), ">= 0.0");
	  }
	  if (gin.umbrella.usref2 < 0.0) {
	    std::stringstream ss;
	    ss << gin.umbrella.usref2;
	    printIO("UMBRELLA", "USREF2", ss.str(), ">= 0.0");
	  }
	}
	//
        // VIRIAL block
	if (gin.virial.found) {
	  if (gin.virial.ntv < 0 || gin.virial.ntv > 1) {
	    std::stringstream ss;
	    ss << gin.virial.ntv;
	    printIO("VIRIAL", "NTV", ss.str(), "0,1");
	  }
	  if (gin.virial.ntvg < 0 || gin.virial.ntvg > 1) {
	    std::stringstream ss;
	    ss << gin.virial.ntvg;
	    printIO("VIRIAL", "NTVG", ss.str(), "0,1,2,3");
	  }
	}
        if (gin.virial.found && gin.virial.ntv != 0) {
          if (gin.readtraj.found && gin.readtraj.ntrd != 0)
            printError("Reading of trajectory (NTRD!=0 in READTRAJ block) does not support recalculation of the virial (NTV in VIRIAL block");
        }
        //
        // WRITETRAJ block
        if (gin.writetraj.found) {
	  if (gin.writetraj.ntwse < 0.0) {
	    std::stringstream ss;
	    ss << gin.writetraj.ntwse;
	    printIO("WRITETRAJ", "NTWSE", ss.str(), ">= 0.0");
	  }
	  if (gin.writetraj.ntwf < 0.0) {
	    std::stringstream ss;
	    ss << gin.writetraj.ntwf;
	    printIO("WRITETRAJ", "NTWF", ss.str(), ">= 0.0");
	  }
	  if (gin.writetraj.ntwe < 0.0) {
	    std::stringstream ss;
	    ss << gin.writetraj.ntwe;
	    printIO("WRITETRAJ", "NTWE", ss.str(), ">= 0.0");
	  }
	  if (gin.writetraj.ntwg < 0.0) {
	    std::stringstream ss;
	    ss << gin.writetraj.ntwg;
	    printIO("WRITETRAJ", "NTWG", ss.str(), ">= 0.0");
	  }
	  if (gin.writetraj.ntwb < 0.0) {
	    std::stringstream ss;
	    ss << gin.writetraj.ntwb;
	    printIO("WRITETRAJ", "NTWB", ss.str(), ">= 0.0");
	  }
          if ((!gin.perturbation.found ||
                  (gin.perturbation.found && gin.perturbation.ntg == 0)) &&
                  gin.writetraj.ntwg != 0)
            printError("You cannot write a free energy trajectorey (NTWG!=0 in WRITETRAJ block) without doing a perturbation");
          if (gin.writetraj.ntwse != 0 && gin.writetraj.ntwx == 0)
            printError("NTWSE!=0 in WRITETRAJ block requires NTWX!=0");
          if (gin.writetraj.ntwse != 0 && gin.writetraj.ntwv != 0)
            printError("NTWSE!=0 in WRITETRAJ block requires NTWV=0");
          if (gin.writetraj.ntwse != 0 && gin.writetraj.ntwf != 0)
            printError("NTWSE!=0 in WRITETRAJ block requires NTWF=0");
          if (gin.writetraj.ntwse != 0 && gin.writetraj.ntwe != 0)
            printError("NTWSE!=0 in WRITETRAJ block requires NTWE=0");
          if (gin.writetraj.ntwse != 0 && gin.writetraj.ntwg != 0)
            printError("NTWSE!=0 in WRITETRAJ block requires NTWG=0");
          if (gin.writetraj.ntwse != 0 && gin.writetraj.ntwb != 0)
            printError("NTWSE!=0 in WRITETRAJ block requires NTWB=0");
          if (gin.energymin.found && gin.energymin.ntem != 0 &&
                  gin.writetraj.ntwv != 0)
            printWarning("NTWV!=0 in WRITETRAJ block for an energy minimisation does not make sense");
	}
	//
	// XRAYRES block
	if (gin.xrayres.found) {
	  if (gin.xrayres.ntxr < -2 || gin.xrayres.ntxr > 3) {
	    std::stringstream ss;
	    ss << gin.xrayres.ntxr;
	    printIO("XRAYRES", "NTXR", ss.str(), "-2,-1,0,1,2,3");
	  }
	  if (gin.xrayres.ntxle < 0 || gin.xrayres.ntxle > 1) {
	    std::stringstream ss;
	    ss << gin.xrayres.ntxr;
	    printIO("XRAYRES", "NTXLE", ss.str(), "0,1");
          }
	  if (gin.xrayres.cxr < 0) {
	    std::stringstream ss;
	    ss << gin.xrayres.cxr;
	    printIO("XRAYRES", "CXR", ss.str(), ">= 0.0");
	  }
	  if (gin.xrayres.ntwxr < 0) {
	    std::stringstream ss;
	    ss << gin.xrayres.ntwxr;
	    printIO("XRAYRES", "NTWXR", ss.str(), ">= 0");
	  }
	  if (gin.xrayres.ntwde < 0 || gin.xrayres.ntwde > 3) {
	    std::stringstream ss;
	    ss << gin.xrayres.ntwde;
	    printIO("XRAYRES", "NTWDE", ss.str(), "0,1,2,3");
	  }
	  if (gin.xrayres.ntwxm < 0) {
	    std::stringstream ss;
	    ss << gin.xrayres.ntwxm;
	    printIO("XRAYRES", "NTWXM", ss.str(), ">= 0");
	  }
	  if (gin.xrayres.cxtau < 0) {
	    std::stringstream ss;
	    ss << gin.xrayres.cxtau;
	    printIO("XRAYRES", "CXTAU", ss.str(), ">= 0.0");
	  }
	  if (gin.xrayres.rdavg < 0 || gin.xrayres.rdavg > 1) {
	    std::stringstream ss;
	    ss << gin.xrayres.rdavg;
	    printIO("XRAYRES", "RDAVG", ss.str(), "0,1");
	  }
	}
      }
      // *********************************************
      // FINISHED THE CROSS CHECKS AND INPUT CHECKS
      // *********************************************

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
                  queue);
          fout << filenames[FILETYPE["coord"]].name(0) << endl;
          filenames[FILETYPE["coord"]].setInfo(systemname,
                  atof(iter->second.param["T"].c_str()),
                  atof(iter->second.param["DELTAT"].c_str()),
                  iter->first,
                  queue);
        }
      }
      // re-analyzing?
      if (gin.readtraj.found && gin.readtraj.ntrd == 1) {
        fout << "ANATRX=${SIMULDIR}/" << filenames[FILETYPE["anatrj"]].name(0) << endl;
      }
      // EVTRL
      if (l_refpos) fout << "REFPOS=${SIMULDIR}/" << s_refpos << endl;
      if (l_posresspec) fout << "POSRESSPEC=${SIMULDIR}/"
              << s_posresspec << endl;
      if (l_xray) fout << "XRAY=${SIMULDIR}/" << s_xray << endl;
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
      if (gin.writetraj.ntwx)
        fout << "OUTPUTTRX="
              << filenames[FILETYPE["outtrx"]].name(0)
        << endl;
      if (gin.writetraj.ntwv)
        fout << "OUTPUTTRV="
              << filenames[FILETYPE["outtrv"]].name(0)
        << endl;
      if (gin.writetraj.ntwf)
        fout << "OUTPUTTRF="
              << filenames[FILETYPE["outtrf"]].name(0)
        << endl;
      if (gin.writetraj.ntwe)
        fout << "OUTPUTTRE="
              << filenames[FILETYPE["outtre"]].name(0)
        << endl;
      if (gin.writetraj.ntwg)
        fout << "OUTPUTTRG="
              << filenames[FILETYPE["outtrg"]].name(0)
        << endl;
      if (gin.writetraj.ntwb)
        fout << "OUTPUTBAE="
              << filenames[FILETYPE["outbae"]].name(0)
        << endl;
      if (gin.polarise.write || gin.jvalueres.write || gin.xrayres.ntwxr)
        fout << "OUTPUTTRS="
              << filenames[FILETYPE["outtrs"]].name(0)
        << endl;

      if (gin.writetraj.ntwb &&
              (gin.perturbation.found && gin.perturbation.ntg > 0))
        fout << "OUTPUTBAG="
              << filenames[FILETYPE["outbag"]].name(0)
        << endl;

      // any additional links?
      for (unsigned int k = 0; k < linkadditions.size(); k++)
        if (linkadditions[k] > 0)
          fout << linknames[k] << "="
                << filenames[numFiletypes + k].name(0) << endl;

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
        if (gin.writetraj.ntwx) fout << setw(25)
          << "ln -s ${OUTPUTTRX}" << " fort.12\n";
        if (gin.writetraj.ntwv) fout << setw(25)
          << "ln -s ${OUTPUTTRV}" << " fort.13\n";
        if (gin.writetraj.ntwf) fout << setw(25)
          << "ln -s ${OUTPUTTRF}" << " fort.13\n";
        if (gin.writetraj.ntwe) fout << setw(25)
          << "ln -s ${OUTPUTTRE}" << " fort.15\n";
        if (gin.writetraj.ntwg) fout << setw(25)
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
        if (gin.readtraj.found && gin.readtraj.ntrd == 1) fout << " \\\n\t"
                << setw(12) << "@anatrj" << " ${ANATRX}";
        if (l_pttopo) fout << " \\\n\t"
                << setw(12) << "@pttopo" << " ${PTTOPO}";
        if (l_posresspec) fout << " \\\n\t"
                << setw(12) << "@posresspec" << " ${POSRESSPEC}";
        if (l_refpos) fout << " \\\n\t"
                << setw(12) << "@refpos" << " ${REFPOS}";
        if (l_xray) fout << " \\\n\t"
                << setw(12) << "@xray" << " ${XRAY}";
        if (l_disres) fout << " \\\n\t"
                << setw(12) << "@distrest" << " ${DISRES}";
        if (l_dihres) fout << " \\\n\t"
                << setw(12) << "@dihrest" << " ${DIHRES}";
        if (l_jvalue) fout << " \\\n\t"
                << setw(12) << "@jval" << " ${JVALUE}";

        fout << " \\\n\t" << setw(12) << "@fin" << " ${OUTPUTCRD}";
        if (gin.writetraj.ntwx) fout << " \\\n\t" << setw(12) << "@trj"
          << " ${OUTPUTTRX}";
        if (gin.writetraj.ntwv) fout << " \\\n\t" << setw(12) << "@trv"
          << " ${OUTPUTTRV}";
        if (gin.writetraj.ntwf) fout << " \\\n\t" << setw(12) << "@trf"
          << " ${OUTPUTTRF}";
        if (gin.writetraj.ntwe) fout << " \\\n\t" << setw(12) << "@tre"
          << " ${OUTPUTTRE}";
        if (gin.writetraj.ntwg) fout << " \\\n\t" << setw(12) << "@trg"
          << " ${OUTPUTTRG}";
        if (gin.writetraj.ntwb) fout << " \\\n\t" << setw(12) << "@bae"
          << " ${OUTPUTBAE}";
        if (gin.polarise.write || gin.jvalueres.write || gin.xrayres.ntwxr)
          fout << " \\\n\t" << setw(12) << "@trs ${OUTPUTTRS}";

        if (gin.writetraj.ntwb > 0 &&
                gin.perturbation.found && gin.perturbation.ntg > 0)
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

        if (gromosXX) {
          fout << "\\\n\t" << setw(12) << ">" << " ${OUNIT}\n";
          fout << "grep \"finished successfully\" ${OUNIT} > /dev/null || MDOK=0";
        } else {
          fout << "\\\n\t" << setw(12) << ">" << " ${OUNIT}     || MDOK=0";
        }
        fout << "\n\n";
      }

      fout << "uname -a >> ${OUNIT}\n";

      if (gin.writetraj.ntwx || gin.writetraj.ntwv || gin.writetraj.ntwf ||
              gin.writetraj.ntwe || gin.writetraj.ntwg)
        fout << "\n# compress some files\n";
      if (gin.writetraj.ntwx) fout << "gzip ${OUTPUTTRX}\n";
      if (gin.writetraj.ntwv) fout << "gzip ${OUTPUTTRV}\n";
      if (gin.writetraj.ntwf) fout << "gzip ${OUTPUTTRF}\n";
      if (gin.writetraj.ntwe) fout << "gzip ${OUTPUTTRE}\n";
      if (gin.writetraj.ntwg) fout << "gzip ${OUTPUTTRG}\n";
      if (gin.writetraj.ntwb) fout << "gzip ${OUTPUTBAE}\n";
      if (gin.writetraj.ntwb &&
              gin.perturbation.found && gin.perturbation.ntg > 0)
        fout << "gzip ${OUTPUTBAG}\n";
      if (gin.polarise.write || gin.jvalueres.write || gin.xrayres.ntwxr)
        fout << "gzip ${OUTPUTTRS}\n";

      fout << "\n# copy the files back\n";
      fout << "OK=1\n";
      fout << setw(25) << "cp ${OUNIT}" << " ${SIMULDIR}";
      if (iter->second.dir != ".") fout << "/" << iter->second.dir;
      fout << " || OK=0\n";
      fout << setw(25) << "cp ${OUTPUTCRD}" << " ${SIMULDIR}";
      if (iter->second.dir != ".") fout << "/" << iter->second.dir;
      fout << " || OK=0\n";
      if (gin.writetraj.ntwx) {
        fout << setw(25) << "cp ${OUTPUTTRX}.gz" << " ${SIMULDIR}";
        if (iter->second.dir != ".") fout << "/" << iter->second.dir;
        fout << " || OK=0\n";
      }
      if (gin.writetraj.ntwv) {
        fout << setw(25) << "cp ${OUTPUTTRV}.gz" << " ${SIMULDIR}";
        if (iter->second.dir != ".") fout << "/" << iter->second.dir;
        fout << " || OK=0\n";
      }
      if (gin.writetraj.ntwf) {
        fout << setw(25) << "cp ${OUTPUTTRF}.gz" << " ${SIMULDIR}";
        if (iter->second.dir != ".") fout << "/" << iter->second.dir;
        fout << " || OK=0\n";
      }
      if (gin.writetraj.ntwe) {
        fout << setw(25) << "cp ${OUTPUTTRE}.gz" << " ${SIMULDIR}";
        if (iter->second.dir != ".") fout << "/" << iter->second.dir;
        fout << " || OK=0\n";
      }
      if (gin.writetraj.ntwg) {
        fout << setw(25) << "cp ${OUTPUTTRG}.gz" << " ${SIMULDIR}";
        if (iter->second.dir != ".") fout << "/" << iter->second.dir;
        fout << " || OK=0\n";
      }
      if (gin.writetraj.ntwb) {
        fout << setw(25) << "cp ${OUTPUTBAE}.gz" << " ${SIMULDIR}";
        if (iter->second.dir != ".") fout << "/" << iter->second.dir;
        fout << " || OK=0\n";

        if (gin.perturbation.found && gin.perturbation.ntg > 0) {

          fout << setw(25) << "cp ${OUTPUTBAG}.gz" << " ${SIMULDIR}";
          if (iter->second.dir != ".") fout << "/" << iter->second.dir;
          fout << " || OK=0\n";
        }
      }
      if (gin.polarise.write || gin.jvalueres.write || gin.xrayres.ntwxr) {
        fout << setw(25) << "cp ${OUTPUTTRS}.gz" << " ${SIMULDIR}";
        if (iter->second.dir != ".") fout << "/" << iter->second.dir;
        fout << " || OK=0\n";
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
      if (misc[4].name(0) != "") {
        fout << "  " << misc[4].name(0) << "\n";
      }
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
                  gin.step.dt * gin.step.nstlim, it->first, queue);
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

void printIO(string block, string variable, string value, string allowed) {
  numErrors++;
  cout << "INPUT ERROR\n";
  cout << numWarnings + numErrors << ". ERROR [IO check] (" << numErrors
          << ")\n";
  cout << "Error in block " << block << "\n";
  cout << "Read " << value << " for " << variable << "\n";
  cout << "Accepted values are " << allowed << endl;
}

void printWarning(string s) {
  numWarnings++;
  cout << numWarnings + numErrors << ". WARNING (" << numWarnings << ")\n";
  cout << s;
  cout << endl;
}

void printError(string s) {
  numErrors++;
  cout << numWarnings + numErrors << ". ERROR (" << numErrors << ")\n";
  cout << s;
  cout << endl;
}

void printInput(string ofile, input gin) {
  ofstream fout(ofile.c_str());
  const time_t t = time(0);
  fout << "TITLE\n";
  fout << "\tAutomatically generated input file\n\t";
  fout << getenv("USER") << " " << ctime(&t);
  fout << "END\n";
  fout << gin;
}

void readJobinfo(string file, map<int, jobinfo> &ji) {
  Ginstream gin(file);
  vector<string> buffer;
  gin.getblock(buffer);
  if (buffer[0] != "JOBSCRIPTS")
    throw gromos::Exception("mk_script", "Reading of jobscript file failed. "
          "No JOBSCRIPTS block");
  istringstream iss(buffer[1]);
  vector<string> head;
  string b;
  while ((iss >> b) != 0) head.push_back(b);
  if (head[0] != "job_id" || head.back() != "run_after"
          || head[head.size() - 2] != "subdir")
    throw gromos::Exception("mk_script", "Reading of jobscript file failed.\n"
          "First line syntax:\n"
          "job_id PARAM PARAM ... subdir run_after");
  if (buffer.back().find("END") != 0)
    throw gromos::Exception("mk_script", "Jobscript file "
          + gin.name() +
          " is corrupted. No END in JOBSCRIPTS"
          " block. Got\n"
          + buffer.back());
  int id = 0;
  for (unsigned int i = 2; i < buffer.size() - 1; i++) {
    vector<string> tmp(head.size());
    iss.clear();
    iss.str(buffer[i]);

    for (unsigned int j = 0; j < head.size(); j++) iss >> tmp[j];
    jobinfo job;
    id = atoi(tmp[0].c_str());
    for (unsigned int j = 1; j < head.size() - 2; j++)
      job.param[head[j]] = tmp[j];
    job.dir = tmp[head.size() - 2];
    job.prev_id = atoi(tmp.back().c_str());
    ji[id] = job;
  }
}

void readLibrary(string file, vector<filename> &names,
        vector<filename> &misc,
        vector<string> &linknames, vector<int> &linkadditions,
        string system, string queue, double t,
        double dt, int ns) {
  // Open the file
  Ginstream templates(file);
  int found_filenames = 0;
  string sdum, temp, first;
  templates.getline(first);

  while (!templates.stream().eof()) {
    vector<string> buffer;
    templates.getblock(buffer);

    if (buffer.size() && first == "FILENAMES") {
      if (buffer[buffer.size() - 1].find("END") != 0)
        throw gromos::Exception("mk_script", "Template file "
              + templates.name() +
              " is corrupted. No END in " + first +
              " block. Got\n"
              + buffer[buffer.size() - 1]);
      for (unsigned int j = 0; j < buffer.size() - 1; j++) {
        found_filenames = 1;
        istringstream iss(buffer[j]);
        iss >> sdum >> temp;
        switch (FILETYPE[sdum]) {
          case inputfile: names[inputfile].setTemplate(temp);
            break;
          case topofile: names[topofile].setTemplate(temp);
            break;
          case coordfile: names[coordfile].setTemplate(temp);
            break;
          case refposfile: names[refposfile].setTemplate(temp);
            break;
          case anatrxfile: names[anatrxfile].setTemplate(temp);
            break;
          case posresspecfile: names[posresspecfile].setTemplate(temp);
            break;
          case xrayfile: names[xrayfile].setTemplate(temp);
            break;
          case disresfile: names[disresfile].setTemplate(temp);
            break;
          case pttopofile: names[pttopofile].setTemplate(temp);
            break;
          case dihresfile: names[dihresfile].setTemplate(temp);
            break;
          case jvaluefile: names[jvaluefile].setTemplate(temp);
            break;
          case ledihfile: names[ledihfile].setTemplate(temp);
            break;
          case outputfile: names[outputfile].setTemplate(temp);
            break;
          case outtrxfile: names[outtrxfile].setTemplate(temp);
            break;
          case outtrvfile: names[outtrvfile].setTemplate(temp);
            break;
          case outtrffile: names[outtrffile].setTemplate(temp);
            break;
          case outtrefile: names[outtrefile].setTemplate(temp);
            break;
          case outtrgfile: names[outtrgfile].setTemplate(temp);
            break;
          case outbaefile: names[outbaefile].setTemplate(temp);
            break;
          case outbagfile: names[outbagfile].setTemplate(temp);
            break;
          case scriptfile: names[scriptfile].setTemplate(temp);
            break;
          case outtrsfile: names[outtrsfile].setTemplate(temp);
            break;
          case unknownfile:
            printWarning("Don't know how to handle template for " + sdum
                    + ". Ingoring");
        }
      }
    }
    if (buffer.size() && first == "MISCELLANEOUS") {
      if (buffer[buffer.size() - 1].find("END") != 0)
        throw gromos::Exception("mk_script", "Template file " +
              templates.name() +
              " is corrupted. No END in " + first +
              " block. Got\n"
              + buffer[buffer.size() - 1]);

      for (unsigned int j = 0; j < buffer.size() - 1; j++) {
        istringstream iss(buffer[j]);
        iss >> sdum;
        if (sdum == "workdir") {
          iss >> temp;
          misc[0].setTemplate(temp);
        }
        if (sdum == "lastcommand") {
          ostringstream os;
          while (!iss.eof()) {
            iss >> sdum;
            os << sdum << " ";
          }
          misc[1].setTemplate(os.str());
        }
        if (sdum == "firstcommand") {
          ostringstream os;
          while (!iss.eof()) {
            iss >> sdum;
            os << sdum << " ";
          }
          misc[2].setTemplate(os.str());
        }
        if (sdum == "mpicommand") {
          ostringstream os;
          while (!iss.eof()) {
            iss >> sdum;
            os << sdum << " ";
          }
          misc[3].setTemplate(os.str());
        }
        if (sdum == "stopcommand") {
          ostringstream os;
          while (!iss.eof()) {
            iss >> sdum;
            os << sdum << " ";
          }
          misc[4].setTemplate(os.str());
        }
      }
    }
    if (buffer.size() && first == "LINKADDITION") {
      if (buffer[buffer.size() - 1].find("END") != 0)
        throw gromos::Exception("mk_script", "Template file " +
              templates.name() +
              " is corrupted. No END in " + first +
              " block. Got\n"
              + buffer[buffer.size() - 1]);

      for (unsigned int j = 0; j < buffer.size() - 1; j++) {
        istringstream iss(buffer[j]);
        int k;
        string varname;

        iss >> sdum >> varname >> temp >> k;
        filename newlink(system, t, dt, ns, queue);
        newlink.setTemplate(temp);
        names.push_back(newlink);
        if (sdum == "input") k *= -1;
        linkadditions.push_back(k);
        linknames.push_back(varname);

      }
    }
    templates.getline(first);
  }
}

void setParam(input &gin, jobinfo const &job) {
  map<string, string>::const_iterator iter = job.param.begin(),
          to = job.param.end();
  for (; iter != to; ++iter) {
    if (iter->first == "ENDTIME")
      ; // do nothing but avoid warning
      // BAROSTAT
    else if (iter->first == "NTP")
      gin.barostat.ntp = atoi(iter->second.c_str());
    else if (gin.barostat.found && iter->first == "COMP")
      gin.barostat.comp = atof(iter->second.c_str());
    else if (iter->first.substr(0, 7) == "PRSBTH[") {
      unsigned int i = atoi(iter->first.substr(7, iter->first.find("]")).c_str());
      if (i <= gin.barostat.pbaths.size())
        gin.barostat.pbaths[i - 1].prsbth = atof(iter->second.c_str());
      else {
        std::stringstream ss;
        ss << iter->second.c_str() << " in joblistfile out of range";
        printError(ss.str());
      }
    }      // BOUNDCOND
    else if (iter->first == "NTB")
      gin.boundcond.ntb = atoi(iter->second.c_str());
    else if (iter->first == "NDFMIN")
      gin.boundcond.ndfmin = atoi(iter->second.c_str());

      // CGRAIN
    else if (iter->first == "NTCGRAN")
      gin.cgrain.ntcgran = atoi(iter->second.c_str());
    else if (iter->first == "EPS")
      gin.cgrain.eps = atof(iter->second.c_str());

      // COMTRANSROT
    else if (iter->first == "NSCM")
      gin.comtransrot.nscm = atoi(iter->second.c_str());

      // CONSISTENCYCHECK
    else if (iter->first == "NTCHK")
      gin.consistencycheck.ntchk = atoi(iter->second.c_str());
    else if (iter->first == "NTCKF")
      gin.consistencycheck.ntckf = atoi(iter->second.c_str());
    else if (iter->first == "FDCKF")
      gin.consistencycheck.fdckf = atof(iter->second.c_str());
    else if (iter->first == "NTCKV")
      gin.consistencycheck.ntckv = atoi(iter->second.c_str());
    else if (iter->first == "FDCKV")
      gin.consistencycheck.fdckv = atof(iter->second.c_str());
    else if (iter->first == "NTCKT")
      gin.consistencycheck.ntckt = atoi(iter->second.c_str());
    else if (iter->first == "NTCKE")
      gin.consistencycheck.ntcke = atoi(iter->second.c_str());
    else if (iter->first == "NTCKR")
      gin.consistencycheck.ntckr = atoi(iter->second.c_str());
    else if (iter->first == "NTCKL")
      gin.consistencycheck.ntckl = atoi(iter->second.c_str());
    else if (iter->first == "FDCKL")
      gin.consistencycheck.fdckl = atof(iter->second.c_str());

      // CONSTRAINT
    else if (iter->first == "NTC")
      gin.constraint.ntc = atoi(iter->second.c_str());
    else if (iter->first == "NTCP")
      gin.constraint.ntcp = atoi(iter->second.c_str());
    else if (iter->first == "NTCS")
      gin.constraint.ntcs = atoi(iter->second.c_str());

      // COVALENTFORM
    else if (iter->first == "NTBBH")
      gin.covalentform.ntbbh = atoi(iter->second.c_str());
    else if (iter->first == "NTBAH")
      gin.covalentform.ntbah = atoi(iter->second.c_str());
    else if (iter->first == "NTBDN")
      gin.covalentform.ntbdn = atoi(iter->second.c_str());

      // DEBUG
      // DIHEDRALRES
    else if (iter->first == "NTDLR")
      gin.dihedralres.ntdlr = atoi(iter->second.c_str());
    else if (iter->first == "CDLR")
      gin.dihedralres.cdlr = atof(iter->second.c_str());
    else if (iter->first == "PHILIN")
      gin.dihedralres.philin = atof(iter->second.c_str());

      //DISTANCERES
    else if (iter->first == "NTDIR")
      gin.distanceres.ntdir = atoi(iter->second.c_str());
    else if (iter->first == "NTDIRA")
      gin.distanceres.ntdira = atoi(iter->second.c_str());
    else if (iter->first == "CDIR")
      gin.distanceres.cdir = atof(iter->second.c_str());
    else if (iter->first == "DIR0")
      gin.distanceres.dir0 = atof(iter->second.c_str());
    else if (iter->first == "TAUDIR")
      gin.distanceres.taudir = atoi(iter->second.c_str());

      // ENERGYMIN
    else if (iter->first == "NTEM")
      gin.energymin.ntem = atoi(iter->second.c_str());
    else if (iter->first == "NCYC")
      gin.energymin.ncyc = atoi(iter->second.c_str());
    else if (iter->first == "DELE")
      gin.energymin.dele = atof(iter->second.c_str());
    else if (iter->first == "DX0")
      gin.energymin.dx0 = atof(iter->second.c_str());
    else if (iter->first == "DXM")
      gin.energymin.dxm = atof(iter->second.c_str());
    else if (iter->first == "NMIN")
      gin.energymin.nmin = atoi(iter->second.c_str());
    else if (iter->first == "FLIM")
      gin.energymin.flim = atof(iter->second.c_str());

      // EWARN
    else if (iter->first == "MAXENER")
      gin.ewarn.maxener = atof(iter->second.c_str());

      // FORCE
    else if (iter->first == "NTF[1]")
      gin.force.ntf[0] = atoi(iter->second.c_str());
    else if (iter->first == "NTF[2]")
      gin.force.ntf[1] = atoi(iter->second.c_str());
    else if (iter->first == "NTF[3]")
      gin.force.ntf[2] = atoi(iter->second.c_str());
    else if (iter->first == "NTF[4]")
      gin.force.ntf[3] = atoi(iter->second.c_str());
    else if (iter->first == "NTF[5]")
      gin.force.ntf[4] = atoi(iter->second.c_str());
    else if (iter->first == "NTF[6]")
      gin.force.ntf[5] = atoi(iter->second.c_str());
    else if (iter->first == "NTF[7]")
      gin.force.ntf[6] = atoi(iter->second.c_str());
    else if (iter->first == "NTF[8]")
      gin.force.ntf[7] = atoi(iter->second.c_str());
    else if (iter->first == "NTF[9]")
      gin.force.ntf[8] = atoi(iter->second.c_str());
    else if (iter->first == "NTF[10]")
      gin.force.ntf[9] = atoi(iter->second.c_str());

      // GEOMCONSTRAINTS
    else if (iter->first == "NTCPH")
      gin.geomconstraints.ntcph = atoi(iter->second.c_str());
    else if (iter->first == "NTCPN")
      gin.geomconstraints.ntcpn = atoi(iter->second.c_str());
    else if (iter->first == "NTCS")
      gin.geomconstraints.ntcs = atoi(iter->second.c_str());
    else if (iter->first == "SHKTOL")
      gin.geomconstraints.shktol = atof(iter->second.c_str());

      // GROMOS96COMPAT
    else if (iter->first == "NTNB96")
      gin.gromos96compat.ntnb96 = atoi(iter->second.c_str());
    else if (iter->first == "NTR96")
      gin.gromos96compat.ntr96 = atoi(iter->second.c_str());
    else if (iter->first == "NTP96")
      gin.gromos96compat.ntp96 = atoi(iter->second.c_str());
    else if (iter->first == "NTG96")
      gin.gromos96compat.ntg96 = atoi(iter->second.c_str());

      // INITIALISE
    else if (iter->first == "NTIVEL")
      gin.initialise.ntivel = atoi(iter->second.c_str());
    else if (iter->first == "NTISHK")
      gin.initialise.ntishk = atoi(iter->second.c_str());
    else if (iter->first == "NTINHT")
      gin.initialise.ntinht = atoi(iter->second.c_str());
    else if (iter->first == "NTINHB")
      gin.initialise.ntinhb = atoi(iter->second.c_str());
    else if (iter->first == "NTISHI")
      gin.initialise.ntishi = atoi(iter->second.c_str());
    else if (iter->first == "NTIRTC")
      gin.initialise.ntirtc = atoi(iter->second.c_str());
    else if (iter->first == "NTICOM")
      gin.initialise.nticom = atoi(iter->second.c_str());
    else if (iter->first == "NTISTI")
      gin.initialise.ntisti = atoi(iter->second.c_str());
    else if (iter->first == "IG")
      gin.initialise.ig = atoi(iter->second.c_str());
    else if (iter->first == "TEMPI")
      gin.initialise.tempi = atof(iter->second.c_str());

      // INNERLOOP
    else if (iter->first == "NTILM")
      gin.innerloop.ntilm = atoi(iter->second.c_str());
    else if (iter->first == "NTILS")
      gin.innerloop.ntils = atoi(iter->second.c_str());

      // INTEGRATE
    else if (iter->first == "NINT")
      gin.integrate.nint = atoi(iter->second.c_str());

      // JVALUERES
    else if (iter->first == "NTJVR")
      gin.jvalueres.ntjvr = atoi(iter->second.c_str());
    else if (iter->first == "NTJVRA")
      gin.jvalueres.ntjvra = atoi(iter->second.c_str());
    else if (iter->first == "CJVR")
      gin.jvalueres.cjvr = atof(iter->second.c_str());
    else if (iter->first == "TAUJVR")
      gin.jvalueres.taujvr = atof(iter->second.c_str());
    else if (iter->first == "LE")
      gin.jvalueres.le = atoi(iter->second.c_str());
    else if (iter->first == "NGRID")
      gin.jvalueres.ngrid = atoi(iter->second.c_str());
    else if (iter->first == "DELTA")
      gin.jvalueres.delta = atof(iter->second.c_str());
    else if (iter->first == "NTWJV")
      gin.jvalueres.write = atoi(iter->second.c_str());

      //LAMBDAS
    else if (iter->first.substr(0, 4) == "ALI[") {
      unsigned int i = atoi(iter->first.substr(4, iter->first.find("]")).c_str());
      if (i <= gin.lambdas.lambints.size())
        gin.lambdas.lambints[i - 1].ali = atof(iter->second.c_str());
      else
        printError(iter->first + " in joblist out of range");
    } else if (iter->first.substr(0, 4) == "BLI[") {
      unsigned int i = atoi(iter->first.substr(4, iter->first.find("]")).c_str());
      if (i <= gin.lambdas.lambints.size())
        gin.lambdas.lambints[i - 1].bli = atof(iter->second.c_str());
      else
        printError(iter->first + " in joblist out of range");
    } else if (iter->first.substr(0, 4) == "CLI[") {
      unsigned int i = atoi(iter->first.substr(4, iter->first.find("]")).c_str());
      if (i <= gin.lambdas.lambints.size())
        gin.lambdas.lambints[i - 1].cli = atof(iter->second.c_str());
      else
        printError(iter->first + " in joblist out of range");
    } else if (iter->first.substr(0, 4) == "DLI[") {
      unsigned int i = atoi(iter->first.substr(4, iter->first.find("]")).c_str());
      if (i <= gin.lambdas.lambints.size())
        gin.lambdas.lambints[i - 1].dli = atof(iter->second.c_str());
      else
        printError(iter->first + " in joblist out of range");
    } else if (iter->first.substr(0, 4) == "ELI[") {
      unsigned int i = atoi(iter->first.substr(4, iter->first.find("]")).c_str());
      if (i <= gin.lambdas.lambints.size())
        gin.lambdas.lambints[i - 1].eli = atof(iter->second.c_str());
      else
        printError(iter->first + " in joblist out of range");
    }

    // LOCALELEV
    else if (iter->first == "NTLES")
      gin.localelev.ntles = atoi(iter->second.c_str());
    else if (iter->first == "NLEPOT")
      gin.localelev.nlepot = atoi(iter->second.c_str());
    else if (iter->first == "NTLESA")
      gin.localelev.ntlesa = atoi(iter->second.c_str());
    else if (iter->first.substr(0,6) == "NLEPID[") {
      int i = atoi(iter->first.substr(6, iter->first.find("]")).c_str());
      if(i > gin.localelev.nlepot) {
        stringstream msg;
        msg << "NLEPID[" << i << "] cannot be set when NLEPOT = " << gin.localelev.nlepot
                << " (LOCALELEV block)";
        printError(msg.str());
      }
      int count = 0;
      map<int, int> tmp;
      for(map<int,int>::iterator it = gin.localelev.nlepid_ntlerf.begin();
      it != gin.localelev.nlepid_ntlerf.end(); ++it){
        if(count == i) {
          tmp.insert( pair<int, int>(atoi(iter->second.c_str()), it->second));
        } else {
          tmp.insert( pair<int, int>(it->first, it->second));
        }
        count++;
      }
      gin.localelev.nlepid_ntlerf = tmp;
      if(int(gin.localelev.nlepid_ntlerf.size()) != gin.localelev.nlepot) {
        printError("NLEPID in LOCALELEV block is ambiguous");
      }
    }
    else if (iter->first.substr(0,6) == "NTLEFR[") {
      int i = atoi(iter->first.substr(6, iter->first.find("]")).c_str());
      int count = 0;
      for(map<int,int>::iterator it = gin.localelev.nlepid_ntlerf.begin();
      it != gin.localelev.nlepid_ntlerf.end(); ++it){
        if(count == i) {
          it->second = atoi(iter->second.c_str());
          if (it->second < 0 || it->second >1) {
            printError("NTLEFR in LOCALELEV must be 0 or 1");
          }
          break;
        }
        count++;
      }
    }
    

      // MULTIBATH
    else if (iter->first == "ALGORITHM")
      gin.multibath.algorithm = atoi(iter->second.c_str());
    else if (iter->first.substr(0, 6) == "TEMP0[") {
      int i = atoi(iter->first.substr(6, iter->first.find("]")).c_str());
      if (i <= gin.multibath.nbaths)
        gin.multibath.temp0[i - 1] = atof(iter->second.c_str());
      else
        printError(iter->first + " in joblist out of range");
    } else if (iter->first.substr(0, 5) == "TAU[") {
      int i = atoi(iter->first.substr(5, iter->first.find("]")).c_str());
      if (i <= gin.multibath.nbaths)
        gin.multibath.tau[i - 1] = atof(iter->second.c_str());
      else
        printError(iter->first + " in joblist out of range");
    } // MULTICELL

    // MULTISTEP
    else if (iter->first == "STEPS")
      gin.multistep.steps = atoi(iter->second.c_str());
    else if (iter->first == "BOOST")
      gin.multistep.boost = atoi(iter->second.c_str());

      // NEIGHBOURLIST
    else if (iter->first == "NMPRPL")
      gin.neighbourlist.nmprpl = atoi(iter->second.c_str());
    else if (iter->first == "NUPRPL")
      gin.neighbourlist.nuprpl = atoi(iter->second.c_str());
    else if (iter->first == "RCPRL")
      gin.neighbourlist.rcprpl = atof(iter->second.c_str());
    else if (iter->first == "GRPRPL")
      gin.neighbourlist.grprpl = atof(iter->second.c_str());
    else if (iter->first == "NMTWPL")
      gin.neighbourlist.nmtwpl = atoi(iter->second.c_str());
    else if (iter->first == "NUTWPL")
      gin.neighbourlist.nutwpl = atoi(iter->second.c_str());
    else if (iter->first == "RSTWPL")
      gin.neighbourlist.rstwpl = atof(iter->second.c_str());
    else if (iter->first == "RLTWPL")
      gin.neighbourlist.rltwpl = atof(iter->second.c_str());
    else if (iter->first == "NUIRIN")
      gin.neighbourlist.nuirin = atoi(iter->second.c_str());
    else if (iter->first == "NUSRIN")
      gin.neighbourlist.nusrin = atoi(iter->second.c_str());
    else if (iter->first == "NMTWIN")
      gin.neighbourlist.nmtwin = atoi(iter->second.c_str());
    else if (iter->first == "RCTWIN")
      gin.neighbourlist.rctwin = atof(iter->second.c_str());
    else if (iter->first == "NCGCEN")
      gin.neighbourlist.ncgcen = atoi(iter->second.c_str());

      // NONBONDED
    else if (iter->first == "NLRELE")
      gin.nonbonded.nlrele = atoi(iter->second.c_str());
    else if (iter->first == "APPAK")
      gin.nonbonded.appak = atof(iter->second.c_str());
    else if (iter->first == "RCRF")
      gin.nonbonded.rcrf = atof(iter->second.c_str());
    else if (iter->first == "EPSRF")
      gin.nonbonded.epsrf = atof(iter->second.c_str());
    else if (iter->first == "NSHAPE")
      gin.nonbonded.nshape = atoi(iter->second.c_str());
    else if (iter->first == "ASHAPE")
      gin.nonbonded.ashape = atof(iter->second.c_str());
    else if (iter->first == "NA2CLC")
      gin.nonbonded.na2clc = atoi(iter->second.c_str());
    else if (iter->first == "TOLA2")
      gin.nonbonded.tola2 = atof(iter->second.c_str());
    else if (iter->first == "EPSLS")
      gin.nonbonded.epsls = atof(iter->second.c_str());
    else if (iter->first == "NKX")
      gin.nonbonded.nkx = atoi(iter->second.c_str());
    else if (iter->first == "NKY")
      gin.nonbonded.nky = atoi(iter->second.c_str());
    else if (iter->first == "NKZ")
      gin.nonbonded.nkz = atoi(iter->second.c_str());
    else if (iter->first == "KCUT")
      gin.nonbonded.kcut = atof(iter->second.c_str());
    else if (iter->first == "NGX")
      gin.nonbonded.ngx = atoi(iter->second.c_str());
    else if (iter->first == "NGY")
      gin.nonbonded.ngy = atoi(iter->second.c_str());
    else if (iter->first == "NGZ")
      gin.nonbonded.ngz = atoi(iter->second.c_str());
    else if (iter->first == "NASORD")
      gin.nonbonded.nasord = atoi(iter->second.c_str());
    else if (iter->first == "NFDORD")
      gin.nonbonded.nfdord = atoi(iter->second.c_str());
    else if (iter->first == "NALIAS")
      gin.nonbonded.nalias = atoi(iter->second.c_str());
    else if (iter->first == "NSPORD")
      gin.nonbonded.nspord = atoi(iter->second.c_str());
    else if (iter->first == "NQEVAL")
      gin.nonbonded.nqeval = atoi(iter->second.c_str());
    else if (iter->first == "FACCUR")
      gin.nonbonded.faccur = atof(iter->second.c_str());
    else if (iter->first == "NRDGRD")
      gin.nonbonded.nrdgrd = atoi(iter->second.c_str());
    else if (iter->first == "NWRGRD")
      gin.nonbonded.nwrgrd = atoi(iter->second.c_str());
    else if (iter->first == "NLRLJ")
      gin.nonbonded.nlrlj = atoi(iter->second.c_str());
    else if (iter->first == "SLVDNS")
      gin.nonbonded.slvdns = atoi(iter->second.c_str());

      // OVERALLTRANSROT
    else if (iter->first == "NCMTR")
      gin.overalltransrot.ncmtr = atoi(iter->second.c_str());
    else if (iter->first == "NCMRO")
      gin.overalltransrot.ncmro = atoi(iter->second.c_str());
    else if (iter->first == "CMAMX")
      gin.overalltransrot.cmamx = atoi(iter->second.c_str());
    else if (iter->first == "CMAMY")
      gin.overalltransrot.cmamy = atoi(iter->second.c_str());
    else if (iter->first == "CMAMZ")
      gin.overalltransrot.cmamz = atoi(iter->second.c_str());

      // PAIRLIST
    else if (iter->first == "NSNB")
      gin.pairlist.nsnb = atoi(iter->second.c_str());
    else if (iter->first == "RCUTP")
      gin.pairlist.rcutp = atof(iter->second.c_str());
    else if (iter->first == "RCUTL")
      gin.pairlist.rcutl = atof(iter->second.c_str());
    else if (iter->first == "TYPE")
      gin.pairlist.type = atoi(iter->second.c_str());

      // PATHINT
    else if (iter->first == "NTPI")
      gin.pathint.ntpi = atoi(iter->second.c_str());

      // PERSCALE
    else if (iter->first == "RESTYPE")
      gin.perscale.restype = atoi(iter->second.c_str());
    else if (iter->first == "KDIH")
      gin.perscale.kdih = atof(iter->second.c_str());
    else if (iter->first == "KJ")
      gin.perscale.kj = atof(iter->second.c_str());
    else if (iter->first == "DIFF")
      gin.perscale.diff = atof(iter->second.c_str());
    else if (iter->first == "RATIO")
      gin.perscale.ratio = atof(iter->second.c_str());
    else if (iter->first == "READ")
      gin.perscale.read = atoi(iter->second.c_str());

      // PERTURBATION
    else if (iter->first == "NTG")
      gin.perturbation.ntg = atoi(iter->second.c_str());
    else if (iter->first == "NRDGL")
      gin.perturbation.nrdgl = atoi(iter->second.c_str());
    else if (iter->first == "RLAM")
      gin.perturbation.rlam = atof(iter->second.c_str());
    else if (iter->first == "DLAMT")
      gin.perturbation.dlamt = atof(iter->second.c_str());
    else if (iter->first == "ALPHLJ")
      gin.perturbation.alphlj = atof(iter->second.c_str());
    else if (iter->first == "ALPHC")
      gin.perturbation.alphc = atof(iter->second.c_str());
    else if (iter->first == "NLAM")
      gin.perturbation.nlam = atoi(iter->second.c_str());
    else if (iter->first == "NSCALE")
      gin.perturbation.nscale = atoi(iter->second.c_str());

      // POLARISE
    else if (iter->first == "COS")
      gin.polarise.cos = atoi(iter->second.c_str());
    else if (iter->first == "EFIELD")
      gin.polarise.efield = atoi(iter->second.c_str());
    else if (iter->first == "MINFIELD")
      gin.polarise.minfield = atof(iter->second.c_str());
    else if (iter->first == "DAMP")
      gin.polarise.damp = atoi(iter->second.c_str());
    else if (iter->first == "WRITE")
      gin.polarise.write = atoi(iter->second.c_str());

      // POSITIONRES
    else if (iter->first == "NTPOR")
      gin.positionres.ntpor = atoi(iter->second.c_str());
    else if (iter->first == "NTPORB")
      gin.positionres.ntporb = atoi(iter->second.c_str());
    else if (iter->first == "NTPORS")
      gin.positionres.ntpors = atoi(iter->second.c_str());
    else if (iter->first == "CPOR")
      gin.positionres.cpor = atof(iter->second.c_str());

      // PRESSURESCALE
    else if (iter->first == "COUPLE")
      gin.pressurescale.couple = atoi(iter->second.c_str());
    else if (iter->first == "SCALE")
      gin.pressurescale.scale = atoi(iter->second.c_str());
    else if (gin.pressurescale.found && iter->first == "COMP")
      gin.pressurescale.comp = atof(iter->second.c_str());
    else if (iter->first == "TAUP")
      gin.pressurescale.taup = atof(iter->second.c_str());
    else if (iter->first == "VIRIAL")
      gin.pressurescale.virial = atoi(iter->second.c_str());
    else if (iter->first == "X_SEMI")
      gin.pressurescale.x_semi = atoi(iter->second.c_str());
    else if (iter->first == "Y_SEMI")
      gin.pressurescale.y_semi = atoi(iter->second.c_str());
    else if (iter->first == "Z_SEMI")
      gin.pressurescale.z_semi = atoi(iter->second.c_str());

      // PRINTOUT
    else if (iter->first == "NTPR")
      gin.printout.ntpr = atoi(iter->second.c_str());
    else if (iter->first == "NTPP")
      gin.printout.ntpp = atoi(iter->second.c_str());

      // RANDOMNUMBERS
    else if (iter->first == "NTRNG")
      gin.randomnumbers.ntrng = atoi(iter->second.c_str());
    else if (iter->first == "NTGSL")
      gin.randomnumbers.ntgsl = atoi(iter->second.c_str());

      // READTRAJ
    else if (iter->first == "NTRD")
      gin.readtraj.ntrd = atoi(iter->second.c_str());
    else if (iter->first == "NTRN")
      gin.readtraj.ntrn = atoi(iter->second.c_str());
    else if (iter->first == "NTRB")
      gin.readtraj.ntrb = atoi(iter->second.c_str());
    else if (iter->first == "NTSHK")
      gin.readtraj.ntshk = atoi(iter->second.c_str());

      // REPLICA
    else if (iter->first.substr(0, 4) == "RET[") {
      unsigned int i = atoi(iter->first.substr(4, iter->first.find("]")).c_str());
      if (i <= gin.replica.ret.size())
        gin.replica.ret[i - 1] = atof(iter->second.c_str());
      else
        printError(iter->first + " in joblist out of range");
    } else if (iter->first.substr(0, 6) == "RELAM[") {
      unsigned int i = atoi(iter->first.substr(6, iter->first.find("]")).c_str());
      if (i <= gin.replica.relam.size())
        gin.replica.relam[i - 1] = atof(iter->second.c_str());
      else
        printError(iter->first + " in joblist out of range");
    } else if (iter->first.substr(0, 5) == "RETS[") {
      unsigned int i = atoi(iter->first.substr(6, iter->first.find("]")).c_str());
      if (i <= gin.replica.rets.size())
        gin.replica.rets[i - 1] = atof(iter->second.c_str());
      else
        printError(iter->first + " in joblist out of range");
    } else if (iter->first == "LRESCALE")
      gin.replica.lrescale = atoi(iter->second.c_str());
    else if (iter->first == "NRETRIAL")
      gin.replica.nretrial = atoi(iter->second.c_str());
    else if (iter->first == "NREQUIL")
      gin.replica.nrequil = atoi(iter->second.c_str());
    else if (iter->first == "NREJOB")
      gin.replica.nrejob = atoi(iter->second.c_str());
    else if (iter->first == "NREWRT")
      gin.replica.nrewrt = atoi(iter->second.c_str());

      // ROTTRANS
    else if (iter->first == "RTC")
      gin.rottrans.rtc = atoi(iter->second.c_str());

      // STEP
    else if (iter->first == "NSTLIM")
      gin.step.nstlim = atoi(iter->second.c_str());
    else if (iter->first == "T")
      gin.step.t = atof(iter->second.c_str());
    else if (iter->first == "DT")
      gin.step.dt = atof(iter->second.c_str());

      // STOCHDYN
    else if (iter->first == "NTSD")
      gin.stochdyn.ntsd = atoi(iter->second.c_str());
    else if (iter->first == "NTFR")
      gin.stochdyn.ntfr = atoi(iter->second.c_str());
    else if (iter->first == "NSFR")
      gin.stochdyn.nsfr = atoi(iter->second.c_str());
    else if (iter->first == "NBREF")
      gin.stochdyn.nbref = atoi(iter->second.c_str());
    else if (iter->first == "RCUTF")
      gin.stochdyn.rcutf = atof(iter->second.c_str());
    else if (iter->first == "CFRIC")
      gin.stochdyn.cfric = atof(iter->second.c_str());
    else if (iter->first == "TEMPSD")
      gin.stochdyn.tempsd = atof(iter->second.c_str());


      // SYSTEM
    else if (iter->first == "NPM")
      gin.system.npm = atoi(iter->second.c_str());
    else if (iter->first == "NSM")
      gin.system.nsm = atoi(iter->second.c_str());

      // THERMOSTAT
    else if (iter->first == "NTT")
      gin.thermostat.ntt = atoi(iter->second.c_str());
    else if (iter->first.substr(0, 7) == "TEMBTH[") {
      unsigned int i = atoi(iter->first.substr(7, iter->first.find("]")).c_str());
      if (i <= gin.thermostat.baths.size())
        gin.thermostat.baths[i - 1].tembth = atof(iter->second.c_str());
      else
        printError(iter->first + " in joblist is out of range");
    }// I don't really dare to let the use change other things...

      // UMBRELLA
    else if (iter->first == "NTUS")
      gin.umbrella.ntus = atoi(iter->second.c_str());
    else if (iter->first == "USCST1")
      gin.umbrella.uscst1 = atof(iter->second.c_str());
    else if (iter->first == "USCST2")
      gin.umbrella.uscst2 = atof(iter->second.c_str());
    else if (iter->first == "USREF1")
      gin.umbrella.usref1 = atof(iter->second.c_str());
    else if (iter->first == "USREF2")
      gin.umbrella.usref2 = atof(iter->second.c_str());

      // VIRIAL
    else if (iter->first == "NTV")
      gin.virial.ntv = atoi(iter->second.c_str());
    else if (iter->first == "NTVG")
      gin.virial.ntvg = atoi(iter->second.c_str());

      //WRITETRAJ
    else if (iter->first == "NTWX")
      gin.writetraj.ntwx = atoi(iter->second.c_str());
    else if (iter->first == "NTWSE")
      gin.writetraj.ntwse = atoi(iter->second.c_str());
    else if (iter->first == "NTWV")
      gin.writetraj.ntwv = atoi(iter->second.c_str());
    else if (iter->first == "NTWF")
      gin.writetraj.ntwf = atoi(iter->second.c_str());
    else if (iter->first == "NTWE")
      gin.writetraj.ntwe = atoi(iter->second.c_str());
    else if (iter->first == "NTWG")
      gin.writetraj.ntwe = atoi(iter->second.c_str());
    else if (iter->first == "NTWB")
      gin.writetraj.ntwb = atoi(iter->second.c_str());


      // XRAYRES
    else if (iter->first == "NTXR")
      gin.xrayres.ntxr = atoi(iter->second.c_str());
    else if (iter->first == "NTXLE")
      gin.xrayres.ntxle = atoi(iter->second.c_str());
    else if (iter->first == "CXR")
      gin.xrayres.cxr = atof(iter->second.c_str());
    else if (iter->first == "NTWXR")
      gin.xrayres.ntwxr = atoi(iter->second.c_str());
    else if (iter->first == "NTWDE")
      gin.xrayres.ntwde = atoi(iter->second.c_str());
    else if (iter->first == "NTWX<")
      gin.xrayres.ntwxm = atoi(iter->second.c_str());
    else if (iter->first == "CXTAU")
      gin.xrayres.cxtau = atof(iter->second.c_str());
    else if (iter->first == "RDAVG")
      gin.xrayres.rdavg = atoi(iter->second.c_str());

    else
    throw gromos::Exception("mk_script", "Cannot automatically change "
            + iter->first + " in input file");
  }
}

