/**
 * @file fit_ener.cc
 * Recalculates interaction energies for a solute molecule superposed on a
 * specific location
 */

/**
 * @page programs Program Documentation
 *
 * @anchor fit_ener
 * @section fit_ener Recalculates interaction energies for a solute molecule 
 * superposed on a specific location
 * @author @ref co
 * @date 16-07-2015
 *
 * THIS PROGRAM IS ONLY INTENDED FOR VERY SPECIFIC USE:
 * IT WILL TAKE THE COORDINATES OF THE SOLUTE SPECIFIED BY @FITTOPO AND ADD
 * THE SOLVENT FROM THE TRAJECTORY. THEN IT CALCULATES THE ENERGY OF THE 
 * SELECTED ATOMS
 *
 * Program ener can recalculate interaction energies over molecular trajectory 
 * files using the interaction parameters specified in the molecular topology 
 * file.
 *
 * Nonbonded interactions are calculated for all selected atoms with all other
 * atoms in the system. Some atoms can be specified as being soft, indicating 
 * that interactions involving any of these atoms have a specified softness
 * parameter, for all other atoms in the system, the softness parameter
 * @f$\alpha = 0@f$. Vanderwaals interactions between particles i and j are
 * calculated as
 * 
 * @f[ V^{vdw}_{ij}=\left[\frac{C_{12}(i,j)}{ ( r_{ij}^6 + \alpha_{LJ} \lambda ^2 C_{126})}-C_6(i,j)\right] \frac{1}{(r_{ij}^6 + \alpha_{LJ} \lambda ^2 C_{126})} @f]
 *
 * with @f$C_{126} = C_{12}/C_6 @f$ for @f$C_{12}@f$ and @f$C_6@f$ unequal 0,
 * @f$C_{126} = 0@f$ otherwise. @f$C_{12}@f$ and @f$C_6@f$ are the interaction
 * parameters taken from the topology, @f$\lambda@f$ and @f$\alpha_{LJ}@f$ are
 * specified by the user. Similarly, the electrostatic interaction, including
 * reaction field contribution for a homogeneous medium outside the cutoff
 * sphere is calculated as 
 *
 * @f[ V^{ele}_{ij}=\frac{q_iq_j}{4\pi\epsilon_0}\left[\frac{1}{(r^2_{ij}+\alpha_{C}\lambda^2)^{1/2}} - \frac{\frac{1}{2}C_{RF}r_{ij}^2}{(R_{RF}^2+\alpha_{C}\lambda^2)^{3/2}} - \frac{(1-\frac{1}{2}C_{RF})}{R_{RF}}\right] @f]
 *
 * where @f$\epsilon_0@f$ is the dielectric permittivity of vacuum and 
 * @f$q_i@f$ and @f$q_j@f$ are the atomic partial charges. @f$R_{RF}@f$ is the
 * reaction field cutoff distance, here assumed to be the same as the
 * interaction cutoff. @f$\alpha_{C}@f$ and @f$\lambda@f$ are again user 
 * specified. @f$C_{RF}@f$ is calculated from the reaction field dielectric
 * constant @f$\epsilon_{RF}@f$ and @f$\kappa_{RF}^{-1}@f$ (user specified) as
 *
 * @f[ C_{RF} = \frac{ (2 - 2 \epsilon_{RF}) (1 + \kappa_{RF}^{-1} R_{RF}) - \epsilon_{RF} (\kappa_{RF}^{-1} R_{RF})^2 }{ (1 + 2 \epsilon_{RF}) (1 + \kappa_{RF}^{-1} R_{RF}) + \epsilon_{RF} (\kappa_{RF}^{-1} R_{RF})^2 } @f]
 *
 * The bonded interactiona are calculated for all specified properties using 
 * the following interaction functions. For bonds we use:
 *
 * @f[ V^{bond}=\frac{1}{4}K_{b_n}\left[b_n^2 - b_{0_n}^2\right]^2 @f]
 *
 * with @f$b_n@f$ the actual bond length, @f$K_{b_n}@f$ and @f$b_{0_n}@f$ the 
 * force constant and optimal bond length, respectively. For angles we use:
 *
 * @f[ V^{angle}=\frac{1}{2}K_{\theta_n}\left[\cos{\theta_n} - \cos{\theta_{0_n}}\right]^2 @f]
 *
 * with @f$\theta_n@f$ the actual bond angle, @f$K_{\theta_n}@f$ and 
 * @f$\theta_{0_n}@f$ the force constant and optimal bond angle respectively.
 * For proper torsional dihedral angle terms we use:
 *
 * @f[ V^{trig}=K_{\phi_n}\left[1+\cos(\delta_n)\cos(m_n\phi_n)\right] @f]
 *
 * with @f$\phi_n@f$ the actual dihedral angle value, @f$K_{\phi_n}@f$ the
 * force constant and @f$\delta_n@f$ and @f$m_n@f$ the phase shift and
 * multiplicity, respectively. Improper dihedral energy contributions are 
 * calculated from:
 * @f[ V^{har}=\frac{1}{2}K_{\xi_n}\left[\xi_n - \xi_{0_n}\right]^2 @f]
 *
 * with @f$\xi_n@f$ the actual dihedral angle value, @f$K_{\xi_n}@f$ and
 * @f$\xi_{0_n}@f$ are the force constant and optimal improper dihedral angle 
 * value.
 *
 * The programme prints out the total bonded and nonbonded energies separately,
 * as well as the overall total energy. It is easily modified to print out more
 * detailed energy contributions as well.
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;boundary type&gt; [&lt;gather method&gt;] </td></tr>
 * <tr><td> \@atoms</td><td>&lt;@ref AtomSpecifier "atoms" for nonbonded interaction&gt; </td></tr>
 * <tr><td> \@props</td><td>&lt;@ref PropertyContainer "properties" to be calculated&gt; </td></tr>
 * <tr><td> [\@time</td><td>&lt;@ref utils::Time "time and dt"&gt;] </td></tr>
 * <tr><td> \@cut</td><td>&lt;cut-off distance&gt; </td></tr>
 * <tr><td> \@eps</td><td>&lt;epsilon for reaction field contribution&gt; </td></tr>
 * <tr><td> \@kap</td><td>&lt;kappa_{RF}^{-1} for reaction field contribution&gt; </td></tr>
 * <tr><td> [\@RFex</td><td>&lt;calculate RF contribution for excluded atoms: on/off&gt;] </td></tr>
 * <tr><td> \@soft</td><td>&lt;soft @ref AtomSpecifier "atoms"&gt; </td></tr>
 * <tr><td> \@softpar</td><td>&lt;lam&gt; &lt;a_lj&gt; &lt;a_crf&gt; </td></tr>
 * <tr><td> \@traj</td><td>&lt;trajectory files&gt; </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  ener
    @topo    ex.top
    @pbc     r
    @atoms   1:3-13
    @props   d%1:1,2 a%1:1,2,3 t%1:1,2,4,6 t%1:4,2,5,6
    [@time    0 0.2]
    @cut     1.4
    @eps     61
    @kap     0.0
    @RFex    on
    @soft    1:4
    @softpar 0.5 1.51 0.5
    @fittopo  solute.top
    @fitcoord solute.cnf
    @fitatoms <atomspecifier> <atomspecifier>
    @traj    ex.tr
 @endverbatim
 *
 * <hr>
 */
#include <cassert>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include "../src/gmath/Vec.h"
#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/gio/InG96.h"
#include "../src/args/OutformatParser.h"
#include "../src/gcore/System.h"
#include "../src/gcore/AtomPair.h"
#include "../src/gcore/GromosForceField.h"
#include "../src/gcore/Box.h"
#include "../src/gio/InTopology.h"
#include "../src/bound/Boundary.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/utils/PropertyContainer.h"
#include "../src/utils/Energy.h"
#include "../src/utils/groTime.h"
#include "../src/fit/PositionUtils.h"
#include "../src/fit/RotationalFit.h"
#include "../src/gio/OutPdb.h"


using namespace std;
using namespace gcore;
using namespace gio;
using namespace bound;
using namespace args;
using namespace utils;
using fit::RotationalFit;

  
int main(int argc, char **argv){

  Argument_List knowns; 
  knowns << "topo" << "pbc" << "atoms" << "props" << "time" << "cut"
         << "eps" << "kap" << "soft" << "softpar" << "RFex" << "traj"
	 << "fittopo" << "fitcoord" << "fitatoms" << "calctopo" << "outformat"
	 << "combsys";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo     <molecular topology file>\n";
  usage += "\t@pbc      <boundary type> [<gather method>]\n";
  usage += "\t@atoms    <atoms for nonbonded (calcsys)>\n";
  usage += "\t@props    <propertyspecifier>\n";
  usage += "\t[@time    <time and dt>]\n";
  usage += "\t@cut      <cut-off distance>\n";
  usage += "\t@eps      <epsilon for reaction field correction>\n";
  usage += "\t@kap      <kappa_{RF}^{-1} for reaction field correction>\n";
  usage += "\t[@RFex    <calculate RF for excluded atoms: on/off>]\n";
  usage += "\t[@soft    <soft atoms>]\n";
  usage += "\t[@softpar <lam> <a_lj> <a_c>]\n";
  usage += "\t@fittopo  <topology for an amino acid>\n";
  usage += "\t@fitcoord <coordinate file for the amino acid>\n";
  usage += "\t@fitatoms <atomspecifier for topo> <atomspecifier for fittopo>\n";
  usage += "\t@calctopo  <topology for chimeric molecule>\n";
  usage += "\t[@outformat  <specify if coordinates should be written out>]\n";
  usage += "\t[@combsys  <combine sys and calcsys in output>]\n";
  usage += "\t@traj    <trajectory files>\n";
  
 
try{
  Arguments args(argc, argv, knowns, usage);

  cerr << "THIS PROGRAM IS ONLY INTENDED FOR VERY SPECIFIC USE:\n"
       << "It will take the coordinates of the amino acid, remove the first \n"
       << "two atoms assuming they are H1 and H2 from the NH3+, and do a \n"
       << "rotational fit onto the trajectory coordinates based on the two sets \n"
       << "of atoms in @fitatoms, afterwards translating to move the CAs of the \n"
       << "fitgroups on top of each other (-> it makes sense to include the CAs in the fitatoms).\n" 
       << "Then a new system is created by combining \n"
       << "the side chain atoms of the amino acid with everything but this sidechain \n"
       << "of the trajectory\n";
  
  // get the @time argument
  utils::Time time(args);

  //  read topologies
  InTopology it(args["topo"]);
  InTopology fitit(args["fittopo"]);
  InTopology calcit(args["calctopo"]);
  
  System sys(it.system());
  System fitsys(fitit.system());
  System calcsys(calcit.system());
  
  if(it.forceField().ForceField() != fitit.forceField().ForceField()){
    throw gromos::Exception("fit_ener", 
			    "force fields for topo and fit_topo are not the same\n");
  }
  
  
  if(it.forceField().ForceField() != calcit.forceField().ForceField()){
    throw gromos::Exception("fit_ener", 
			    "force fields for topo and calc_topo are not the same\n");
  }
    
  GromosForceField gff(it.forceField());

  // parse boundary conditions
  Boundary *pbc = BoundaryParser::boundary(calcsys, args);
  Boundary *pbcsys = BoundaryParser::boundary(sys, args);
  Boundary *pbcfit = BoundaryParser::boundary(fitsys, args);
  //parse gather method
  Boundary::MemPtr gathmethod = args::GatherParser::parse(sys, sys, args);

  // declare the energy class
  Energy en(calcsys, gff, *pbc);

  //  set atoms
  AtomSpecifier atoms(calcsys);
  {
    Arguments::const_iterator iter=args.lower_bound("atoms");
    Arguments::const_iterator to=args.upper_bound("atoms");
    for(;iter!=to;iter++){
      string spec=iter->second.c_str();
      atoms.addSpecifier(spec);
    }
  }
  en.setAtoms(atoms);

  // get fitatoms
  AtomSpecifier fit_atoms_sys(sys);
  AtomSpecifier fit_atoms_fit(fitsys);
  {
    Arguments::const_iterator iter=args.lower_bound("fitatoms");
    Arguments::const_iterator to=args.upper_bound("fitatoms");
    string spec=iter->second.c_str();
    fit_atoms_sys.addSpecifier(spec);
    iter++;
    spec=iter->second.c_str();
    fit_atoms_fit.addSpecifier(spec);
    iter++;
    if(iter!=to){
      throw gromos::Exception("fit_ener", "you have to give two arguments for fitatoms");
    }
    if(fit_atoms_sys.size() != fit_atoms_fit.size()){
      throw gromos::Exception("fit_ener", "fitatom groups for trajectory and fittopo have to have same size!");
      
    }   
  }

  // read in coordinates 
  {
    InG96 fitic;
    if (args.count("fitcoord") > 0) {
      fitic.open(args["fitcoord"].c_str());
      fitic >> fitsys;
      (*pbcfit.*gathmethod)();
    } else {
      throw gromos::Exception("fit_ener", "no fit-coordinates given");
    }
  }
  
    
  // RF for excluded atoms?
  {
    std::string s=args.getValue<string>("RFex",1,"on");
    if(s=="off")
      en.setRFexclusions(false);
    else
      en.setRFexclusions(true);
  }

  // set properties
  PropertyContainer props(calcsys, pbc);
  {
    Arguments::const_iterator iter=args.lower_bound("props");
    Arguments::const_iterator to=args.upper_bound("props");
    for(;iter!=to;iter++){
      string p=iter->second.c_str();
      props.addSpecifier(p);
    }
  }
  en.setProperties(props);

  // set non-bonded parameters
  //   get cut-off distance
  {
    double cut = args.getValue<double>("cut", false, 1.4);
    en.setCutOff(cut);
  }
  //  get epsilon and kappa_{RF}^{-1}
  {
    double eps = args.getValue<double>("eps", false, 1.0);
    double kap = args.getValue<double>("kap", false, 0.0);
    en.setRF(eps, kap);
  }
  // get soft atom list
  AtomSpecifier soft(calcsys);
  {
    bool lsoft=false;
    Arguments::const_iterator iter=args.lower_bound("soft");
    Arguments::const_iterator to=args.upper_bound("soft");
    for(;iter!=to;iter++){
      string spec=iter->second.c_str();
      soft.addSpecifier(spec);
      lsoft=true;
    }
    //  get softness parameter
    std::vector<double> softpar = args.getValues<double>("softpar", 3, lsoft, 
            Arguments::Default<double>() << 0.0 << 0.0 << 0.0);
    if (lsoft)
      en.setSoft(soft, softpar[0], softpar[1], softpar[2]);
  }
  
  // see if we want to write coordinates
  bool write_coord=false;
  bool combsys = false;
  if (args.count("combsys") >= 0) {
        combsys = true;
  }
  
  // prepare for coordinate writing
  string ext;
  // original system
  OutCoordinates *oref = OutformatParser::parse(args, ext);
  ostringstream pdbName;
  pdbName << "osref" << ext;
  ofstream osref;
  
  //  system with replaced amino acid
  OutCoordinates *ocalc = OutformatParser::parse(args, ext);
  ostringstream pdbName3;
  pdbName3 << "oscalc"  << ext;
  ofstream oscalc;
  
  if (args.count("outformat") >= 0) {
    write_coord=true;
          osref.open(pdbName.str().c_str());
          oref->open(osref);
          oref->select("SOLUTE");
          oref->writeTitle(pdbName.str());
  if (!combsys) {
          oscalc.open(pdbName3.str().c_str());
          ocalc->open(oscalc);
          ocalc->select("SOLUTE");
          ocalc->writeTitle(pdbName3.str());
    }
  }
  
  // reference system used only for final rotational fit for outcoordinates
  System refsys(it.system()); 
  if (args.count("traj") > 0) {
    InG96 ic;
    ic.open(args.lower_bound("traj")->second);
    ic.select("ALL");
    ic >> refsys;
    ic.close();
  } else {
     throw gromos::Exception("fit_ener", "no trajectory specified (@traj)");
  }
  AtomSpecifier refbb;
  refbb.setSystem(refsys);
  std::string s="a:C,N,CA,O";
  refbb.addSpecifier(s);
  RotationalFit refrf(refbb);
 
  // define input coordinate
  InG96 ic;

  // print titles
  cout << "# Time"
       << "              covalent"
       << "             nonbonded"
       << "                 Total"
       << endl;

  // declare some variables for averaging
  int num_frames=0;
  double cov=0.0;
  double nb=0.0;
  double tot=0.0;
  
  // atom specifiers for the backbone atoms of the residue
  // that is going to be replaced in sys and calcsys
  stringstream ss;
  ss << fit_atoms_sys.mol(0)+1 ; 
  ss << ":res("<<fit_atoms_sys.resnum(0)+1 << ":N,H,CA,C,O)";
  AtomSpecifier backbonesys(sys);
  AtomSpecifier backbonecalcsys(calcsys);
  backbonesys.addSpecifier(ss.str());
  backbonecalcsys.addSpecifier(ss.str());
  ss.str ("");
  
  // getting the number of atoms of the fit residue in sys - is there a less ugly way?
  AtomSpecifier fitresidue(sys);
  ss << fit_atoms_sys.mol(0)+1 << ":res("<<fit_atoms_sys.resnum(0)+1 << ":a)";
  fitresidue.addSpecifier(ss.str());
  int nresatoms=fitresidue.size();
  ss.str ("");
  
  
  
  // loop over all trajectories
  for(Arguments::const_iterator 
        iter=args.lower_bound("traj"),
        to=args.upper_bound("traj");
      iter!=to; ++iter){

      // open file
    ic.open((iter->second).c_str());
    ic.select("ALL");
    
      // loop over single trajectory
    while(!ic.eof()){
      ic >> sys >> time;    
  
      (*pbcsys.*gathmethod)();
      
      RotationalFit rf(fit_atoms_sys);
      rf.fit(fit_atoms_sys,fit_atoms_fit);

      // move fitsys so that its Ca is exactly on the Ca of the sys
      Vec dist(0.0, 0.0, 0.0);
      AtomSpecifier ca_sys(sys);
      AtomSpecifier ca_fitsys(fitsys);
      ss << fit_atoms_sys.mol(0)+1 << ":res("<<fit_atoms_sys.resnum(0)+1 << ":CA)";
      ca_sys.addSpecifier(ss.str());
      ss.str ("");
      ss << fit_atoms_fit.mol(0)+1 << ":res("<<fit_atoms_fit.resnum(0)+1 << ":CA)";
      ca_fitsys.addSpecifier(ss.str());
      dist=ca_sys.pos(0)-ca_fitsys.pos(0);
      ss.str ("");
      
      fit::PositionUtils::translate(&fitsys, dist);

      // calcsys gets solute from sys, only one residue from fitsys;
      // which residue is replaced is determined by the mol- and residuenum
      // of the first atom in the fit_atoms selection
      for(int m=0; m< calcsys.numMolecules(); m++) {
      int fitsyscnt=2; // starting at 2 because we skip the two additional Hs of the NH3+ group
      int syscnt=0;
	    for(int a=0; a < calcsys.mol(m).numAtoms(); a++) {
     
	      if (calcsys.mol(m).topology().resNum(a) == fit_atoms_sys.resnum(0) && m == fit_atoms_sys.mol(0)) {
	        calcsys.mol(m).pos(a) = fitsys.mol(fit_atoms_fit.mol(0)).pos(fitsyscnt);
	        fitsyscnt+=1;
	        syscnt=nresatoms-(fitsyscnt-2);
	      }
	      else {
	        calcsys.mol(m).pos(a) = sys.mol(m).pos(a+syscnt);
	      }
	    }
	  }
	  // now replace backbone positions of the one exchanged residue with the positions from sys
	  for (int a=0; a<backbonecalcsys.size(); a++) {
	    calcsys.mol(backbonecalcsys.mol(a)).pos(backbonecalcsys.atom(a)) = backbonesys.pos(a);
	  }
	  
      // calcsys gets the solvent and the box from sys
      calcsys.sol(0) = sys.sol(0);
      calcsys.box() = sys.box();
      calcsys.hasBox = true;
      
      // we have to gather with any method to get covalent interactions 
      // and charge-groups connected
      //pbc->gathergr();

      // calculate the energies
      en.calc();

      // print any ouput you like
      cout.precision(10);
      cout.setf(ios::right, ios::adjustfield);
      cout << time
	       << ' ' << setw(21) << en.cov()
           << ' ' << setw(21) << en.nb()
           << ' ' << setw(21) << en.tot()
	       << endl;

      //store some averages
      cov+=en.cov();
      nb+=en.nb();
      tot+=en.tot();

      // write coordinates if required
      if (write_coord) { 
          // do rotational fit on backbone of the reference
          // to make it easier to look at if molecules are moving
          // around between frames
          AtomSpecifier calcbb(calcsys);
          AtomSpecifier sysbb(sys);
          ss << "a:C,N,CA,O";
          calcbb.addSpecifier(ss.str());
          sysbb.addSpecifier(ss.str());
          ss.str ("");
          refrf.fit(refbb,calcbb);
          refrf.fit(refbb,sysbb);
          
          System comsys(sys);
          if (combsys) {
            for (int m=0; m<calcsys.numMolecules(); m++)
              // very rough: just put it on top, all coordinates will be duplicated
              comsys.addMolecule(calcsys.mol(m));
          } else {          
          ocalc->writeTimestep(time.steps(), time.time());
          *ocalc << calcsys; 
          }          

          oref->writeTimestep(time.steps(), time.time());
          *oref << comsys; 
      }      
      num_frames++;
    }
  }
      if (write_coord) { 
        if (!combsys) oscalc.close();
        osref.close();
  }
  // print out averages
  if(num_frames>1){
    cout.precision(10);
    cout.setf(ios::right, ios::adjustfield);
    cout << endl << "# ave."
         << ' ' << setw(21) << cov/num_frames 
         << ' ' << setw(21) << nb/num_frames
         << ' ' << setw(21) << tot/num_frames
         << endl;
  }
}
 
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}

