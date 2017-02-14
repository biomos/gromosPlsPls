/**
 * @file fit_ener.cc
 * Recalculates interaction energies for a solute molecule superposed on a
 * specific location
 */

/**
 * @page programs Program Documentation
 *
 * @anchor fit_ener
 * @section ener Recalculates interaction energies for a solute molecule 
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
#include <time.h>
#include <cassert>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include "../src/gmath/Vec.h"
#include "../src/gmath/Matrix.h"
#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
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
#include "../src/gio/OutPdb.h"

//#ifdef OMP
//#include <omp.h>
//#endif

using namespace std;
using namespace gcore;
using namespace gio;
using namespace bound;
using namespace args;
using namespace utils;

  
vector<int> max_dist_atoms(AtomSpecifier &as);

int main(int argc, char **argv){
  clock_t t;
  t=clock();

  Argument_List knowns; 
  knowns << "topo" << "pbc" << "atoms" << "props" << "time" << "cut"
         << "eps" << "kap" << "soft" << "softpar" << "RFex" << "traj"
     << "fittopo" << "fitcoord" << "fitatoms" << "outformat" << "alignaxes"
     << "rotate";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo     <molecular topology file>\n";
  usage += "\t@pbc      <boundary type> [<gather method>]\n";
  usage += "\t@atoms    <atoms for nonbonded>\n";
  usage += "\t@props    <propertyspecifier>\n";
  usage += "\t[@time    <time and dt>]\n";
  usage += "\t@cut      <cut-off distance>\n";
  usage += "\t@eps      <epsilon for reaction field correction>\n";
  usage += "\t@kap      <kappa_{RF}^{-1} for reaction field correction>\n";
  usage += "\t[@RFex    <calculate RF for excluded atoms: on/off>]\n";
  usage += "\t[@soft    <soft atoms>]\n";
  usage += "\t[@softpar <lam> <a_lj> <a_c>]\n";
  usage += "\t[@fittopo  <topology for molecule to superpose>]\n";
  usage += "\t[@fitcoord <coordinate file for molecule to superpose>]\n";
  usage += "\t[@fitatoms <atomspecifier for topo> <atomspecifier for fittopo>]\n";
  usage += "\t[@outformat <write coordinates of fitted system in specified format; slows down a lot!!>]\n";
  usage += "\t[@alignaxes <align longest axes of the two fitatom sets>]\n";
  usage += "\t[@rotate <number of steps for 360deg rotation along longest axis> [<random>]]\n";
  usage += "\t@traj    <trajectory files>\n";
  
 
try{
  Arguments args(argc, argv, knowns, usage);

  cerr << "THIS PROGRAM IS ONLY INTENDED FOR VERY SPECIFIC USE:\n"
       << "  IT WILL TAKE THE COORDINATES OF THE SOLUTE SPECIFIED BY @FITTOPO AND ADD\n"
       << "  THE SOLVENT FROM THE TRAJECTORY. THEN IT CALCULATES THE ENERGY OF THE\n"
       << "  SELECTED ATOMS\n";
  
  // get the @time argument
  utils::Time time(args);

  //  read topologies
  InTopology it(args["topo"]);
  InTopology fitit(args["fittopo"]);
  
  System sys(it.system());
  System fitsys(fitit.system());
  System calcsys(fitit.system());
  
  if(it.forceField().ForceField() != fitit.forceField().ForceField()){
    throw gromos::Exception("fit_ener", 
                "force fields for topo and fit_topo are not the same\n");
  }
    
  GromosForceField gff(it.forceField());

  // parse boundary conditions
  Boundary *pbc = BoundaryParser::boundary(calcsys, args);
  Boundary *pbcfit = BoundaryParser::boundary(fitsys, args);
  Boundary *pbcsys = BoundaryParser::boundary(sys, args);

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
    if(fit_atoms_sys.size() < 2){
      throw gromos::Exception("fit_ener", "give at least 2 atoms each for fitatoms (trajectory)");
      
    }    
    if(fit_atoms_fit.size() < 2){
      throw gromos::Exception("fit_ener", "give at least 2 atoms each for fitatoms (fittopo)");
    }
  }

  // read in coordinates 
  {
    InG96 fitic;
    fitic.open(args["fitcoord"].c_str());
    fitic >> fitsys;
    pbcfit->gathergr();
  }
  vector<int> maxdistatoms_fit=max_dist_atoms(fit_atoms_fit);
  // cerr << "maxdistatoms_fit " << maxdistatoms_fit[0] << " " << maxdistatoms_fit[1] << endl;
    
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
  
  // prepare for coordinate writing
  string ext;
  OutCoordinates *ofit = OutformatParser::parse(args, ext);
  ostringstream pdbName;
  pdbName << "fitsys" << ext;
  ofstream osfit;
  OutCoordinates *oref = OutformatParser::parse(args, ext);
  ostringstream pdbNamer;
  pdbNamer << "refsys" << ext;
  ofstream osref;
  if (args.count("outformat") >= 0) {
    write_coord=true;
    osfit.open(pdbName.str().c_str());
    osref.open(pdbNamer.str().c_str());
    
          ofit->open(osfit);
          ofit->select("ALL");
          ofit->writeTitle(pdbName.str());
          oref->open(osref);
          oref->select("ALL");
          oref->writeTitle(pdbNamer.str());
  }
  
  bool alignaxes=false;
  if (args.count("alignaxes") >= 0) alignaxes=true;
  
  int rotsteps=1;
  bool randomrot=false;
  if (args.count("rotate") > 0) {
    if (!alignaxes)
      cerr << "\nWARNING: rotating around the axis without aligning along it first"
           << " -- this is probably not what you want to do, you might want to \n"
           << " you might want to add @alignaxes !!\n\n";
    
    rotsteps = args.getValue<int>("rotate", false, 1);
    
    Arguments::const_iterator iter=args.lower_bound("rotate");
    Arguments::const_iterator to=args.upper_bound("rotate");
    iter++;
    if (iter != to) {
        string r=iter->second.c_str();
        if (r == "random") randomrot = true;
        else cerr << "\nWARNING: if a second argument is given for @rotate it should be 'random',\n"
                  << "         otherwise you are not using random rotation!!";
    }
  }
 
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
  
  //  #ifdef OMP
  //  omp_set_num_threads(4); //set the number of cpus for the parallel section
  //  #endif // OMP
  
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
      pbcsys->gathergr();
      
      // calcsys gets the solvent and the box from sys
      calcsys.sol(0) = sys.sol(0);
      calcsys.box() = sys.box();
      calcsys.hasBox = true;
      
      
      //calculate norm and angle of the selections' longest axes
      vector<int> maxdistatoms_sys=max_dist_atoms(fit_atoms_sys);
      Vec v1=fit_atoms_fit.pos(maxdistatoms_fit[1])-fit_atoms_fit.pos(maxdistatoms_fit[0]);
      Vec v2=fit_atoms_sys.pos(maxdistatoms_sys[1])-fit_atoms_sys.pos(maxdistatoms_sys[0]);
      Vec cross=v1.cross(v2);
      double angle=acos((v1.dot(v2)) / (v1.abs() * v2.abs()))*180 / M_PI;
          
      // rotate fitsys to align longest axes 
      if (alignaxes) {           
        gmath::Matrix mat=fit::PositionUtils::rotateAround(cross, angle);
        fit::PositionUtils::rotate(&fitsys, mat);
      }
      
      // create pairlist only once for all rotations -- 
      // en.calcPairlist();
      
      double ang_inc=360/rotsteps;
   // #ifdef OMP
   // #pragma omp parallel for firstprivate(fitsys, calcsys,en) reduction(+:nb,cov,tot)
   // #endif
      for (int a=0; a<rotsteps; a++) {
      
      // move center of the longest axes onto each other
      //fit::PositionUtils::translate(&fitsys, 
        // fit_atoms_sys.pos(maxdistatoms_sys[0])-fit_atoms_fit.pos(maxdistatoms_fit[0])-0.5*(fit_atoms_sys.pos(maxdistatoms_sys[1])-fit_atoms_sys.pos(maxdistatoms_sys[0])));
      // maybe better: move cogs of fitatoms and sysatoms onto each other
     // fit::PositionUtils::translate(&fitsys, 
      //              fit::PositionUtils::cog(sys,fit_atoms_sys)-fit::PositionUtils::cog(fitsys,fit_atoms_fit));

      // rotate around longest axis
      if (a) {
        gmath::Matrix mat_rot;
        if (randomrot) {
          double rnum=rand();
          double randomang = 360.0 / RAND_MAX * rnum;
          //cerr << "randomang " << randomang << " " << rnum << " " <<RAND_MAX<<endl;
          mat_rot=fit::PositionUtils::rotateAround(v2, randomang);
        } else {
          mat_rot=fit::PositionUtils::rotateAround(v2, ang_inc);
        }
        fit::PositionUtils::rotate(&fitsys, mat_rot);
      }
      fit::PositionUtils::translate(&fitsys, 
                    fit::PositionUtils::cog(sys,fit_atoms_sys)-fit::PositionUtils::cog(fitsys,fit_atoms_fit));
                    

      // calcsys gets the solute from fitsys
      for(int m=0; m< calcsys.numMolecules(); m++)
        for(int a=0; a < calcsys.mol(m).numAtoms(); a++)
          calcsys.mol(m).pos(a) = fitsys.mol(m).pos(a);
      
      
      // we have to gather with any method to get covalent interactions 
      // and charge-groups connected -- I already gathered solute of fitsys and sys, so calcsys should be ok? MP
      //pbc->gathergr();

      // calculate pairlist and energies
      en.calc();
      // calculate the energies assuming pairlist is already there
      // en.calcNb_interactions();
      // en.calcCov();

      // print any ouput you like
      cout.precision(10);
      cout.setf(ios::right, ios::adjustfield);
      cout << time
       << ' ' << setw(21) << en.cov()
           << ' ' << setw(21) << en.nb()
           << ' ' << setw(21) << en.tot()
       << endl;

      //store some averages
      cov+=en.cov()/rotsteps;
      nb+=en.nb()/rotsteps;
      tot+=en.tot()/rotsteps;
      
      if (write_coord) {
          oref->writeTimestep(time.steps(), time.time());
          *oref << sys; //prot;  
          ofit->writeTimestep(time.steps(), time.time());
          *ofit << calcsys; //prot;          
      }
      }
      
      num_frames++;
    }
  }
      if (write_coord) { 
  osref.close();
  osfit.close();
  }
  // print out averages
  if(num_frames * rotsteps>1){
    cout.precision(10);
    cout.setf(ios::right, ios::adjustfield);
    cout << endl << "# ave."
         << ' ' << setw(21) << cov/(num_frames)
         << ' ' << setw(21) << nb/(num_frames)
         << ' ' << setw(21) << tot/(num_frames)
         << endl;
  }
  t=clock()-t;
  printf ("# Runtime: %f seconds\n",((float)t)/CLOCKS_PER_SEC); 
}
 
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}


std::vector<int> max_dist_atoms(AtomSpecifier &as)
{
  // calculate the atoms furthest apart from each other
  // gathering has to be done elsewhere 
  
  std::vector<int> maxdistatoms (2,-99);
  if (as.size() < 2) 
    throw gromos::Exception("axis_rotate", 
                "AtomSpecifier has less than 2 atoms\n");
  
  double d=0;
  double dmax=0;
  int i_max=-1, j_max=-1;
  for(int i=0; i < as.size(); i++) {
    for(int j=i+1; j<as.size(); j++) {
      Vec v=as.pos(i)-as.pos(j);
      d=v.abs();
      if (d>dmax) {
        dmax=d;
        i_max=i;
        j_max=j;
      }
    }
  }
  maxdistatoms[0]=i_max;
  maxdistatoms[1]=j_max;
  return maxdistatoms;
}




