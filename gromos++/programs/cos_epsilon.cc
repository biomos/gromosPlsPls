/**
 * @file epsilon.cc
 * Calculate the relative permittivity over a trajectory
 */

/**
 * @page programs Program Documentation
 *
 * @anchor epsilon
 * @section epsilon Calculate the relative permittivity over a trajectory
 * @author @ref co
 * @date 13-6-07
 *
 * Program epsilon estimates the relative dielectric permittivity, 
 * @f$\epsilon(0)@f$, of simulation box from a Kirkwood - Fr&ouml;hlich type of
 * equation, as derived by Neumann [Mol. Phys. 50, 841 (1983)]
 *
 * @f[ (\epsilon(0) - 1) \frac{2\epsilon_{RF}+1}{2\epsilon_{RF}+\epsilon(0)} = \frac{<\vec{M}^2> - <\vec{M}>^2}{3\epsilon_0 Vk_BT} @f]
 *
 * where @f$\vec{M}@f$ is the total dipole moment of the system, 
 * @f$\epsilon_0@f$ is the dielectric permittivity of vacuum,
 * @f$\epsilon_{RF}@f$ is a reaction-field epsilon value, @f$V@f$ is the volume
 * and @f$k_BT@f$ is the absolute temperature multiplied by the Boltzmann 
 * constant. 
 * 
 * Note that the total dipole moment of the system is only well-defined for 
 * systems consisting of neutral molecules. If the system carries a net-charge,
 * the dipole moment will depend on the position of the origin. In cases where
 * the overall system is neutral but contains ions, the dipole moment will
 * depend on which periodic copy of the ions is taken. In these cases, epsilon
 * issues a warning that results will critically depend on the choice of
 * gathering method.
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;boundary type&gt; [&lt;gather method&gt;] </td></tr>
 * <tr><td> [\@time</td><td>&lt;@ref utils::Time "time and dt"&gt;] </td></tr>
 * <tr><td> [\@e_rf</td><td>&lt;reaction field epsilon&gt;] </td></tr>
 * <tr><td> \@temp</td><td>&lt;temperature&gt; </td></tr>
 * <tr><td> \@traj</td><td>&lt;trajectory files&gt; </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  epsilon
    @topo  ex.top
    @pbc   r
    [@time  0 0.2]
    @e_rf  61
    @temp  300
    @traj  ex.tr
 @endverbatim
 *
 * <hr>
 */

#include <vector>
#include <iomanip>
#include <iostream>
#include <cassert>
#include <locale>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/fit/Reference.h"
#include "../src/gio/InG96.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/Box.h"
#include "../src/gio/InTopology.h"
#include "../src/bound/Boundary.h"
#include "../src/fit/PositionUtils.h"
#include "../src/gmath/Vec.h"
#include "../src/gmath/Physics.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/utils/groTime.h"

using namespace std;
using namespace fit;
using namespace gcore;
using namespace gio;
using namespace bound;
using namespace args;

int main(int argc, char **argv){

  Argument_List knowns;
  knowns << "topo" << "pbc" << "time" << "e_rf" << "temp" << "traj" << "trs";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo <molecular topology file>\n";
  usage += "\t@pbc    <boundary type> [<gather method>]\n";
  usage += "\t[@time   <time and dt>]\n";
  usage += "\t[@e_rf   <reaction field epsilon>]\n";
  usage += "\t@temp   <temperature>\n";
  usage += "\t@traj   <trajectory files>\n";
  usage += "\t@trs    <special traj with cosdiscplacement>\n";
  
 
  try{
    Arguments args(argc, argv, knowns, usage);

    // get the @time argument
    utils::Time time(args);
  
    // read the temperature
    double temp = args.getValue<double>("temp");
    
    //  read topology
    InTopology it(args["topo"]);
    System sys(it.system());

    System refSys(it.system());

    // create an atomspecifier that contains all atoms (including the solvent)
    utils::AtomSpecifier atoms(sys);
    for(int m=0; m< sys.numMolecules(); ++m)
      for(int a=0; a<sys.mol(m).numAtoms(); ++a)
	atoms.addAtom(m,a);
    atoms.addSpecifier("s:a");
    
    //determine net charge
    double ncharge=0;
    double nchargepa=0;
    for(int i=0; i < atoms.size(); i++){
      ncharge+=atoms.charge(i);
    }
    nchargepa=ncharge/atoms.size();
    
    // read e_rf
    double e_rf= args.getValue<double>("e_rf", false, 1.0);
    
    
    // parse boundary conditions
    Boundary *pbc = BoundaryParser::boundary(sys, args);
    // parse gather method
    Boundary::MemPtr gathmethod = args::GatherParser::parse(sys,refSys,args);
    
    // write a title
    cout << "#\n";
    if(ncharge!=0.0){
      cout << "# WARNING the system carries a net charge ( "
	   << ncharge << ")\n"
	   << "#         this means that the dipole depends on the origin\n"
	   << "#         and your estimate of eps is probably wrong\n"
	   << "#\n";
    } 
    else {
      // we also want to know if there are any ions, we'll restrict that
      // possibility to the solute
      int ion_count=0;
      for(int m=0; m<sys.numMolecules(); ++m) {
	double tc=0.0;
	for(int a=0; a< sys.mol(m).numAtoms(); ++a)
	  tc+=sys.mol(m).topology().atom(a).charge();
	if(tc!=0.0) ion_count++;
      }
      if(ion_count!=0)
	cout << "# WARNING although the system is neutral overall, there are\n"
	     << "#         "  << ion_count << " molecules that carry a charge\n"
	     << "#         the system-dipole will depend on their positions in "
	     << "the periodic box\n"
	     << "#         this is likely to give random results\n";
    }
    cout << "#\n";
    cout << "#     time    dipole       epsilon\n";

    // prepare the calculation of the average volume
    double vol=0,sum_vol=0, vcorr=1.0;
    Arguments::const_iterator iter=args.lower_bound("pbc");
    if(iter!=args.upper_bound("pbc"))
      if(iter->second[0]=='t') vcorr=0.5;

    // and values to calculate the dipole fluctuation
    Vec sum_dip(0.0,0.0,0.0);
    double sum_dip2 = 0.0, fluc = 0.0;
    int numFrames=0;
    
    // define input coordinate
    InG96 ic,is;
    
    // we also need these factors    
    double fac, a, b, eps;
    double f=3.0 * gmath::physConst.get_eps0() * gmath::physConst.get_boltzmann() * temp*(2*e_rf + 1.0);
    
    // loop over all trajectories
    for(Arguments::const_iterator 
	  iter=args.lower_bound("traj"),
	  to=args.upper_bound("traj"),
          siter=args.lower_bound("trs"),
          sto=args.upper_bound("trs");
	iter!=to || siter!=sto; ++iter,++siter){
      
      // open file
      ic.open((iter->second).c_str());
      ic.select("ALL");
      
      //open special traj
      is.open((siter->second).c_str());
      
      // loop over single trajectory
      while(!ic.eof() && !is.eof()){
        is >> sys >> time;
        
        //save trime for checkign later
        double time_trs=time.time();
        
	ic >> sys >> time;
        
        //hier check if the time are the same
        if(time_trs != time.time())
        {
            cerr << "Time is not corresponding between trajectories\n" <<
                    "Specieal traj:\t" << time_trs << "\n" <<
                    "Position traj:\t" << time.time() << endl;
            return -1;
        }
        
	// we have to reconnect the molecules
	// in the case of neutral molecules, this will be sufficient
	// non-neutral molecules will behave weirdly.
	(*pbc.*gathmethod)();

	// calculate the volume
	sys.box().update_triclinic();
        vol=vcorr*sys.box().K_L_M();
	sum_vol+=vol;
	
	Vec dipole(0.0,0.0,0.0);
        double sum_mol_dip=0;
 
        for(int m=0; m<sys.numMolecules(); m++) 
        {
	  Vec mol_dip(0,0,0);
            //cout << "Amount of atoms in molecule:\t" << sys.mol(m).numAtoms() << endl;
            for(int a=0; a<sys.mol(m).numAtoms(); a++)
            {
	      
	        if(sys.mol(m).topology().atom(a).isPolarisable())
		{
		  if(sys.mol(m).topology().atom(a).poloffsiteGamma()>0.000001)
                  {
		    //cout << "Number: atom i :\t" << sys.mol(m).topology().atom(a).poloffsiteI() <<endl;
                    //cout << "Number: atom j :\t" << sys.mol(m).topology().atom(a).poloffsiteJ() <<endl;
                    //cout << "Position i:\t" << sys.mol(m).pos(sys.mol(m).topology().atom(a).poloffsiteI()+a)[0] << "\t"
		    // << sys.mol(m).pos(sys.mol(m).topology().atom(a).poloffsiteI()+a)[1] << "\t"
		    //	 << sys.mol(m).pos(sys.mol(m).topology().atom(a).poloffsiteI()+a)[2] << endl;
                    //cout << "Position j:\t" << sys.mol(m).pos(sys.mol(m).topology().atom(a).poloffsiteJ()+a)[0] << "\t"
		    //	 << sys.mol(m).pos(sys.mol(m).topology().atom(a).poloffsiteJ()+a)[1] << "\t"
		    //     << sys.mol(m).pos(sys.mol(m).topology().atom(a).poloffsiteJ()+a)[2] << endl;
		    //cout <<  sys.mol(m).pos(a+1)[0] << "\t" <<  sys.mol(m).pos(a+1)[1] <<  sys.mol(m).pos(a+1)[2] << endl;
		    //cout <<  sys.mol(m).pos(a+2)[0] << "\t" <<  sys.mol(m).pos(a+2)[1] <<  sys.mol(m).pos(a+2)[2] <<endl;
                            
                    //position of the offset atom
                    Vec ofset_pos = sys.mol(m).pos(a) + sys.mol(m).topology().atom(a).poloffsiteGamma()*
                                    (sys.mol(m).pos(sys.mol(m).topology().atom(a).poloffsiteI()+a)
                                    +sys.mol(m).pos(sys.mol(m).topology().atom(a).poloffsiteJ()+a)
                                    -2*sys.mol(m).pos(a));
		    //dipole contribution of the offset atom
                    dipole+=ofset_pos*(sys.mol(m).topology().atom(a).charge()+(sys.mol(m).topology().atom(a).cosCharge()*-1));

		    //sum of the ofsetpos and the cosdicplacment
                    Vec cos_pos = ofset_pos + sys.mol(m).cosDisplacement(a);
                    //dipole of the virrtual cos particle
                    dipole += cos_pos*sys.mol(m).topology().atom(a).cosCharge();

		    //calculation for the molecular dipol
		    mol_dip += ofset_pos*(sys.mol(m).topology().atom(a).charge()+(sys.mol(m).topology().atom(a).cosCharge()*-1));
		    mol_dip += cos_pos*sys.mol(m).topology().atom(a).cosCharge();
 
		  }
		  else
		  {
		    //contribution to the box dipole no offset but polarisable
		    dipole+=sys.mol(m).pos(a)*(sys.mol(m).topology().atom(a).charge()+(sys.mol(m).topology().atom(a).cosCharge()*-1));
		    Vec cos_pos = sys.mol(m).pos(a) + sys.mol(m).cosDisplacement(a);
		    //contribution of the virtual cos particle
		    dipole += cos_pos*sys.mol(m).topology().atom(a).cosCharge();

		    //calculation for the molecular dipol
		      mol_dip += sys.mol(m).pos(a)*(sys.mol(m).topology().atom(a).charge()+(sys.mol(m).topology().atom(a).cosCharge()*-1));
		      mol_dip += cos_pos*sys.mol(m).topology().atom(a).cosCharge();    
		  }
		}
		else
		{
		  //case no polarisation and offset
		  dipole += (sys.mol(m).topology().atom(a).charge()-nchargepa)*sys.mol(m).pos(a);
		  mol_dip += (sys.mol(m).topology().atom(a).charge()-nchargepa)*sys.mol(m).pos(a);
		}   
	    }
	    sum_mol_dip += mol_dip.abs(); 
        }
        sum_mol_dip = sum_mol_dip/sys.numMolecules();
	

	// do some bookkeeping
	numFrames++;

        sum_dip2+= dipole.abs2();
	sum_dip += dipole;

	// calculate the fluctuations of the dipole moment
        fluc=sum_dip2 / numFrames - (sum_dip / numFrames).abs2();

	// calculate the current estimate of eps
	fac = f * sum_vol/numFrames;
        a = 2 * e_rf * fluc + fac;
	b= -fluc + fac;
	eps = a/b;
	
	cout << time
	     << setw(15) << setprecision(8) << dipole.abs()
	     << setw(15) << setprecision(8) << sum_mol_dip*48.032045
	     << setw(15) << setprecision(8) << eps
	     << endl;
      }
      ic.close();
    }
  }
  
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}

