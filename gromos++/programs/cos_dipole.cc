/**
 * @file cos_dipole.cc
 * Calculate the induced and the fix dipole
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
  knowns << "topo" << "pbc" << "time" << "traj" << "trs"<< "molecules";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo <molecular topology file>\n";
  usage += "\t@pbc    <boundary type> [<gather method>]\n";
  usage += "\t[@time   <time and dt>]\n";
  usage += "\t@molecules <molecule numbers>(specify molecules for averaging)\n";
  usage += "\t@traj   <trajectory files>\n";
  usage += "\t@trs    <special traj with cosdiscplacement>\n";
  
 
  try{
    Arguments args(argc, argv, knowns, usage);

    // get the @time argument
    utils::Time time(args);
     
    //  read topology
    InTopology it(args["topo"]);
    System sys(it.system());

    System refSys(it.system());

    // create an atomspecifier that contains all atoms (including the solvent)
    //or just that ones that got specified by at molecules
    utils::AtomSpecifier atoms(sys);
    int start_mol=0;
    int stop_mol=0;
    int num_mol=0;
    
    for(int m=0; m< sys.numMolecules(); ++m)
	for(int a=0; a<sys.mol(m).numAtoms(); ++a)
	  atoms.addAtom(m,a);
    atoms.addSpecifier("s:a");

    if(args.count("molecules") == -1)
    {

      start_mol=0;
      stop_mol=sys.numMolecules();
      num_mol=sys.numMolecules();
      cout << " I am here" << endl; 
    }
    else
    {
      vector<int> molecules = args.getValues<int>("molecules", 2, false,
            Arguments::Default<int>() << 1 << 2);
      cout << "start molecules\t" << molecules[0] << endl;
      cout << "stop molecules\t" << molecules[1] << endl;
      //cout << "Upper bond:\t" << args.upper_bound("molecules") << endl;
      // cout << "lower bond:\t" << args.lower_bound("molecules") << endl;
      start_mol=molecules[0]-1;
      stop_mol=molecules[1];
      num_mol=molecules[1]-molecules[0]+1;
    }
    
    //determine net charge
    double ncharge=0;
    double nchargepa=0;
    for(int i=0; i < atoms.size(); i++){
      ncharge+=atoms.charge(i);
    }
    nchargepa=ncharge/atoms.size();
    
     
    
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
    cout << "#     time           total [D]       fix[D]           ind[D]\n";

    // prepare the calculation of the average volume
    double vol=0,sum_vol=0, vcorr=1.0;
    Arguments::const_iterator iter=args.lower_bound("pbc");
    if(iter!=args.upper_bound("pbc"))
      if(iter->second[0]=='t') vcorr=0.5;

    // and values to calculate the dipole fluctuation
    Vec sum_dip(0.0,0.0,0.0);
    double sum_dip2 = 0.0, fluc = 0.0;
    int numFrames=0;

    //Variable tos save thetime series of the dipole of the frame
    vector<double> tot_mol_dip;
    vector<double> tot_fix_dip;
    vector<double> tot_ind_dip;
    
    // define input coordinate
    InG96 ic,is;
    
     
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
	double sum_fix_dip=0;
	double sum_ind_dip=0;

	int pol_count=0;

	cout << "#amount of selected Molecules:\t" << num_mol << endl;

        for(int m=start_mol; m<stop_mol; m++) 
        {
	  Vec mol_dip(0,0,0);
	  Vec fix_dip(0,0,0);
	  Vec ind_dip(0,0,0);
            //cout << "#Amount of atoms in molecule:\t" << sys.mol(m).numAtoms() << endl;
            for(int a=0; a<sys.mol(m).numAtoms(); a++)
            {
	      
	        if(sys.mol(m).topology().atom(a).isPolarisable())
		{
		  pol_count++;
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
		    //sum of the ofsetpos and the cosdicplacment
		    Vec cos_pos = ofset_pos + sys.mol(m).cosDisplacement(a);
		    //calculation for the molecular dipol
		    mol_dip += ofset_pos*(sys.mol(m).topology().atom(a).charge()+(sys.mol(m).topology().atom(a).cosCharge()*-1));
		    mol_dip += cos_pos*sys.mol(m).topology().atom(a).cosCharge();
		    fix_dip += ofset_pos*(sys.mol(m).topology().atom(a).charge());
		    ind_dip += sys.mol(m).cosDisplacement(a)*sys.mol(m).topology().atom(a).cosCharge();
 
		  }
		  else
		  {
		    Vec cos_pos = sys.mol(m).pos(a) + sys.mol(m).cosDisplacement(a);

		    //calculation for the molecular dipol
		    mol_dip += sys.mol(m).pos(a)*(sys.mol(m).topology().atom(a).charge()+(sys.mol(m).topology().atom(a).cosCharge()*-1));
		    mol_dip += cos_pos*sys.mol(m).topology().atom(a).cosCharge();
		    fix_dip += sys.mol(m).pos(a)*(sys.mol(m).topology().atom(a).charge());
		    ind_dip += sys.mol(m).cosDisplacement(a)*sys.mol(m).topology().atom(a).cosCharge();
		  }
		}
		else
		{
		  //case no polarisation and offset
		  mol_dip += (sys.mol(m).topology().atom(a).charge()-nchargepa)*sys.mol(m).pos(a);
		  fix_dip += (sys.mol(m).topology().atom(a).charge()-nchargepa)*sys.mol(m).pos(a);
		}   
	    }
	    //cout << m << ":\t" << mol_dip.abs() << endl;
	    sum_mol_dip += mol_dip.abs();
	    sum_fix_dip += fix_dip.abs();
	    sum_ind_dip += ind_dip.abs();
        }
 
	tot_mol_dip.push_back(sum_mol_dip/num_mol);
	tot_fix_dip.push_back(sum_fix_dip/num_mol);
	if(pol_count==0)
	{
	  tot_ind_dip.push_back(0);
	}
	else
	{
	  tot_ind_dip.push_back(sum_ind_dip/pol_count);
	}
	// do some bookkeeping
	numFrames++;

	//cout << num_mol << "\t\t"<< sum_mol_dip << endl;
	cout << time
	     << setw(15) << setprecision(8) << sum_mol_dip*48.032045/num_mol
	     << setw(15) << setprecision(8) << sum_fix_dip*48.032045/num_mol
	     << setw(15) << setprecision(8) << sum_ind_dip*48.032045/pol_count
	     << endl;
      }
      ic.close();
    }
    //now print avergae over the systems for all three Values
    double sum1=0,sum2=0,sum3=0;
    for(int i=0; i < tot_mol_dip.size(); i++)
    {
      sum1+=tot_mol_dip[i];
      sum2+=tot_fix_dip[i];
      sum3+=tot_ind_dip[i];
    }
    cout << "#Averages:\n" << setprecision(8) << sum1*48.032045/numFrames
	 << setw(15) << setprecision(8) << sum2*48.032045/numFrames
	 << setw(15) << setprecision(8) << sum3*48.032045/numFrames
	 << endl;

  }
  
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}

