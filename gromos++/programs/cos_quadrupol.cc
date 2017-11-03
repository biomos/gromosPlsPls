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
  knowns << "topo" << "pbc" << "time" << "offset" << "atoms" << "traj" << "trs";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo <molecular topology file>\n";
  usage += "\t@pbc    <boundary type> [<gather method>]\n";
  usage += "\t[@time   <time and dt>]\n";
  usage += "\t[@offset (offset of the cose site in gamma/2)]\n";
  usage += "\t[@atoms   <atoms to include> (default all solute)]\n";
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

    utils::AtomSpecifier atoms(sys);
    
    if(args.count("atoms")>0)
    {
        for(Arguments::const_iterator iter=args.lower_bound("atoms"), 
		    to=args.upper_bound("atoms"); iter!=to; ++iter)
        {
	  atoms.addSpecifier(iter->second);
	}   
    } 
    else 
    {
	for(int m=0; m<sys.numMolecules(); m++)
          for(int a=0; a<sys.mol(0).numAtoms(); a++)
		    atoms.addAtom(m,a);

    }
    
    //get the COS sites to be included
    utils::AtomSpecifier cosatoms(sys);
    bool ltrs=false;
    if(args.count("trs") >= 0) ltrs=true;
    if(ltrs) 
    {
        //	if(args.count("cosatoms")>0) {
	for(Arguments::const_iterator iter=args.lower_bound("atoms"), 
            to=args.upper_bound("atoms"); iter!=to; ++iter)
        {
                atoms.addSpecifier(iter->second);
	}
    } 
    else 
    {
        for(int m=0; m<sys.numMolecules(); m++)
	    for(int a=0; a<sys.mol(0).numAtoms(); a++)
		    if(sys.mol(m).topology().atom(a).cosCharge())
			cosatoms.addAtom(m,a);
    }
       
    //read offset when neccesary othewise it is 0
    double off_set_gamma=0;
    if(args.count("offset")==1)
    {
        off_set_gamma=args.getValue<double>("offset", true);
    }
    
    // parse boundary conditions
    Boundary *pbc = BoundaryParser::boundary(sys, args);
    // parse gather method
    Boundary::MemPtr gathmethod = args::GatherParser::parse(sys,refSys,args);
    
       
    //header for the time series of the dipole and qudrupol moments
    cout << "# Time\t\t  dipole [D]\t\t fix dipole [D]\t\tQ_xx [D nm]\t\tQ_yy [D nm]\t\tQ_zz[D nm]" << endl;

    // and values to calculate the dipole fluctuation
    double tot_dip=0;
    double fix_dip=0;
    
    //the value of the sum of qudrupol moments
    double tot_qxx=0;
    double tot_qyy=0;
    double tot_qzz=0;
    
    //counter for frames
    int numFrames=0;
    
    // define input coordinate
    InG96 ic,is;
    
   
    
    // loop over all trajectories
    for(Arguments::const_iterator 
	  iter=args.lower_bound("traj"),
	  to=args.upper_bound("traj"),
          siter=args.lower_bound("trs"),
          sto=args.upper_bound("trs");
	iter!=to || siter!=sto; ++iter,++siter)
    {
      // open file
      ic.open((iter->second).c_str());
      ic.select("ALL");
      //open special traj
      is.open((siter->second).c_str());
      
      // loop over single trajectory
      while(!ic.eof() && !is.eof())
      {
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

	
        //variable to hold the add up of the single mol. dipol moments of the frame
        double tot_dip_frame=0;
	double fix_dip_frame=0;
        
        //variable to add up the three quadrupol moments
        double tot_qxx_frame=0;
        double tot_qyy_frame=0;
        double tot_qzz_frame=0;
        
        //variable to count how many molecules
        int count=0;
        
        //cout << "Amount of molecules:\t" << sys.numMolecules() << endl;
        for(int m=0; m<sys.numMolecules(); m++) 
        {
            //variable to keep the dipol for the molecule sum over all atom
            Vec dipole(0,0,0);
	    Vec fix_dipole(0,0,0);
            
            //variable to save the molecular components of the three quadrupols
            double qxx=0;
            double qyy=0;
            double qzz=0;
            
            //cout << "Amount of atoms in molecule:\t" << sys.mol(m).numAtoms() << endl;
            for(int a=0; a<sys.mol(m).numAtoms(); a++)
            {
                if(off_set_gamma==0)
                {
		  
                    dipole += (sys.mol(m).topology().atom(a).charge())*sys.mol(m).pos(a);
		    fix_dipole += (sys.mol(m).topology().atom(a).charge())*sys.mol(m).pos(a);
                    //cout << "atom charge\t" << sys.mol(m).topology().atom(a).charge() << endl;
                    //cout << "distance suqred:\t" << sys.mol(m).pos(a).abs2() << endl;
                    //cout << "3*pos_x*pos_x:\t" << 3*sys.mol(m).pos(a)[0]*sys.mol(m).pos(a)[0] << endl;
                    qxx += (3*sys.mol(m).pos(a)[0]*sys.mol(m).pos(a)[0]-sys.mol(m).pos(a).abs2())*sys.mol(m).topology().atom(a).charge();
                    //cout << "qxx:\t" << qxx << endl;
                    qyy += (3*sys.mol(m).pos(a)[1]*sys.mol(m).pos(a)[1]-sys.mol(m).pos(a).abs2())*sys.mol(m).topology().atom(a).charge();
                    qzz += (3*sys.mol(m).pos(a)[2]*sys.mol(m).pos(a)[2]-sys.mol(m).pos(a).abs2())*sys.mol(m).topology().atom(a).charge();
                }
                else
                {
		  
                    if(sys.mol(m).cosDisplacement(a).abs2()!=0)
                    {
		      
                        Vec ofset_pos = sys.mol(m).pos(a) + off_set_gamma*
                                        (sys.mol(m).pos(a+1)+sys.mol(m).pos(a+2)-2*sys.mol(m).pos(a));
			
                        dipole+=ofset_pos*(sys.mol(m).topology().atom(a).charge()+(sys.mol(m).topology().atom(a).cosCharge()*-1));
			fix_dipole+=ofset_pos*(sys.mol(m).topology().atom(a).charge());
                        Vec cos_pos = ofset_pos + sys.mol(m).cosDisplacement(a);
                        dipole += cos_pos*sys.mol(m).topology().atom(a).cosCharge();
                        
                        //now adding the two contibution of the offset charge and the cos charge
                        qxx += (3*ofset_pos[0]*ofset_pos[0]-ofset_pos.abs2())*(sys.mol(m).topology().atom(a).charge()+(sys.mol(m).topology().atom(a).cosCharge()*-1));
                        qyy += (3*ofset_pos[1]*ofset_pos[1]-ofset_pos.abs2())*(sys.mol(m).topology().atom(a).charge()+(sys.mol(m).topology().atom(a).cosCharge()*-1));
                        qxx += (3*ofset_pos[2]*ofset_pos[2]-ofset_pos.abs2())*(sys.mol(m).topology().atom(a).charge()+(sys.mol(m).topology().atom(a).cosCharge()*-1));
                        //now the cos charges with cos position
                        qxx += (3*cos_pos[0]*cos_pos[0]-cos_pos.abs2())*sys.mol(m).topology().atom(a).cosCharge();
                        qyy += (3*cos_pos[1]*cos_pos[1]-cos_pos.abs2())*sys.mol(m).topology().atom(a).cosCharge();
                        qzz += (3*cos_pos[2]*cos_pos[2]-cos_pos.abs2())*sys.mol(m).topology().atom(a).cosCharge();
                    }
                    else
                    {
                        dipole += (sys.mol(m).topology().atom(a).charge())*sys.mol(m).pos(a);
			fix_dipole += (sys.mol(m).topology().atom(a).charge())*sys.mol(m).pos(a);
                        qxx += (3*sys.mol(m).pos(a)[0]*sys.mol(m).pos(a)[0]-sys.mol(m).pos(a).abs2())*sys.mol(m).topology().atom(a).charge();
                        qyy += (3*sys.mol(m).pos(a)[1]*sys.mol(m).pos(a)[1]-sys.mol(m).pos(a).abs2())*sys.mol(m).topology().atom(a).charge();
                        qzz += (3*sys.mol(m).pos(a)[2]*sys.mol(m).pos(a)[2]-sys.mol(m).pos(a).abs2())*sys.mol(m).topology().atom(a).charge();
                    }
                }
             }
             //molecular dipol is calculate add to total of frame
             tot_dip_frame+=dipole.abs();
	     fix_dip_frame+=fix_dipole.abs();
             
             //add up th emolecular quadrupol moments
             tot_qxx_frame+=qxx;
             tot_qyy_frame+=qyy;
             tot_qzz_frame+=qzz;
             
             //cout << numFrames << "\t" << tot_dip_frame << endl; 
             count++;
        }
	//add the average dipol of the frame
        tot_dip+=tot_dip_frame/count;
	fix_dip+=fix_dip_frame/count;
        tot_qxx+=tot_qxx_frame/count;
        tot_qyy+=tot_qyy_frame/count;
        tot_qzz+=tot_qzz_frame/count;

	// do some bookkeeping
	numFrames++;

	cout << time
	     << setw(15) << setprecision(8) << (tot_dip_frame/count)*48.032045
             << setw(15) << setprecision(8) << (fix_dip_frame/count)*48.032045
             //<< setw(15) << setprecision(8) << (tot_qxx_frame/count)*48.032045
             //<< setw(15) << setprecision(8) << (tot_qyy_frame/count)*48.032045
             //<< setw(15) << setprecision(8) << (tot_qzz_frame/count)*48.032045
              << setw(15) << setprecision(8) << 0
             << setw(15) << setprecision(8) << 0
             << setw(15) << setprecision(8) << 0 
	     << endl;
      }
      ic.close();
      is.close();
    }
    //cout of the averge results
    cout <<"# Average: " 
         << setw(15) << setprecision(8) << (tot_dip/numFrames)* 48.032045
	 << setw(15) << setprecision(8) << (fix_dip/numFrames)* 48.032045
         //<< setw(15) << setprecision(8) << (tot_qxx/numFrames)*48.032045
         //<< setw(15) << setprecision(8) << (tot_qyy/numFrames)*48.032045
         //<< setw(15) << setprecision(8) << (tot_qzz/numFrames)*48.032045
            << setw(15) << setprecision(8) << 0
         << setw(15) << setprecision(8) << 0
         << setw(15) << setprecision(8) << 0
         << endl;
  }
  
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}

