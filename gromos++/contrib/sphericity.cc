/**
 * @file sphericity.cc
 * Calculates the sphericity of a specified set of atoms using the moment of inertia as criterion.
 */

/**
 * @page programs Program Documentation
 *
 * @anchor sphericity
 * @section sphericity Calculates the sphericity of a specified set of atoms
 * @author @ref ms
 * @date 20-01-2018
 *
 * Calculates the sphericity of a specified set of atoms (<b>\@atoms</b>) using the moment of inertia (I) as criterion.
 *
 * Sphericity is calculated from Salaniwal et al. "Molecular Simulation of a Dichain Surfactant/Water/Carbon Dioxide System. 1. Structural Properties of Aggregates", Langmuir, 2001:
 *
 * @f[ S=1-\frac{min(I_x,I_y,I_z)}{average(I_x,I_y,I_z)} @f]
 *
 * Thereby, 0 represents a perfect sphere and 1 a highly non-spherical shape.
 *
 * Since the moments of inertia (Ix, Iy, Iz) depend heavily on the alignment of the body along the Carthesian axes, the user can specify 
 * which atoms should be used to calculate the necessary rotation to align <b>\@atoms</b> along the z- and y-axis with <b>\@alignatoms</b>
 * (largest extend in 1. z-direction and 2. y-direction). By default <b>\@alignatoms</b> = <b>\@atoms</b>.
 *
 * Sphericity is OpenMP-parallelised by trajectory files (each thread works on one <b>\@traj</b> file). Specify number of threads by <b>\@cpus</b> flag.
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;boundary type&gt; </td></tr>
 * <tr><td> \@atoms</td><td>&lt;@ref AtomSpecifier "atoms" to calculate sphericity for&gt; </td></tr>
 * <tr><td> [\@alignatoms</td><td>&lt;@ref AtomSpecifier "atoms" used to align \@atoms along Carthesian z- and y-axes (by rotation). Default: \@atoms&gt;]</td></tr>
 * <tr><td> [\@time</td><td>&lt;@ref utils::Time "time and dt"&gt;]</td></tr>
 * <tr><td> [\@cpus</td><td>&lt;number of CPUs to use (default: 1)&gt;]</td></tr>
 * <tr><td> \@traj</td><td>&lt;trajectory files&gt; </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  rot_rel
    @topo  ex.top
    @pbc   r cog
    @time  0 1
    @atoms  1:a
    @alignatoms 1:res(3:a)
    @cpus 5
    @traj ex1.trc ex2.trc ex3.trc ex4.trc ex5.trc
 @endverbatim
 *
 * <hr>
 */


#include <cassert>
#include <vector>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <fstream>

#ifdef OMP
#include <omp.h>
#endif

#include "../src/utils/AtomSpecifier.h"
#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/gio/InG96.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/LJException.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gcore/Box.h"
#include "../src/gio/InTopology.h"
#include "../src/bound/Boundary.h"
#include "../src/fit/PositionUtils.h"
#include "../src/utils/Value.h"
#include "../src/utils/VectorSpecifier.h"
#include "../src/gmath/Vec.h"
#include "../src/utils/groTime.h"
#include "../src/gmath/Matrix.h"
#include "../src/args/OutformatParser.h"
#include "../src/gio/OutCoordinates.h"

using namespace fit;
using namespace gcore;
using namespace gio;
using namespace bound;
using namespace args;
using namespace std;
using namespace utils;
using namespace gmath;

void rotate_solute(System &sys, AtomSpecifier &atoms, AtomSpecifier &rotation_atoms);
Vec calc_max_vector(AtomSpecifier &as_solute, int dim);  // modified from sim_box.cc

int main(int argc, char **argv) {

    Argument_List knowns;
    knowns << "topo" << "pbc" << "time" << "traj" << "atoms" << "alignatoms" << "cpus"; // << "outformat";

    string usage = argv[0];
    usage += "\n\t@topo    <molecular topology file>\n";
    usage += "\t@pbc     <boundary type [gather method]>\n";
    usage += "\t@atoms     <atoms to include in sphericity calculation>\n";
    usage += "\t[@alignatoms     <atoms used to align @atoms along Carthesian z- and y-axis. Default: @atoms]>\n";
    usage += "\t[@time    <time and dt>]\n";
    usage += "\t[@cpus    <numbers of CPUs to use. Default: 1>]\n";
    usage += "\t@traj    <trajectory files>\n";

    try {
        Arguments args(argc, argv, knowns, usage);

        // read topology
        InTopology it(args["topo"]);
        System sys(it.system());


        // get simulation time
        Time time(args);
        double time_dt = time.dt();
        double time_start = time.start_time();

        // parse boundary conditions



        //@traj
        const int traj_size = args.count("traj"); //number of trajectory files
        Arguments::const_iterator it_arg = args.lower_bound("traj");
        vector<Arguments::const_iterator> traj_array(traj_size);//array with pointers to trajectories:only way to go with omp

        for(int i = 0; i < traj_size; ++it_arg, ++i) {
            traj_array[i] = it_arg;
        }
        typedef map<int, vector< vector<double> > > TrajMap;
        TrajMap traj_map;

        //@cpus
        int num_cpus = 1;

        it_arg = args.lower_bound("cpus");
        if(it_arg != args.upper_bound("cpus")) {
            std::istringstream is(it_arg->second);
            is >> num_cpus;
            if(num_cpus <= 0)
                throw gromos::Exception("sphericity", "You must specify a number >0 for @cpus\n\n" + usage);
#ifdef OMP
            if(num_cpus > traj_size) {
                if(traj_size > omp_get_max_threads())
                    num_cpus = omp_get_max_threads();
                else
                    num_cpus = traj_size;
                cerr << "# Number of threads > number of trajectory files: not feasible. Corrected to " << num_cpus << " threads." << endl;
            }

            if(num_cpus > omp_get_max_threads()) {
                cerr << "# You specified " << num_cpus << " number of threads. There are only " << omp_get_max_threads() << " threads available." << endl;
                num_cpus = omp_get_max_threads();
            }
#else
            if(num_cpus != 1)
                throw gromos::Exception("sphericity", "@cpus: Your compilation does not support multiple threads. Use --enable-openmp for compilation.\n\n" + usage);
#endif

        }
#ifdef OMP
        omp_set_num_threads(num_cpus); //set the number of cpus for the parallel section
#endif // OMP
        cout << "# Number of threads: " << num_cpus << endl;


        // define input coordinate
//        string ext;
//        ofstream os;
//        OutCoordinates *oc = OutformatParser::parse(args, ext);
//        string inc = "SOLUTE";
//        ostringstream pdbName;
//        pdbName << "sphericity_coords" << ext;
//        string file=pdbName.str();
//        os.open(file.c_str());
//        oc->open(os);
//        oc->select(inc);
//        oc->writeTitle(file);

#ifdef OMP
        #pragma omp parallel for schedule (dynamic,1) firstprivate(sys, time)
#endif
        for(int traj = 0 ; traj < traj_size; ++traj) {
            double frame_time = time_start - time_dt;
            vector<double> time_vec;
            vector<double> spher_vec;

#ifdef OMP
            #pragma omp critical
#endif
            cout << "# Processing file: " << traj_array[traj]->second << endl;

            AtomSpecifier atoms(sys);
            for(Arguments::const_iterator
                    iter = args.lower_bound("atoms");
                    iter != args.upper_bound("atoms");
                    ++iter)
                atoms.addSpecifier(iter->second);

            AtomSpecifier alignatoms(sys);
            if(args.count("alignatoms") <= 0)
                alignatoms = atoms;
            else
                for(Arguments::const_iterator
                        iter = args.lower_bound("alignatoms");
                        iter != args.upper_bound("alignatoms");
                        ++iter)
                    alignatoms.addSpecifier(iter->second);
#ifdef OMP
            #pragma omp critical
#endif
            {
                if(atoms.empty())
                    throw gromos::Exception("sphericity", "No atoms in @atoms");
                if(alignatoms.empty())
                    throw gromos::Exception("sphericity", "No atoms in @alignatoms");
            }
            // open file
            InG96 ic;
            ic.open(traj_array[traj]->second);
            Boundary::MemPtr gathmethod; // must define variable first, otherwise compiler complains due to critical section.
            // must be critical, otherwise stdout is mangled with the couts of the gathermethod:
#ifdef OMP
            #pragma omp critical
#endif
            gathmethod = args::GatherParser::parse(sys, sys, args);
            Boundary *pbc = BoundaryParser::boundary(sys, args);
            vector<double> all_I;
            all_I.resize(3);

            // loop over single trajectory
            while(!ic.eof()) {
//                ic.select(inc);
                if(time.read()) {
                    ic >> sys >> time;//read coordinates & time from file
                    frame_time = time.time(); //get current time
                } else {
                    ic >> sys;
                    frame_time += time_dt; //numbering starts at 0 for every traj, correct overall times are generated in printstatistics
                }
                (*pbc.*gathmethod)();

                // rotate system
                rotate_solute(sys, atoms, alignatoms);

                double xx, yy, zz, mass, Imin, spher;
                double Ix = 0, Iy = 0, Iz = 0;

                for(int i = 0; i < atoms.size(); ++i) {
                    Vec pos = atoms.pos(i);
                    xx = pos[0] * pos[0];
                    yy = pos[1] * pos[1];
                    zz = pos[2] * pos[2];
                    mass = atoms.mass(i);
                    Iz += mass * (xx + yy);
                    Iy += mass * (xx + zz);
                    Ix += mass * (zz + yy);
                }
                all_I[0] = Ix;
                all_I[1] = Iy;
                all_I[2] = Iz;

                Imin = all_I[0];
                for(int i = 1; i < all_I.size(); ++i)
                    Imin = min(Imin, all_I[i]);

                spher = 1 - Imin / ((Ix + Iy + Iz) / 3.0);

                spher_vec.push_back(spher);
                time_vec.push_back(frame_time);

                // write new coordinates

//                oc->writeTimestep(time.steps(), time.time());
                //                *oc << sys;
            }
            ic.close();
            vector< vector<double> > traj_output;
            traj_output.push_back(time_vec);
            traj_output.push_back(spher_vec);
#ifdef OMP
            #pragma omp critical
#endif
            traj_map[traj] = traj_output;
        } // loop over traj end
        //        os.close();
        cout << "# perfect round shape: sphericity = 0. maximum deviation: sphericity=1" << endl;
        cout << setw(10) << "# time" << " " << setw(10) << "sphericity" << endl;

        map<int, double > start_tme;
        for(int trj = 0; trj < traj_map.size(); ++trj) {
            if(!traj_map.count(trj))
                continue;
            start_tme[trj] = traj_map[trj][0].size() * time_dt;
        }
        double start_time = time_start;
        double tme;

        for(int trj = 0; trj < traj_map.size(); ++trj) {
            if(!traj_map.count(trj))
                continue;
            for(int n = 0; n < traj_map[trj][0].size(); ++n) {
                tme = traj_map[trj][0][n];
                if(!time.read())
                    tme += start_time;
                cout << setw(10) << tme
                     << " "
                     << setw(10) << traj_map[trj][1][n]
                     << endl;
            }
            start_time += start_tme[trj];
        }
        cout << "# sphericity finished" << endl;
    }

    catch(const gromos::Exception &e) {
        cerr << e.what() << endl;
        exit(1);
    }
    return 0;
}


void rotate_solute(System &sys, AtomSpecifier &atoms, AtomSpecifier &rotation_atoms) {
    // this will be a two stage function.
    // First calculate the maximum distance between any two solute atoms
    // rotate the solute such that the atoms in as are along the z-axis
    Vec v = calc_max_vector(rotation_atoms, 3);
    AtomSpecifier all_atoms = atoms + rotation_atoms;

    double r = v.abs();
    double r_yz = sqrt(v[1] * v[1] + v[2] * v[2]);

    // the rotation matrix is the product of two rotations
    // 1. around the x-axis by theta
    //    with sin(theta) = v[1]/r_yz; cos(theta) = v[2]/r_yz
    // 2. around the y-axis by phi
    //    with sin(phi) = v[0]/r; cos(phi) = r_yz / r
    if(r == 0.0 || r_yz == 0.0) {
        throw gromos::Exception("sphericity",
                                "rotation failed (z-dimension).");
    }

    Matrix rot1(Vec(r_yz / r         ,  0         , v[0] / r),
                Vec(-v[0]*v[1] / r / r_yz ,  v[2] / r_yz , v[1] / r),
                Vec(-v[0]*v[2] / r / r_yz , -v[1] / r_yz , v[2] / r));

    for(int a = 0; a < all_atoms.size(); a++)
        sys.mol(all_atoms.mol(a)).pos(all_atoms.atom(a)) = rot1 * all_atoms.pos(a);
    
    // calculate the maximum distance in the x-y-plane
    v = calc_max_vector(rotation_atoms, 2);

    // rotate the solute around the z-axis, such that the atoms in as are
    // along the y-axis, this is done by a rotation around psi with
    // sin(psi) = x/r_xy; cos(psi) = y/r_xy;
    double r_xy = sqrt(v[0] * v[0] + v[1] * v[1]);

    if(r_xy == 0.0) {
        throw gromos::Exception("sphericity",
                                "rotation failed (y-dimension).");
    }

    Matrix rot2(Vec(+v[1] / r_xy ,  v[0] / r_xy , 0),
                Vec(-v[0] / r_xy ,  v[1] / r_xy , 0),
                Vec(0         ,  0         , 1));

    for(int a = 0; a < all_atoms.size(); a++)
        sys.mol(all_atoms.mol(a)).pos(all_atoms.atom(a)) = rot2 * all_atoms.pos(a);
    
    // finally align the com with the origin:
    Vec com = PositionUtils::com(sys, atoms);
    for(int a = 0; a < atoms.size(); a++)
        sys.mol(atoms.mol(a)).pos(atoms.atom(a)) = atoms.pos(a) - com;
//    for(int m = 0; m < sys.numMolecules(); m++)
//        for(int a = 0; a < sys.mol(m).numAtoms(); a++)
//            sys.mol(m).pos(a) = sys.mol(m).pos(a) - com;
}


Vec calc_max_vector(AtomSpecifier &as_rot, int dim) {
    // calculate the longest distance between selected atoms considering
    // the first dim dimensions.
    double max2 = 0.0, d2 = 0;
    int max_as1 = 0, max_as2 = 0;
    Vec pos_as1, pos_as2, ret;

    for(int as1 = 0; as1 < as_rot.size() - 1; as1++) {
        //for(int a1 = 0; a1 < sys.mol(m1).numAtoms(); a1++) {
        pos_as1 = as_rot.pos(as1);
        for(int as2 = as1 + 1; as2 < as_rot.size(); as2++) {
            //int start = 0;
            //if(m1 == m2) start = a1;
            //for(int a2 = start; a2 < sys.mol(m2).numAtoms(); a2++) {
            d2 = 0.0;
            pos_as2 = as_rot.pos(as2);
            for(int i = 0; i < dim; i++) {
                d2 += (pos_as1[i] - pos_as2[i]) * (pos_as1[i] - pos_as2[i]);
            }
            if(d2 > max2) {
                max2 = d2;
                max_as1 = as1;
                max_as2 = as2;
            }
            //}
        }
        //}
    }
    ret = as_rot.pos(max_as1) - as_rot.pos(max_as2);
    return ret;
}
