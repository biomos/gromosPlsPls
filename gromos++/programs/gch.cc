/**
 * @file gch.cc
 * Generate coordinates for explicit hydrogens
 */

/**
 * @page programs Program Documentation
 *
 * @anchor gch
 * @section gch Generate coordinates for explicit hydrogens
 * @author @ref co
 * @date 11-6-07
 *
 * In the standard GROMOS force fields, part of the hydrogen atoms (polar,
 * aromatic) are explicitly treated, whereas other hydrogen atoms (aliphatic,
 * some aromatic) are implicitly treated by incorporation into the (carbon)atom
 * to which they are attached. Depending on the presence or absence of hydrogen
 * atom coordinates in a molecular configuration file, hydrogen atom
 * coordinates may have to be recalculated.
 *
 * Program gch calculates optimal positions for hydrogen atoms for which the 
 * connecting bond shows a relative deviation from the zero-energy value larger
 * than a user specified threshold. Coordinates for all hydrogen atoms that are
 * explicitly listed in the topology should already be contained in the
 * coordinate file. Program @ref pdb2g96 e.g. will include atoms for which no
 * coordinates were present in the pdb-file with coordinates set to zero. If
 * defined, gch uses topological information on bonds, bond angles and dihedral
 * angles to place hydrogen atoms at the optimal location. In cases where the
 * necessary angular parameters are not provided in the topology, gch uses
 * 109.5 degrees for tetrahedral centers and 120 degrees for planar centers.
 *
 * Eight types of geometries can be handled when generating hydrogen atom 
 * coordinates:
 * <ol>
 * <li> An atom (a) is bonded to one hydrogen (H) and one other heavy atom 
 *      (nh). A fourth atom (f) is searched for which is bounded to nh and 
 *      preferably is used to define the dihedral around the nh-a bond. The 
 *      coordinates of H are generated in such a way that the dihedral 
 *      f-nh-a-H is trans and that the angle nh-a-H and bond length a-H
 *      correspond to their minimum energy values.
 * <li> An atom (a) is bonded to one hydrogen (H) and two other heavy atoms
 *      (nh1 and nh2). The coordinates of H are generated to be in the plane 
 *      through nh1, nh2 and a, on the line bisecting the nh1-a-nh2 angle and
 *      with an a-H bond length corresponding to the minimum energy value in
 *      the topology, such that the nh1-a-H and nh2-a-H angles are larger than
 *      90 degrees.
 * <li> An atom (a) is bonded to two hydrogens (H1 and H2) and one other heavy
 *      atom (nh). A fourth atom (f) is searched for which is bounded to nh
 *      and preferably is used to define the dihedral around the nh-a bond. The
 *      coordinates of H1 are generated in such a way that the dihedral 
 *      f-nh-a-H1 is trans and that the angle nh-a-H1 and bond length a-H1 
 *      correspond to their minimum energy values. The coordinates of H2 are 
 *      generated to have the angles nh-a-H2 and H1-a-H2 as well as the bond
 *      length a-H2 at their minimum energy values. If this does not result in
 *      a planar configuration around a, the improper dihedral a-nh-H1-H2 will
 *      be positive.
 * <li> An atom (a) is bonded to three hydrogens (H1, H2 and H3) and one other
 *      heavy atom (nh). A fourth atom (f) is searched for wich is bounded to
 *      nh and preferably is used to define the dihedral around the nh-a bond.
 *      The coordinates of H1 are generated in such a way that the dihedral
 *      f-nh-a-H1 is trans and that the angle nh-a-H1 and bond length a-H1
 *      correspond to their minimum energy values. The coordinates of H2 are 
 *      such that the angles nh-a-H2 and H1-a-H2 and the bond length a-H2 are 
 *      at their minimum energy values, and the improper dihedral a-nh-H1-H2 is
 *      positive. The coordinates of H3 are such that the angles nh-a-H3 and
 *      H1-a-H3 and the bond length a-H3 are at their minimum energy values and
 *      the improper dihedral a-nh-H1-H3 has a negative value.
 * <li> An atom (a) is bonded to one hydrogen atom (H) and three other heavy
 *      atoms (nh1, nh2, nh3). The coordinates of H are generated along the line
 *      going through atom a and a point corresponding to the average position
 *      of nh1, nh2 and nh3, such that the bond length a-H is at its minimum
 *      energy value and the angles nh1-a-H, nh2-a-H and nh3-a-H are larger than
 *      90 degree.
 * <li> An atom (a) is bonded to two hydrogen atoms (H1 and H2) and two other
 *      heavy atoms (nh1 and nh2). The coordinates of H1and H2 are placed above
 *      and below the plane going through atoms nh1, nh2 and a, in such a way
 *      that the a-H1 and a-H2 bond lengths and the H1-a-H2 bond angle are at
 *      their minimum energy values. The improper dihedral angle a,nh1,nh2,H1
 *      will be positive.
 * <li> An atom (a) is bonded to two hydrogen atoms (H1 and H2), but to no
 *      heavy atoms. This is likely to be a (crystallographic) water molecule. 
 *      First a molecule is generated having the a-H1 aligned in the 
 *      z-direction and the a-H2 in the z-y plane with the angle H1-a-H2 and
 *      bond lengths a-H1 and a-H2 according to their minimum energy values. 
 *      This molecule is then rotated around x, y and z by three random angles.
 * <li> An atom (a) is bonded to four hydrogen atoms (H1, H2, H3 and H4), but
 *      to no heavy atoms. A molecule is generated with all bond lengths at 
 *      their minimum energy value, the a-H1 aligned in the z-directionm H2 in
 *      the x-z plane and H3 such that the improper a-H1-H2-H3 is positive and 
 *      H4 such that the improper a-H1-H2-H4 is negative. The complete molecule
 *      is then rotated by three random angles around x, y and z.
 * </ol>
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@pos</td><td>&lt;input coordinate file&gt; </td></tr>
 * <tr><td> \@tol</td><td>&lt;tolerance (default 0.1 %)&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;boundary type&gt; &lt;gather method&gt; </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  gch
    @topo  ex.top
    @pos   exref.coo
    @tol   0.1
    @pbc   r
 @endverbatim
 *
 * <hr>
 */


#include <cassert>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/bound/Boundary.h"
#include "../src/gio/InG96.h"
#include "../src/gio/OutG96S.h"
#include "../src/gcore/AtomPair.h"
#include "../src/gcore/GromosForceField.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/LJException.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/Solvent.h"
#include "../src/gcore/SolventTopology.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gcore/Bond.h"
#include "../src/gcore/BondType.h"
#include "../src/gcore/Constraint.h"
#include "../src/gcore/Angle.h"
#include "../src/gcore/AngleType.h"
#include "../src/gcore/Dihedral.h"
#include "../src/gio/InTopology.h"
#include "../src/utils/Neighbours.h"
#include "../src/gmath/Physics.h"
#include "../src/gmath/Vec.h"
#include "../src/gmath/Matrix.h"
#include <vector>
#include <iomanip>
#include <fstream>
#include <iostream>
#include <sstream>

using namespace std;
using namespace gcore;
using namespace gio;
using namespace args;
using namespace utils;
using namespace gmath;
using namespace bound;



int generate_coordinate(System *sys, GromosForceField *gff, int m, int a,
        vector<int> h, vector<int> nh, int geom, double eps);

double find_bond(System *sys, GromosForceField *gff, int m, Bond b, double guess);

double find_angle(System *sys, GromosForceField *gff, int m, Angle a, double guess);
int find_dihedral(System *sys, int m, int i, int j, vector<int> h);

int main(int argc, char **argv) {

  Argument_List knowns;
  knowns << "topo" << "pos" << "tol" << "pbc";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo  <molecular topology file>\n";
  usage += "\t@pos   <input coordinate file>\n";
  usage += "\t@tol   <tolerance (default 0.1 %)>\n";
  usage += "\t@pbc   <boundary type> <gather method>\n";


  try {
    Arguments args(argc, argv, knowns, usage);

    // read topology make system and force field
    InTopology it(args["topo"]);
    System sys(it.system());
    GromosForceField gff(it.forceField());

    System refSys(it.system());

    // read in coordinates
    InG96 ic(args["pos"]);
    ic.select("ALL");
    ic >> sys;
    ic.close();
    System outSys(sys);

    // read in the accuracy
    double eps = args.getValue<double>("tol", false, 0.1) / 100.0;

    // parse boundary conditions
    Boundary *pbc = BoundaryParser::boundary(sys, args);
    // parse gather method
    Boundary::MemPtr gathmethod = args::GatherParser::parse(sys, refSys, args);

    // gather the system!
    (*pbc.*gathmethod)();

    // a bit of an ugly hack for solvent molecules. The whole program has
    // been written for solutes, but sometimes we would have crystallographic
    // waters for which we want to do the same thing.
    if (sys.sol(0).numPos()) {
      MoleculeTopology mt;
      for (int a = 0; a < sys.sol(0).topology().numAtoms(); a++) {
        mt.addAtom(sys.sol(0).topology().atom(a));
        mt.setResNum(a, 0);
      }
      mt.setResName(0, "SOLV");

      ConstraintIterator ci(sys.sol(0).topology());
      for (; ci; ++ci) {
        gff.addBondType(BondType(1, ci().dist()));
        Bond b(ci()[0], ci()[1], false);
        b.setType(gff.numBondTypes() - 1);
        mt.addBond(b);
      }

      // add every solvent molecule as a solute
      int numSolvent = sys.sol(0).numPos() / sys.sol(0).topology().numAtoms();
      for (int i = 0; i < numSolvent; i++) {
        Molecule m(mt);
        m.initPos();
        for (int a = 0; a < mt.numAtoms(); a++) {
          m.pos(a) = sys.sol(0).pos(i * sys.sol(0).topology().numAtoms() + a);
        }
        sys.addMolecule(m);
      }
      // and remove the original solvent
      sys.sol(0).setNumPos(0);
    }

    int dummyType = it.forceField().dummyAtomType();
    // initialize two counters
    int replaced = 0, kept = 0;

    // loop over all atoms
    for (int m = 0; m < sys.numMolecules(); m++) {

      // flag the atoms with mass 1.008 as hydrogens
      sys.mol(m).topology().setHmass(1.008);
      for (int a = 0; a < sys.mol(m).numAtoms(); a++) {

        // find the neighbours of this atom
        utils::Neighbours n(sys, m, a);

        // divide into hydrogens and non-hydrogens
        vector<int> h;
        vector<int> nh;

        for (unsigned int i = 0; i < n.size(); i++) {
          if (sys.mol(m).topology().atom(n[i]).isH()) {
            h.push_back(n[i]);
          } else {
            // only add it if is not a dummy atom
            if (dummyType != -1 &&
                    sys.mol(m).topology().atom(n[i]).iac() != dummyType)
              nh.push_back(n[i]);
          }
        }


        // determine what kind of geometry they should be
        int geom = 0;

        // only continue if we have hydrogens
        int numH = h.size();
        int numNH = nh.size();

        if (numH == 1 && numNH == 1) geom = 1;
        if (numH == 1 && numNH == 2) geom = 2;
        if (numH == 2 && numNH == 1) geom = 3;
        if (numH == 3 && numNH == 1) geom = 4;
        // crystallographic water
        if (numH == 2 && numNH == 0) geom = 5;
        // nh4+
        if (numH == 4 && numNH == 0) geom = 6;
        if (numH == 1 && numNH == 3) geom = 7;
        if (numH == 2 && numNH == 2) geom = 8;

        if (numH && !geom) {
          ostringstream os;
          os << "Unexpected geometry for hydrogen bound to atom: "
                  << m + 1 << ":" << a + 1 << endl;
          throw (gromos::Exception("gch", os.str()));
        }
        // we have to have a geometry (this means that there are hydrogens)
        // and a should not be a hydrogen itself. (in the case of H2O we have
        // e.g. H-H bonds. These are only treated via the oxygen.
        if (geom && !sys.mol(m).topology().atom(a).isH()) {
          int r = generate_coordinate(&sys, &gff, m, a, h, nh, geom, eps);
          replaced += r;
          kept += (numH - r);
        }
      }
    }

    // fix the solvent hack
    int solventIndex = 0;
    for (int m = 0; m < sys.numMolecules(); ++m) {
      if (m < outSys.numMolecules()) {
        // solute
        for (int a = 0; a < sys.mol(m).numAtoms(); ++a) {
          outSys.mol(m).pos(a) = sys.mol(m).pos(a);
        }
      } else {
        // solvent
        for (int a = 0; a < sys.mol(m).numAtoms(); ++a, ++solventIndex) {
          outSys.sol(0).pos(solventIndex) = sys.mol(m).pos(a);
        }
      }
    }

    OutG96S oc(cout);
    oc.select("ALL");
    ostringstream os;
    os << "gch found " << replaced + kept << " hydrogen atoms in "
            << args["pos"] << endl;
    os << kept << " were within " << eps * 100 << "% of minimum energy bond "
            << "length" << endl;
    os << replaced << " were assigned new coordinates based on geometry";


    oc.writeTitle(os.str());

    oc << outSys;
    oc.close();

  }  catch (const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);

  }
  return 0;
}

int generate_coordinate(System *sys, GromosForceField *gff, int m, int a,
        vector<int> h, vector<int> nh, int geom, double eps) {
  int count = 0;
  switch (geom) {
    case(1):
    {

      double bond = find_bond(sys, gff, m, Bond(a, h[0], false), 0.1);
      Vec v0 = sys->mol(m).pos(a) - sys->mol(m).pos(h[0]);
      if (fabs(v0.abs() - bond) / bond > eps) {

        double angle = M_PI / 180.0 * find_angle(sys, gff, m,
                Angle(nh[0], a, h[0]), 109.5);

        // in case of ARGN, HH1 will be in cis-conformation instead of anti
        // (dihedral NE-CZ-NH1-HH1)
        // Claudio M. Soares, private communication, May 2011
        if (sys->mol(m).topology().resName(sys->mol(m).topology().resNum(nh[0])) == "ARGN") {
          angle *= -1;
        }

        int fourth = find_dihedral(sys, m, nh[0], a, h);

        Vec v1 = (sys->mol(m).pos(fourth) - sys->mol(m).pos(nh[0])).normalize();
        Vec v2 = (sys->mol(m).pos(nh[0]) - sys->mol(m).pos(a)).normalize();
        Vec v4 = (v1.cross(v2)).normalize();
        Vec v5 = (v2.cross(v4)).normalize();
        Vec v6 = bond * cos(angle) * v2 - bond * sin(angle) * v5;

        sys->mol(m).pos(h[0]) = sys->mol(m).pos(a) + v6;
        count++;
      }
      break;
    }

    case(2):
    {
      double bond = find_bond(sys, gff, m, Bond(a, h[0], false), 0.1);
      Vec v0 = sys->mol(m).pos(a) - sys->mol(m).pos(h[0]);
      if (fabs(v0.abs() - bond) / bond > eps) {

        Vec v1 = (sys->mol(m).pos(nh[0]) - sys->mol(m).pos(a)).normalize();
        Vec v2 = (sys->mol(m).pos(nh[1]) - sys->mol(m).pos(a)).normalize();
        Vec v3 = (v1 + v2).normalize();
        sys->mol(m).pos(h[0]) = sys->mol(m).pos(a) - bond*v3;
        count++;
      }
      break;
    }
    case(3):
    {
      double bond1 = find_bond(sys, gff, m, Bond(a, h[0], false), 0.1);
      double bond2 = find_bond(sys, gff, m, Bond(a, h[1], false), 0.1);
      Vec v01 = sys->mol(m).pos(a) - sys->mol(m).pos(h[0]);
      Vec v02 = sys->mol(m).pos(a) - sys->mol(m).pos(h[1]);
      if (fabs(v01.abs() - bond1) / bond1 > eps ||
              fabs(v02.abs() - bond2) / bond2 > eps) {

        double angle1 = M_PI / 180.0 * find_angle(sys, gff, m,
                Angle(nh[0], a, h[0], false), 120);
        double angle2 = M_PI / 180.0 * find_angle(sys, gff, m,
                Angle(nh[0], a, h[1], false), 120);
        double angle3 = M_PI / 180.0 * find_angle(sys, gff, m,
                Angle(h[0], a, h[1], false), 120);
        int fourth = find_dihedral(sys, m, nh[0], a, h);

        Vec v1 = (sys->mol(m).pos(fourth) - sys->mol(m).pos(nh[0])).normalize();
        Vec v2 = (sys->mol(m).pos(nh[0]) - sys->mol(m).pos(a)).normalize();
        Vec v4 = v1.cross(v2).normalize();
        Vec v5 = v2.cross(v4).normalize();
        Vec v6 = bond1 * cos(angle1) * v2
                - bond1 * sin(angle1) * v5;
        double A = bond2 * cos(angle2);
        double B = (bond1 * bond2 * cos(angle3) - A * v2.dot(v6)) / v5.dot(v6);
        double C = sqrt(bond2 * bond2 - A * A - B * B);
        Vec v7 = A * v2 + B * v5 + C * v4;

        if (fabs(v01.abs() - bond1) / bond1 > eps) {
          sys->mol(m).pos(h[0]) = sys->mol(m).pos(a) + v6;
          count++;
        }
        if (fabs(v02.abs() - bond2) / bond2 > eps) {
          sys->mol(m).pos(h[1]) = sys->mol(m).pos(a) + v7;
          count++;
        }
      }
      break;
    }
    case(4):
    {
      // very similar to the non-planar type of above
      double bond1 = find_bond(sys, gff, m, Bond(a, h[0], false), 0.1);
      double bond2 = find_bond(sys, gff, m, Bond(a, h[1], false), 0.1);
      double bond3 = find_bond(sys, gff, m, Bond(a, h[2], false), 0.1);
      Vec v01 = sys->mol(m).pos(a) - sys->mol(m).pos(h[0]);
      Vec v02 = sys->mol(m).pos(a) - sys->mol(m).pos(h[1]);
      Vec v03 = sys->mol(m).pos(a) - sys->mol(m).pos(h[2]);

      if (fabs(v01.abs() - bond1) / bond1 > eps ||
              fabs(v02.abs() - bond2) / bond2 > eps ||
              fabs(v03.abs() - bond3) / bond3 > eps) {

        double angle1 = M_PI / 180.0 * find_angle(sys, gff, m,
                Angle(nh[0], a, h[0], false), 109.5);
        double angle2 = M_PI / 180.0 * find_angle(sys, gff, m,
                Angle(nh[0], a, h[1], false), 109.5);
        double angle3 = M_PI / 180.0 * find_angle(sys, gff, m,
                Angle(nh[0], a, h[2], false), 109.5);
        double angle4 = M_PI / 180.0 * find_angle(sys, gff, m,
                Angle(h[0], a, h[1], false), 109.5);
        double angle5 = M_PI / 180.0 * find_angle(sys, gff, m,
                Angle(h[0], a, h[2], false), 109.5);
        int fourth = find_dihedral(sys, m, nh[0], a, h);

        Vec v1 = (sys->mol(m).pos(fourth) - sys->mol(m).pos(nh[0])).normalize();
        Vec v2 = (sys->mol(m).pos(nh[0]) - sys->mol(m).pos(a)).normalize();
        Vec v3 = (sys->mol(m).pos(a) - sys->mol(m).pos(h[0])).normalize();
        Vec v4 = v2.cross(v1);
        if(v4.abs() < 1e-10) {
          if(v2[0] != 0.0 || v2[1] != 0.0) {
            v4[0] = -v2[1];
            v4[1] = v2[0];
            v4[2] = 0;
          } else {
            v4[0] = 0;
            v4[1] = -v2[2];
            v4[2] = v2[1];
          }
          v4 = v4.normalize();
        } else {
          v4 = v4.normalize();
        }
        Vec v5 = v2.cross(v4).normalize();
        Vec v6 = bond1 * cos(angle1) * v2 - bond1 * sin(angle1) * v5;

        double A = bond2 * cos(angle2);
        double B = (bond1 * bond2 * cos(angle4) - A * v2.dot(v6)) / v5.dot(v6);
        double C = sqrt(bond2 * bond2 - A * A - B * B);
        Vec v7 = A * v2 + B * v5 + C * v4;
        A = bond3 * cos(angle3);
        B = (bond1 * bond3 * cos(angle5) - A * v2.dot(v6)) / v5.dot(v6);
        C = sqrt(bond3 * bond3 - A * A - B * B);
        Vec v8 = A * v2 + B * v5 - C * v4;
        if (fabs(v01.abs() - bond1) / bond1 > eps) {
          sys->mol(m).pos(h[0]) = sys->mol(m).pos(a) + v6;
          count++;
        }
        if (fabs(v02.abs() - bond2) / bond2 > eps) {
          sys->mol(m).pos(h[1]) = sys->mol(m).pos(a) + v7;
          count++;
        }
        if (fabs(v03.abs() - bond3) / bond3 > eps) {
          sys->mol(m).pos(h[2]) = sys->mol(m).pos(a) + v8;
          count++;
        }
      }
      break;
    }
    case(5):
    {
      // likely to be a water molecule. Here we have to come up with some
      // random orientation.
      double bond1 = find_bond(sys, gff, m, Bond(a, h[0], false), 0.1);
      double bond2 = find_bond(sys, gff, m, Bond(a, h[1], false), 0.1);
      Vec v01 = sys->mol(m).pos(a) - sys->mol(m).pos(h[0]);
      Vec v02 = sys->mol(m).pos(a) - sys->mol(m).pos(h[1]);
      if (fabs(v01.abs() - bond1) / bond1 > eps ||
              fabs(v02.abs() - bond2) / bond2 > eps) {
        // first generate a standard molecule. If it really is a water 
        // molecule, there is probably no angle defined, but putting them 
        // at 109.5 degrees is not so bad.
        double angle = M_PI / 180.0 * find_angle(sys, gff, m,
                Angle(h[0], a, h[1], false), 109.5);
        Vec v1(0.0, 0.0, bond1);
        Vec v2(0.0, bond2 * sin(angle), bond2 * cos(angle));

        //get three random numbers for the angles
        //calculate sin and cos of these angles

        Vec angle_cos, angle_sin;
        for (int i = 0; i < 3; i++) {
          double ang = 2.0 * M_PI / RAND_MAX * double(rand());
          angle_cos[i] = cos(ang);
          angle_sin[i] = sin(ang);
        }

        // prepare a matrix to perform three random rotations  
        // The product of three rotation matrices about three axes
        /*
         * (  1.0   0.0   0.0)   ( cosy   0.0   siny)   ( cosx  -sinx   0.0)
         * (  0.0  cosz -sinz) X (  0.0   1.0    0.0) X ( sinx   cosx   0.0)
         * (  0.0  sinz  cosz)   (-siny   0.0   cosy)   (  0.0    0.0   1.0)
         */
        gmath::Matrix rot(3, 3);
        rot(0, 0) = angle_cos[0] * angle_cos[1];
        rot(1, 0) = angle_sin[0] * angle_cos[2]
                + angle_cos[0] * angle_sin[1] * angle_sin[2];
        rot(2, 0) = angle_sin[0] * angle_sin[2]
                - angle_cos[0] * angle_sin[1] * angle_cos[2];
        rot(0, 1) = -angle_sin[0] * angle_cos[1];
        rot(1, 1) = angle_cos[0] * angle_cos[2]
                - angle_sin[0] * angle_sin[1] * angle_sin[2];
        rot(2, 1) = angle_cos[0] * angle_sin[2]
                + angle_sin[0] * angle_sin[1] * angle_cos[2];
        rot(0, 2) = angle_sin[1];
        rot(1, 2) = -angle_cos[1] * angle_sin[2];
        rot(2, 2) = angle_cos[1] * angle_cos[2];

        // rotate the hydrogens and put the coordinates
        if (fabs(v01.abs() - bond1) / bond1 > eps) {
          sys->mol(m).pos(h[0]) = sys->mol(m).pos(a) + rot*v1;
          count++;
        }
        if (fabs(v02.abs() - bond2) / bond2 > eps) {
          sys->mol(m).pos(h[1]) = sys->mol(m).pos(a) + rot*v2;
          count++;
        }
      }
      break;
    }
    case(6):
    {
      // nh4+, simliar to case 5
      double bond1 = find_bond(sys, gff, m, Bond(a, h[0], false), 0.1);
      double bond2 = find_bond(sys, gff, m, Bond(a, h[1], false), 0.1);
      double bond3 = find_bond(sys, gff, m, Bond(a, h[2], false), 0.1);
      double bond4 = find_bond(sys, gff, m, Bond(a, h[3], false), 0.1);
      Vec v01 = sys->mol(m).pos(a) - sys->mol(m).pos(h[0]);
      Vec v02 = sys->mol(m).pos(a) - sys->mol(m).pos(h[1]);
      Vec v03 = sys->mol(m).pos(a) - sys->mol(m).pos(h[2]);
      Vec v04 = sys->mol(m).pos(a) - sys->mol(m).pos(h[3]);
      if (fabs(v01.abs() - bond1) / bond1 > eps ||
              fabs(v02.abs() - bond2) / bond2 > eps ||
              fabs(v03.abs() - bond3) / bond3 > eps ||
              fabs(v04.abs() - bond4) / bond4 > eps) {

        // no angle search, theta is 109.5
        // phi is 0, 120, 240
        double angle = M_PI / 180.0 * 109.5;
        Vec v1(0.0, 0.0, bond1);
        Vec v2(bond2 * sin(angle) * cos(0.0),
                bond2 * sin(angle) * sin(0.0),
                bond2 * cos(angle));
        Vec v3(bond3 * sin(angle) * cos(M_PI / 180.0 * 120.0),
                bond3 * sin(angle) * sin(M_PI / 180.0 * 120.0),
                bond3 * cos(angle));
        Vec v4(bond4 * sin(angle) * cos(M_PI / 180.0 * 240.0),
                bond4 * sin(angle) * sin(M_PI / 180.0 * 240.0),
                bond4 * cos(angle));

        //get three random numbers for the angles
        //calculate sin and cos of these angles

        Vec angle_cos, angle_sin;
        for (int i = 0; i < 3; i++) {
          double ang = 2.0 * M_PI / RAND_MAX * double(rand());
          angle_cos[i] = cos(ang);
          angle_sin[i] = sin(ang);
        }

        gmath::Matrix rot(3, 3);
        rot(0, 0) = angle_cos[0] * angle_cos[1];
        rot(1, 0) = angle_sin[0] * angle_cos[2]
                + angle_cos[0] * angle_sin[1] * angle_sin[2];
        rot(2, 0) = angle_sin[0] * angle_sin[2]
                - angle_cos[0] * angle_sin[1] * angle_cos[2];
        rot(0, 1) = -angle_sin[0] * angle_cos[1];
        rot(1, 1) = angle_cos[0] * angle_cos[2]
                - angle_sin[0] * angle_sin[1] * angle_sin[2];
        rot(2, 1) = angle_cos[0] * angle_sin[2]
                + angle_sin[0] * angle_sin[1] * angle_cos[2];
        rot(0, 2) = angle_sin[1];
        rot(1, 2) = -angle_cos[1] * angle_sin[2];
        rot(2, 2) = angle_cos[1] * angle_cos[2];



        // rotate the hydrogens and put the coordinates
        if (fabs(v01.abs() - bond1) / bond1 > eps) {
          sys->mol(m).pos(h[0]) = sys->mol(m).pos(a) + rot*v1;
          count++;
        }
        if (fabs(v02.abs() - bond2) / bond2 > eps) {
          sys->mol(m).pos(h[1]) = sys->mol(m).pos(a) + rot*v2;
          count++;
        }
        if (fabs(v03.abs() - bond3) / bond3 > eps) {
          sys->mol(m).pos(h[2]) = sys->mol(m).pos(a) + rot*v3;
          count++;
        }
        if (fabs(v04.abs() - bond4) / bond4 > eps) {
          sys->mol(m).pos(h[3]) = sys->mol(m).pos(a) + rot*v4;
          count++;
        }
      }
      break;
    }
    case(7):
    {
      // charged NH, connected to 3 NH atoms.
      double bond = find_bond(sys, gff, m, Bond(a, h[0], false), 0.1);
      Vec v0 = sys->mol(m).pos(a) - sys->mol(m).pos(h[0]);
      if (fabs(v0.abs() - bond) / bond > eps) {

        Vec v1 = sys->mol(m).pos(nh[0]) - sys->mol(m).pos(a);
        Vec v2 = sys->mol(m).pos(nh[1]) - sys->mol(m).pos(a);
        Vec v3 = sys->mol(m).pos(nh[2]) - sys->mol(m).pos(a);
        Vec v4 = (v1 + v2 + v3).normalize();
        sys->mol(m).pos(h[0]) = sys->mol(m).pos(a) - bond*v4;
        count++;
      }
      break;
    }
    case(8):
    {
      // charged NH2, connected to 2 NH atoms
      double bond1 = find_bond(sys, gff, m, Bond(a, h[0], false), 0.1);
      double bond2 = find_bond(sys, gff, m, Bond(a, h[1], false), 0.1);
      Vec v0 = sys->mol(m).pos(a) - sys->mol(m).pos(h[0]);
      Vec v1 = sys->mol(m).pos(a) - sys->mol(m).pos(h[1]);
      if (fabs(v0.abs() - bond1) / bond1 > eps ||
              fabs(v1.abs() - bond2) / bond2 > eps) {

        Vec v2 = sys->mol(m).pos(nh[0]) - sys->mol(m).pos(a);
        Vec v3 = sys->mol(m).pos(nh[1]) - sys->mol(m).pos(a);
        Vec v4 = -(v2 + v3).normalize();
        Vec v5 = v2.cross(v3).normalize();
        double angle = M_PI / 180.0 * find_angle(sys, gff, m,
                Angle(h[0], a, h[1], false), 109.5);

        sys->mol(m).pos(h[0]) = sys->mol(m).pos(a) +
                bond1 * sin(0.5 * angle) * v5 +
                bond1 * cos(0.5 * angle) * v4;
        sys->mol(m).pos(h[1]) = sys->mol(m).pos(a) -
                bond2 * sin(0.5 * angle) * v5 +
                bond2 * cos(0.5 * angle) * v4;
        count++;
      }
      break;
    }


  }
  return count;
}

double find_bond(System *sys, GromosForceField *gff, int m, Bond b, double guess) {
  BondIterator bi(sys->mol(m).topology());
  double value = 0.0;
  for (; bi; ++bi) {
    if (bi()[0] == b[0] && bi()[1] == b[1]) {
      value = gff->bondType(bi().type()).b0();
      break;
    }
  }
  if (value != 0.0)
    return value;
  else
    return guess;
}

double find_angle(System *sys, GromosForceField *gff, int m, Angle a, double guess) {
  AngleIterator ai(sys->mol(m).topology());
  double value = 0.0;

  for (; ai; ++ai) {
    if (ai()[0] == a[0] && ai()[1] == a[1] && ai()[2] == a[2]) {
      value = gff->angleType(ai().type()).t0();
      break;
    }
  }
  if (value != 0.0)
    return value;
  else
    return guess;
}

int find_dihedral(System *sys, int m, int i, int j, vector<int> h) {

  if (j < i) {
    int t = i;
    i = j;
    j = t;
  }

  DihedralIterator di(sys->mol(m).topology());
  int fourth = -1;

  for (; di; ++di) {
    if (di()[1] == i && di()[2] == j) {
      for (unsigned int k = 0; k < h.size(); k++) {
        if (di()[0] == h[k]) fourth = di()[3];
        if (di()[3] == h[k]) fourth = di()[0];
      }
      if (fourth != -1) break;
    }
  }
  if (fourth == -1) {

    // cannot find a dihedral, then just take the first atom bounded to i, whicih is not a hydrogen
    Neighbours n(*sys, m, i);
    for (unsigned int k = 0; k < n.size(); k++) {
      int hydrogen = 0;
      for (unsigned int l = 0; l < h.size(); l++)
        if (n[k] == h[l]) hydrogen = 1;
      if (!hydrogen && n[k] != j) fourth = n[k];
    }
    if (fourth == -1) {
      Neighbours n(*sys, m, j);
      for (unsigned int k = 0; k < n.size(); k++) {
        int hydrogen = 0;
        for (unsigned int l = 0; l < h.size(); l++)
          if (n[k] == h[l]) hydrogen = 1;
        if (!hydrogen && n[k] != i) fourth = n[k];
      }
    }
    if (fourth == -1)
      throw (gromos::Exception("find_dihedral", "undefined position"));
  }

  return fourth;

}

