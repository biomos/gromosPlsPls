/*
 * This file is part of GROMOS.
 * 
 * Copyright (c) 2011, 2012, 2016, 2018, 2021, 2023 Biomos b.v.
 * See <https://www.gromos.net> for details.
 * 
 * GROMOS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

// gio_OutCif.cc

#include <iostream>
#include <iomanip>
#include "OutCif.h"
#include "../gromos/Exception.h"
#include "../gcore/System.h"
#include "../gcore/Molecule.h"
#include "../gcore/Solvent.h"
#include "../gmath/Vec.h"
#include "../gcore/Box.h"
#include "../utils/AtomSpecifier.h"

#include "../gcore/Remd.h"
#include "../args/Arguments.h"
#include "../gmath/Matrix.h"

using namespace gio;
using namespace gcore;
using namespace utils;

class gio::OutCif_i {
    friend class gio::OutCif;
    ostream &d_os;
    int d_count, d_res_off, d_switch;
    double d_factor;

    OutCif_i ( ostream &os ) :
        d_os{os}, d_count{0}, d_res_off{1}, d_switch{0}, d_factor{10.0} {
        cout << "I'm OutCIF_i constructor\n";
    }

    ~OutCif_i () {
        cout << "I'm OutCIF_i destructor\n";
    }

    void writeSingleM ( const Molecule &mol, const int mn );
    void writeSingleV ( const gcore::System &sys );
    void writeSingleS ( const Solvent &sol );
    void writeCell ( const Box &box );
    void writeBox ( const Box &box );
    void writeTriclinicBox ( const Box &box );
    void writeGenBox ( const Box &box );
    void writeAtomSpecifier ( const AtomSpecifier & atoms );
};

OutCif::OutCif() :
    OutCoordinates(), d_this{nullptr} {
    cout << "I'm OutCIF default constructor\n";
}

OutCif::OutCif( ostream &os ) :
    OutCoordinates(), d_this{new OutCif_i(os)} {
    cout << "I'm OutCIF default constructor\n";
}

OutCif::~OutCif() {
    if ( d_this ) delete d_this;
    cout << "I'm OutCIF destructor\n";
}

void OutCif::select ( const string &thing ) {
    if (thing == "ALL") {
        d_this->d_switch = 1;
    } else if (thing == "SOLVENT") {
        d_this->d_switch = 2;
    } else if (thing == "SOLUTEV") {
        d_this->d_switch = 3;
    } else if (thing == "ALLV") {
        d_this->d_switch = 4;
    } else if (thing == "SOLVENTV") {
        d_this->d_switch = 5;
    } else {
        d_this->d_switch = 0;
    }
}

void OutCif::open ( ostream &os ) {
    if ( d_this ) delete d_this;

    d_this = new OutCif_i( os );
}

void OutCif::close() {
    if ( d_this ) delete d_this;

    d_this = nullptr;
}

void OutCif::writeTitle ( const string &title ) {
    d_this->d_os << "_struct.title\n" << title << "\n_end\n";
}

void OutCif::writeTimestep ( const int step, const double time ) {
    d_this->d_os.precision(9);
    d_this->d_os.setf( ios::fixed, ios::floatfield );

    d_this->d_os << "_timestep\n"
          << std::setw(18)
          << step
          << std::setw(20)
          << time
          //<< "\n#if @time flag is used the value for step refers to the"
          //<< "\n#step-th configuration in the original trajectory file"          
          << "\n_end\n";
}

OutCif& OutCif::operator<< ( const gcore::System &sys ) {

    d_this->d_os << "_model\n";
    d_this->d_count = 0;
    d_this->d_res_off = 1;
    if ( d_this->d_switch == 0 || d_this->d_switch == 1 || 
         d_this->d_switch == 3 || d_this->d_switch == 4 ) {
        for ( auto i = 0; i != sys.numMolecules(); ++i ) {
            d_this->writeSingleM( sys.mol(i), i+1 );
        }
    }
    if ( d_this->d_switch == 3 || d_this->d_switch == 4 || 
         d_this->d_switch == 5 ) {
        d_this->writeSingleV( sys );
    }
    if ( d_this->d_switch == 1 || d_this->d_switch == 2 ||
         d_this->d_switch == 4 || d_this->d_switch == 5) {
        for ( auto i = 0; i != sys.numSolvents(); ++i ) {
            d_this->writeSingleS( sys.sol(i) );
        }
    }

    //TODO: write connect

    return *this;
}

OutCif& OutCif::operator<< ( const utils::AtomSpecifier &atoms ) {

    d_this->d_os << "_model\n";

    // TODO: write box cryst
    d_this->writeAtomSpecifier( atoms );

    return *this;
}

void gio::OutCif_i::writeCell ( const Box &box ) {
    if ( box.ntb() == Box::vacuum || box.ntb() == Box::truncoct ) {
        return;
    }

    ++d_count;
    d_os.setf( ios::unitbuf );
    d_os.setf( ios::left, ios::adjustfield );
    d_os << setw(6) << "_cryst1";
    d_os.setf( ios::fixed | ios::right );
    d_os.precision(3);
    d_os << setw(9) << box.K().abs()*d_factor<< setw(9) << box.L().abs()*d_factor << setw(9) << box.M().abs()*d_factor;
    d_os.precision(2);
    d_os << setw(7) << box.alpha() << setw(7) << box.beta() << setw(7) << box.gamma();
    d_os.setf( ios::left, ios::adjustfield );
    d_os << setw(11) << " P 1";
    d_os.setf( ios::fixed | ios::right );
    d_os << setw(4) << 1;
    d_os << endl;
}

void gio::OutCif_i::writeBox ( const Box &box ) {
    d_os.setf(ios::fixed, ios::floatfield);
    d_os.precision(9);

    d_os << setw(15) << box.K()[0]
        << setw(15) << box.L()[1]
        << setw(15) << box.M()[2] << endl;
}

void gio::OutCif_i::writeTriclinicBox(const Box &box) {
    d_os.setf(ios::fixed, ios::floatfield);
    d_os.precision(9);

    d_os << setw(8) << box.ntb() << endl;
    for (int i = 0; i < 3; ++i) {
        d_os << setw(15) << box.K()[i]
            << setw(15) << box.L()[i]
            << setw(15) << box.M()[i] << endl;
    }

}

void gio::OutCif_i::writeSingleM(const Molecule &mol, const int mn) {
    d_os.setf(ios::fixed, ios::floatfield);
    d_os.setf(ios::unitbuf);
    d_os.precision(3);

    for ( auto i = 0; i != mol.numAtoms(); ++i ) {
        ++d_count;
        int res = mol.topology().resNum(i);

        d_os << "_atom";
        d_os.setf( ios::right, ios::adjustfield );
        d_os << setw(7) << d_count;
        d_os.setf( ios::left, ios::adjustfield );
        if ( mol.topology().atom(i).name().length() == 4 ) {
            d_os << ' ' << setw(5) << mol.topology().atom(i).name().c_str();
        } else {
            d_os << ' ' << setw(4) << mol.topology().atom(i).name().c_str();
        }
        d_os << setw(4) << mol.topology().resName(res).c_str();

        //TODO: write chains

        d_os.setf( ios::right, ios::adjustfield );
        d_os << setw(4) << res + d_res_off << "    "
             << setw(8) << mol.pos(i)[0]*d_factor
             << setw(8) << mol.pos(i)[1]*d_factor
             << setw(8) << mol.pos(i)[2]*d_factor;
    }
    d_os << "_ter\n";
    d_res_off += mol.topology().numRes();
}

void gio::OutCif_i::writeSingleV ( const gcore::System &sys ) {
    d_os.setf( ios::fixed, ios::floatfield );
    d_os.setf( ios::unitbuf );
    d_os.precision(3);

    int res = 0;
    int mn = 1;
    if ( d_count != 0 ) {
        mn = sys.numMolecules();
    }

    for ( auto i = 0; i != sys.vas().numVirtualAtoms(); ++i ) {
        ++d_count;

        d_os << "_atom";
        d_os.setf( ios::right, ios::adjustfield );
        d_os << setw(7) << d_count;
        d_os.setf( ios::left, ios::adjustfield );
        d_os << ' ' << setw(5) << "_virt";
        d_os << setw(4) << "_vrt";

        //TODO: write chain

        d_os.setf( ios::right, ios::adjustfield );
        d_os << setw(4) << res + d_res_off
             << setw(8) << sys.vas().atom(i).pos()[0]*d_factor
             << setw(8) << sys.vas().atom(i).pos()[1]*d_factor
             << setw(8) << sys.vas().atom(i).pos()[2]*d_factor;
        }
    d_os << "_ter\n";
}

void gio::OutCif_i::writeSingleS ( const Solvent &sol ) {

    int na = sol.topology().numAtoms();
    d_os.setf( ios::fixed, ios::floatfield );
    d_os.setf( ios::unitbuf );
    d_os.precision(3);
    for ( auto i = 0; i != sol.numPos(); ++i ) {
        ++d_count;
        int res = i / na;
        int nameid = i % na;

        d_os << "_atom";
        d_os.setf( ios::right, ios::adjustfield );
        d_os << setw(7) << d_count;
        d_os.setf( ios::left, ios::adjustfield );
        if ( sol.topology().atom(nameid).name().length() == 4 ) {
            d_os << ' ' << setw(5) << sol.topology().atom(nameid).name().substr(0, 4).c_str();
        } else {
            d_os << ' ' << setw(4) << sol.topology().atom(nameid).name().substr(0, 3).c_str();
        }
        d_os << setw(4) << sol.topology().solvName().c_str();
        d_os.setf( ios::right, ios::adjustfield );
        d_os << setw(5) << res + 1 << "    "
             << setw(8) << sol.pos(i)[0]*d_factor
             << setw(8) << sol.pos(i)[1]*d_factor
             << setw(8) << sol.pos(i)[2]*d_factor;
    }
    d_os << "_ter\n";
    d_res_off += sol.numPos()/na;
}

void gio::OutCif_i::writeGenBox(const Box &box) {
    d_os.setf(ios::fixed, ios::floatfield);
    d_os.precision(9);
    const double k = box.K().abs();
    const double l = box.L().abs();
    const double m = box.M().abs();
    d_os << setw(8) << box.ntb() << endl;
    if (box.ntb() == gcore::Box::vacuum) {
        d_os << setw(15) << 0.0 << setw(15) << 0.0 << setw(15) << 0.0 << endl
            << setw(15) << 0.0 << setw(15) << 0.0 << setw(15) << 0.0 << endl
            << setw(15) << 0.0 << setw(15) << 0.0 << setw(15) << 0.0 << endl
            << setw(15) << box.X() << setw(15) << box.Y() << setw(15) << box.Z() << endl;
    } else {
        d_os << setw(15) << k
            << setw(15) << l
            << setw(15) << m << endl;
        d_os << setw(15) << acos(box.L().dot(box.M()) / (l * m))*180 / M_PI
            << setw(15) << acos(box.K().dot(box.M()) / (k * m))*180 / M_PI
            << setw(15) << acos(box.K().dot(box.L()) / (k * l))*180 / M_PI << endl;

        // calculate the Euler rotation angles as described in Phils manuscript:
        // "GROMOS01: Description of the changes", Philippe Huenenberger, October 5, 2004
        gmath::Vec x = box.K().normalize();
        gmath::Vec y = (box.L() - (box.L().dot(x) * x)).normalize();
        gmath::Vec z = x.cross(y);

        gmath::Matrix R_(x, y, z);
        double R11R21 = R_(0,0) * R_(0,0) + R_(1,0) * R_(1,0);
        double theta, psi, phi;
        if(R11R21 == 0.0) {
            int sign = 1;
            if(R_(2,0)<0) sign = -1;
            theta = -sign*M_PI/2;
            psi = 0.0;
            sign = 1;
            if(R_(0,1)<0) sign = -1;
            phi = -sign*acos(R_(1,1));
        } else {
            int sign =1;
            if(R_(2,0)<0) sign = -1;
            theta = -sign*acos(sqrt(R_(0,0)*R_(0,0)+R_(1,0)*R_(1,0)));
            sign = 1;
            if((R_(2,1)/cos(theta))<0) sign = -1;
            psi = sign*acos(R_(2,2)/cos(theta));
            sign = 1;
            if((R_(1,0)/cos(theta))<0) sign = -1;
            phi = sign*acos(R_(0,0)/cos(theta));
        }

        d_os << setw(15) << phi/M_PI*180
            << setw(15) << theta/M_PI*180
            << setw(15) << psi/M_PI*180 << endl;

        d_os << setw(15) << box.X()
            << setw(15) << box.Y()
            << setw(15) << box.Z() << endl;
    }
}

void gio::OutCif_i::writeAtomSpecifier( const AtomSpecifier &atoms ) {

    d_os.setf( ios::fixed, ios::floatfield );
    d_os.setf( ios::unitbuf );
    d_os.precision(3);
    int res = 0;
    int count = 0;
    int resoff = 0;

    for ( auto i = 0; i != atoms.size(); ++i ) {
        int maxmol = atoms.mol(i);
        if (maxmol < 0) maxmol = atoms.sys()->numMolecules();
        count = atoms.atom(i);
        resoff = 0;
        for ( auto j = 0; j != maxmol; ++j ) {
            count += atoms.sys()->mol(j).numAtoms();
            resoff += atoms.sys()->mol(j).topology().numRes();
        }

        if ( atoms.mol(i) < 0 ) res = atoms.atom(i) / atoms.sys()->sol(0).topology().numAtoms();
        else d_os << setw(4) << atoms.sys()->mol(atoms.mol(i)).topology().resName(res).substr(0, 4).c_str();

        d_os << "_atom";
        d_os.setf( ios::right, ios::adjustfield );
        d_os << setw(7) << count + 1;
        d_os.setf(ios::left, ios::adjustfield);
        d_os << "  " << setw(4) << atoms.name(i).substr(0, 3).c_str();
        if (atoms.mol(i) < 0) d_os << setw(4) << "SOLV";
        else d_os << setw(4) << atoms.sys()->mol(atoms.mol(i)).topology().resName(res).substr(0, 4).c_str();
        d_os.setf(ios::right, ios::adjustfield);

        int resn = res + 1;
        resn += resoff;
        d_os << setw(4) << resn << "    "
             << setw(8) << atoms.pos(i)[0]*d_factor
             << setw(8) << atoms.pos(i)[1]*d_factor
             << setw(8) << atoms.pos(i)[2]*d_factor;
    }
    d_os << "_end" << endl;
}

