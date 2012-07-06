/**
 * @file temperature.cc
 * 
 * Implementation of Temperature
 */

#include <cassert>
#include <iostream>
#include <sstream>
#include <cstdio>
#include <string>
#include <set>
#include "../utils/AtomSpecifier.h"
#include "../gcore/System.h"
#include "../gmath/Vec.h"
#include "../gmath/Physics.h"
#include "../utils/Temperature.h"

utils::Temperature::Temperature (const AtomSpecifier &as, double dof)
    : m_as(as), dof(dof)
{
    //std::cout << "# Degree of freedom: = " << dof 
    //          << " Boltzmann = " << gmath::physConst.get_boltzmann() <<std::endl;
}

double
utils::Temperature::temperature(const gcore::System &sys){
    double e_kin = 0.0;
    int num_atoms = m_as.size();
    for (int i = 0; i < num_atoms; i++){
        double mass = m_as.mass(i);
        gmath::Vec vel = m_as.vel(i);
        e_kin += mass * vel.abs2();
    }
    e_kin /= (dof * gmath::physConst.get_boltzmann());
    //e_kin *= 0.5;
    return e_kin;
}
