// gmath_Physics.cc
#include <cmath>
#include <string>
#include <iostream>

#include "Physics.h"

using namespace std;

namespace gmath {

  PhysConst::PhysConst() {

    // some mathematical constants and factors
    pi = 4 * atan(1.0);
    kilo = 1000;
    nano = 1e-9;
    pico = 1e-12;
    radian2degree = 180.0 / pi;
    degree2radian = pi / 180.0;

    // other physical constants than the most basic ones (see below)
    elementary_charge = 1.60217733e-19; // unit: C
    atomic_mass_unit = 1.6605402e-27;   // unit: kg
    avogadro = 6.0221367e23;            // unit: 1/mol
    euler = 2.718281828459;

    // The most basic physical constants, which should be read from the topology
    // file. Here just the values from literature ([1], see Physics.h) to
    // inertialize them in case there is no topology (or an old toplogy format).
    // The corresponding bool variables indicat if the variable was initialized
    // from literature (false) or replaced by the value from the topology
    // file (true)
    four_pi_eps_i = 138.935484611;      // unit: (kJ nm)/(mol e^2)
    hbar = 6.35078077e-2;               // unit: (kJ ps)/mol
    speed_of_light = 2.99792458e5;      // unit: nm/ps
    boltzmann = 1.380658e-23 * avogadro / kilo;        // unit: kJ/(K mol)

    // set all variables to inertialised from literature (bbo = false)
    four_pi_eps_i_ = false;
    hbar_ = false;
    speed_of_light_ = false;
    boltzmann_ = false;
    elementary_charge_ = false;
    atomic_mass_unit_ = false;
    avogadro_ = false;
    euler_ = false;
    eps0_ = false;
    mu0_ = false;
    h_ = false;

    // calculate all dependent variables out of the ones above
    eps0 = 1.0 / (4.0 * pi * four_pi_eps_i);
    mu0 = 1.0 / (eps0 * speed_of_light * speed_of_light);
    h = hbar * 2 * pi;

  }

  void PhysConst::calc() {
    if(four_pi_eps_i_) {
      eps0_ = true;
    }
    eps0 = 1.0 / (4.0 * pi * four_pi_eps_i);
    if(eps0_ && speed_of_light_) {
      mu0_ = true;
    }
    mu0 = 1.0 / (eps0 * speed_of_light * speed_of_light);
    if(hbar_) {
      h_ = true;
    }
    h = hbar * 2 * get_pi();
  }

  double PhysConst::get_pi() {
    return pi;
  }

  double PhysConst::get_radian2degree() {
    return radian2degree;
  }

  double PhysConst::get_degree2radian() {
    return degree2radian;
  }

  double PhysConst::get_kilo() {
    return kilo;
  }
  
  double PhysConst::get_nano() {
    return nano;
  }

  double PhysConst::get_pico() {
    return pico;
  }

  double PhysConst::get_avogadro() {
    if(!avogadro_) {
      printWarning("avogadro", avogadro, avogadro_);
    }
    return avogadro;
  }

  double PhysConst::get_elementary_charge() {
    if(!elementary_charge_) {
      printWarning("elementary_charge", elementary_charge, elementary_charge_);
    }
    return elementary_charge;
  }
  
  double PhysConst::get_atomic_mass_unit() {
    if(!atomic_mass_unit_) {
      printWarning("atomic_mass_unit", atomic_mass_unit, atomic_mass_unit_);
    }
    return atomic_mass_unit;
  }
  
  double PhysConst::get_eps0() {
    if(!eps0_) {
      printWarning("eps0", eps0, eps0_);
    }
    return eps0;
  }

  double PhysConst::get_euler() {
    if(!euler_) {
      printWarning("euler", euler, euler_);
    }
    return euler;
  }

  double PhysConst::get_mu0() {
    if(!mu0_) {
      printWarning("mu0", mu0, mu0_);
    }
    return mu0;
  }

  double PhysConst::get_h() {
    if(!h_) {
      printWarning("h", h, h_);
    }
    return h;
  }

  double PhysConst::get_speed_of_light() {
    if(!speed_of_light_) {
      printWarning("speed_of_light", speed_of_light, speed_of_light_);
    }
    return speed_of_light;
  }

  double PhysConst::get_four_pi_eps_i() {
    if(!four_pi_eps_i_) {
      printWarning("four_pi_eps_i", four_pi_eps_i, four_pi_eps_i_);
    }
    return four_pi_eps_i;
  }

  double PhysConst::get_boltzmann() {
    if(!boltzmann_) {
      printWarning("boltzmann", boltzmann, boltzmann_);
    }
    return boltzmann;
  }

  double PhysConst::get_hbar() {
    if(!hbar_) {
      printWarning("hbar", hbar, hbar_);
    }
    return hbar;
  }

  void PhysConst::set_speed_of_light(double value) {
    speed_of_light = value;
    speed_of_light_ = true;
    calc();
  }

  double PhysConst::set_four_pi_eps_i(double value) {
    four_pi_eps_i = value;
    four_pi_eps_i_ = true;
    calc();
  }

  double PhysConst::set_boltzmann(double value) {
    boltzmann = value;
    boltzmann_ = true;
    calc();
  }

  double PhysConst::set_hbar(double value) {
    hbar = value;
    hbar_ = true;
    calc();
  }

  void PhysConst::printWarning(std::string name, double & value, bool & b) {
    cerr << "WARNING: hardcoded value was used for:\n";
    cerr << "         " << name << " = " << value << endl;
    // we only want the warning once
    b = true;
  }

  // the extern declared (Physics.h) class to hold/handle the physical constants
  PhysConst physConst;

}
